using MPI

mutable struct ReplicaExchange{T}
    # All temperatures
    temps::Vector{Float64}
    num_attemps::UInt64
    num_accepted::Vector{UInt64}
    start_idx::UInt64
    end_idx::UInt64
    #recv_buffer::Vector{T}
    recv_buffer::Vector{Float64}
    send_buffer::Vector{Float64}
    BufferType::Type
end

function ReplicaExchange{T}(temps_init::Vector{Float64}, start_idx, end_idx, num_spins) where T
    num_temps = length(temps_init)
    num_temps_local = end_idx - start_idx + 1
    if !all(temps_init[1:num_temps-1] .< temps_init[2:num_temps]) && !all(temps_init[1:num_temps-1] .> temps_init[2:num_temps])
        error("Temperatures must be given either in ascending order or in descending order!")
    end
    if any(temps_init .== 0.0)
        error("Zero temperature is not allowed!")
    end
    ReplicaExchange{T}(copy(temps_init), 0, zeros(UInt64, num_temps_local),
        start_idx, end_idx,
        #Vector{T}(undef, num_spins),
        collect(reinterpret(Float64, Vector{T}(undef, num_spins))),
        collect(reinterpret(Float64, Vector{T}(undef, num_spins))),
        Float64
    )
end

function swap_configs(energy_local::Vector{Float64}, config_local::Vector{Vector{T}}, i, j) where T
    energy_local[i], energy_local[j] = energy_local[j], energy_local[i]
    for idx in 1:length(config_local[i])
        config_local[i][idx], config_local[j][idx] = config_local[j][idx], config_local[i][idx]
    end
end

function reset_stats!(rex::ReplicaExchange)
    rex.num_attemps = 0
    rex.num_accepted[:] .= 0
end

# Update the distribution of temperatures using an idea similar to the one described in Sec. 3 of K. Hukushima and K. Nemoto (1996)
# Note: Equation (11) looks wrong.
function update_temps_dist!(rex::ReplicaExchange, comm)
    MPI.Barrier(comm)

    # mixing parameter
    alpha = 0.4
    rank = MPI.Comm_rank(comm)
    num_accepted = MPI.Allgather(rex.num_accepted, comm)
    num_temps = length(rex.temps)
    @assert length(num_accepted) == num_temps

    # Compute acceptance rates of exchange
    # plus one is for avoinding singularity.
    # Discard the last element.
    pm = (num_accepted ./ rex.num_attemps)[1:end-1]
    pm = pm .+ 1e-8
    betas = 1 ./ rex.temps
    c = (betas[end]-betas[1])/sum(pm .* (betas[2:end]-betas[1:end-1]))
    
    # Update beta
    new_betas = similar(rex.temps)
    new_betas[1] = betas[1]
    for i in 2:num_temps
        #if rank == 0
            #println(i, " ", (betas[i] - betas[i-1]) * pm[i-1] * c, " ", pm[i-1] * c)
        #end
        new_betas[i] = new_betas[i-1] + (betas[i] - betas[i-1]) * pm[i-1] * c
    end
    @assert all(new_betas .> 0)
    @assert isapprox(new_betas[end], betas[end])

    new_betas = alpha * new_betas + (1-alpha) * betas

    # Broadcast new temps
    new_temps = 1 ./ new_betas
    new_temps[1], new_temps[end] = rex.temps[1], rex.temps[end]
    rex.temps[:] = MPI.bcast(new_temps, 0, comm)

    reset_stats!(rex)
end

function perform!(rex::ReplicaExchange, config_local::Array{Array{T,1},1}, energy_local::Array{Float64}, comm) where T
    MPI.Barrier(comm)

    rex.num_attemps += 1

    rank = MPI.Comm_rank(comm)
    num_proc = MPI.Comm_size(comm)
    start_idx = rex.start_idx
    end_idx = rex.end_idx
    offset = rex.start_idx-1
    num_temps_local = rex.end_idx - rex.start_idx + 1
    temps_local = rex.temps[start_idx:end_idx]

    # Intra process
    for it in 1:num_temps_local-1
        dbeta = 1/temps_local[it+1] - 1/temps_local[it]
        dE = energy_local[it+1] - energy_local[it]
        if exp(dbeta * dE) > rand()
            swap_configs(energy_local, config_local, it, it+1)
            rex.num_accepted[it] += 1
        end
    end

    if num_proc == 1
        return
    end

    # Interprocess exchange
    # iex = 1: 0<->1, 2<->3, 4<->5, ..
    # iex = 2: 1<->2, 3<->4, 5<->6, ..
    for iex in 1:2
        if mod(rank, 2) == 0
            target_rank = rank - (-1)^iex
        else
            target_rank = rank + (-1)^iex
        end
        @assert abs(target_rank - rank) == 1

        if target_rank < 0 || target_rank >= num_proc
            continue
        end

        my_idx = target_rank > rank ? num_temps_local : 1

        # Send config from even-number process to odd-number process
        to_be_updated = false
        if mod(rank, 2) == 0
            (ene_target, _) = MPI.recv(target_rank, 1000, comm)
            (temp_target, _) = MPI.recv(target_rank, 1001, comm)
            prob = exp((1/temp_target-1/temps_local[my_idx])*(ene_target - energy_local[my_idx]))
            to_be_updated = prob > rand()
            MPI.send(to_be_updated, target_rank, 1002, comm)
        else
            MPI.send(energy_local[my_idx], target_rank, 1000, comm)
            MPI.send(temps_local[my_idx], target_rank, 1001, comm)
            (to_be_updated, _) = MPI.recv(target_rank, 1002, comm)
        end

        if rank < target_rank
            rex.num_accepted[my_idx] += to_be_updated ? 1 : 0
        end

        if !to_be_updated
            continue
        end

        # Swap configs
        rex.send_buffer[:] = reinterpret(rex.BufferType, config_local[my_idx])
        if mod(rank, 2) == 0
            rreq = MPI.Irecv!(rex.recv_buffer, target_rank, 2000, comm)
            sreq = MPI.Isend(rex.send_buffer, target_rank, 2002, comm)
            stats = MPI.Waitall!([rreq, sreq])
            config_local[my_idx][:] = reinterpret(T, rex.recv_buffer)

            (new_ene, _) = MPI.recv(target_rank, 2001, comm)
            MPI.send(energy_local[my_idx], target_rank, 2003, comm)
            energy_local[my_idx] = new_ene
        else
            sreq = MPI.Isend(rex.send_buffer, target_rank, 2000, comm)
            rreq = MPI.Irecv!(rex.recv_buffer, target_rank, 2002, comm)
            stats = MPI.Waitall!([rreq, sreq])
            config_local[my_idx][:] = reinterpret(T, rex.recv_buffer)

            MPI.send(energy_local[my_idx], target_rank, 2001, comm)
            (new_ene, _) = MPI.recv(target_rank, 2003, comm)
            energy_local[my_idx] = new_ene
        end
        #if rank == 0
           #println("DEBUGTT $(tt2-tt1) $(tt3-tt2) $(tt4-tt3) $(tt5-tt4)")
        #end
    end
    #if rank == 0
       #println("DEBUG $(t2-t1) $(t3-t2) $(t4-t3)")
    #end
end

function print_stat(rex::ReplicaExchange, comm, outf=stdout)
    rank = MPI.Comm_rank(comm)
    num_accepted = MPI.Gather(rex.num_accepted, 0, comm)
    if MPI.Comm_rank(comm) == 0
        println(outf, "")
        println(outf, "Statistics of replica exchange Monte Carlo")
        for it in 1:length(rex.temps)-1
            println(outf, "itemp ", it, " <=> ", it+1, " acceptance_rate= ", num_accepted[it]/rex.num_attemps)
        end
    end
end
