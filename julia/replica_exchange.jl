using MPI

struct ReplicaExchange
    #num_temps::UInt64
    temps::Array{Float64}
    num_attemps::Array{UInt64}
    num_accepted::Array{UInt64}
    start_idx::UInt64
    end_idx::UInt64
end

function ReplicaExchange(temps_init::Array{Float64}, start_idx, end_idx)
    num_temps = length(temps_init)
    num_temps_local = end_idx - start_idx + 1
    ReplicaExchange(temps_init, zeros(UInt64, num_temps_local), zeros(UInt64, num_temps_local), start_idx, end_idx)
end

function swap_temps(a, i, j)
    a[i], a[j] = a[j], a[i]
end

function perform!(rex::ReplicaExchange, config_local, energy_local, comm)
    rank = MPI.Comm_rank(comm)
    num_proc = MPI.Comm_size(comm)
    start_idx = rex.start_idx
    end_idx = rex.end_idx
    offset = rex.start_idx-1
    num_temps_local = rex.end_idx - rex.start_idx + 1
    temps_local = rex.temps[start_idx:end_idx]

    # Intra process
    for it in 1:num_temps_local-1
        rex.num_attemps[it] += 1
        dbeta = 1/temps_local[it+1] - 1/temps_local[it]
        dE = energy_local[it+1] - energy_local[it]
        if exp(dbeta * dE) > rand()
            swap_temps(energy_local, it, it+1)
            swap_temps(config_local, it, it+1)
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
        MPI.Barrier(comm)

        if mod(rank, 2) == 0
            target_rank = rank - (-1)^iex
        else
            target_rank = rank + (-1)^iex
        end

        if target_rank < 0 || target_rank >= num_proc
            continue
        end

        my_config_idx = target_rank > rank ? num_temps_local : 1

        # Send config from even-number process to odd-number process
        to_be_updated = false
        if mod(rank, 2) == 0
            (ene_target, _) = MPI.recv(target_rank, 1000, comm)
            (temp_target, _) = MPI.recv(target_rank, 1001, comm)
            prob = exp((1/temp_target-1/rex.temps[my_config_idx])*(ene_target - energy_local[my_config_idx]))
            to_be_updated = prob > rand()
            MPI.send(to_be_updated, target_rank, 1002, comm)
        else
            MPI.send(energy_local[my_config_idx], target_rank, 1000, comm)
            MPI.send(temps_local[my_config_idx], target_rank, 1001, comm)
            (to_be_updated, _) = MPI.recv(target_rank, 1002, comm)
        end

        if rank < target_rank
            rex.num_attemps[my_config_idx] += 1
            rex.num_accepted[my_config_idx] += to_be_updated ? 1 : 0
        end

        if !to_be_updated
            continue
        end

        # Swap configs
        if mod(rank, 2) == 0
            (new_config, _) = MPI.recv(target_rank, 2000, comm)
            (new_ene, _) = MPI.recv(target_rank, 2001, comm)
            MPI.send(config_local[my_config_idx], target_rank, 2002, comm)
            MPI.send(energy_local[my_config_idx], target_rank, 2003, comm)
            config_local[my_config_idx] = new_config
            energy_local[my_config_idx] = new_ene
        else
            MPI.send(config_local[my_config_idx], target_rank, 2000, comm)
            MPI.send(energy_local[my_config_idx], target_rank, 2001, comm)
            (new_config, _) = MPI.recv(target_rank, 2002, comm)
            (new_ene, _) = MPI.recv(target_rank, 2003, comm)
            config_local[my_config_idx] = new_config
            energy_local[my_config_idx] = new_ene
        end

    end
end

function print_stat(rex::ReplicaExchange)
    println("Statistics of replica exchange Monte Carlo")
    num_attemps = MPI.Gather(rex.num_attemps, 0, comm)
    num_accepted = MPI.Gather(rex.num_accepted, 0, comm)
    if MPI.Comm_rank(comm) == 0
        for it in 1:length(num_attemps)-1
            println("itemp ", it, " <=> ", it+1, " acceptance_rate= ", num_accepted[it]/num_attemps[it])
        end
    end
end
