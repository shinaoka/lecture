module ClassicalMC

export solve,get_param,read_temps,read_Jij,read_spin_config

include("mcmc.jl")
include("accumulator.jl")
include("replica_exchange.jl")
include("loop_update.jl")
include("measure_mc.jl")

using Random
using ConfParser
using ArgParse
using MPI
using Test
using CPUTime
using HDF5
using Profile


# Read a list of temperatures
function read_temps(temperature_file::String)
    temps = Vector{Float64}(undef, 0)
    num_temps = 0
    open(temperature_file) do file
        num_temps = parse(Int64, readline(file))
        for l in 1:num_temps
            temp = parse(Float64, readline(file))
            push!(temps, temp)
        end
    end

    # Check if temperatures are in ascending order or descending order
    if !all(temps[1:num_temps-1] .< temps[2:num_temps]) && !all(temps[1:num_temps-1] .> temps[2:num_temps])
        error("Temperatures must be given either in ascending order or in descending order!")
    end

    return temps
end

# Distribute temperature over MPI processes
function distribute_temps(rank, num_temps, num_proc)
    num_temps_local = fill(trunc(Int, num_temps/num_proc), (num_proc,))
    left_over = mod(num_temps, num_proc)
    num_temps_local[1:left_over] .+= 1
    @test sum(num_temps_local) == num_temps
    start_idx = sum(num_temps_local[1:rank]) + 1
    end_idx = start_idx + num_temps_local[rank+1] - 1

    return start_idx, end_idx
end

# Read non-zero elements in the right-upper triangle part of Jij
function read_Jij(Jij_file::String,num_spins::Int64)
    Jij = Vector{Tuple{SpinIndex,SpinIndex,Float64,Float64,Float64,Int64}}(undef, 0)
    open(Jij_file, "r" ) do fp
        
        @assert num_spins == parse(Int64, readline(fp)) "!match num_spins. See 2d.ini and head of Jij.txt"

        num_Jij_elems = parse(Int64, readline(fp))
        for i in 1:num_Jij_elems
            str     = split(readline(fp))
            i       = parse(SpinIndex, str[1])
            j       = parse(SpinIndex, str[2])
            val_x   = parse(Float64,   str[3])
            val_y   = parse(Float64,   str[4])
            val_z   = parse(Float64,   str[5])
            flag_nn = parse(Int64,     str[6])
            if i >= j
                error("Only right-upper triangle part must be given.")
            end
            if i > num_spins || i < 0  || j > num_spins || j < 0
                error("i or j is out of the range [1, num_spins].")
            end
            push!(Jij, (i, j, val_x, val_y, val_z,flag_nn))
        end
    end
    return Jij
end

function read_spin_config(file_name::String,num_spins::Int64)

    spins = fill((0.,0.,0.),num_spins)

    open(file_name,"r") do fp

        @assert num_spins == parse(Int64, readline(fp)) "!match num_spins. See 2d.ini and head of spin_config.txt "

        for i in 1:num_spins
            str = split(readline(fp))
            sx = parse(Float64, str[1])
            sy = parse(Float64, str[2])
            sz = parse(Float64, str[3])
            spins[i] = (sx,sy,sz)
        end
    end

    return spins
end

function write_spin_config(file_name::String,spins)
    
    num_spins = length(spins)
    open(file_name,"w") do fp
        println(fp,num_spins)
        for i in 1:num_spins 
            sx,sy,sz = spins[i]
            println(fp,sx," ",sy," ",sz)
        end
    end

end
         

function compute_magnetization(acc::Accumulator,num_spins::Int64,spins::Vector{Vector{Tuple{Float64,Float64,Float64}}},num_temps::Int64)
    
    mx = zeros(Float64, num_temps)
    my = zeros(Float64, num_temps)
    mz = zeros(Float64, num_temps)
    m2 = zeros(Float64, num_temps)

    for i in 1:num_temps
      temp_mx = 0.0    
      temp_my = 0.0    
      temp_mz = 0.0    
      for j in 1:num_spins
        temp_mx += spins[i][j][1]
        temp_my += spins[i][j][2]
        temp_mz += spins[i][j][3]
      end
      mx[i] = temp_mx
      my[i] = temp_my
      mz[i] = temp_mz
      m2[i] = float(mx[i]^2 + my[i]^2+ mz[i]^2)
    end
    
    add!(acc, "Mx2",mx.^2)    
    add!(acc, "My2",my.^2)    
    add!(acc, "Mz2",mz.^2)
    #add!(acc, "M2", m2)
    #add!(acc, "M4", m2.^2)
   
 
end


function get_param(type, conf, block, key, default_value)
    if haskey(conf, block, key)
        return parse(type, retrieve(conf, block, key))
    else
        return convert(type, default_value)
    end
end


function solve(input_file::String, comm)
    if !isfile(input_file)
        error("$input_file does not exists!")
    end
    conf = ConfParse(input_file)
    parse_conf!(conf)

    rank = MPI.Comm_rank(comm)
    num_proc = MPI.Comm_size(comm)

    num_spins = get_param(Int64, conf, "model", "num_spins", 0)
    Jij_file = retrieve(conf, "model", "Jij")
    temperature_file = retrieve(conf, "model", "temperatures")
    is_xy = get_param(Bool, conf, "model", "xy_spins", false)
    if rank == 0 && is_xy
        println("Using XY spins")
    end

    num_sweeps       = parse(Int64, retrieve(conf, "simulation", "num_sweeps"))
    num_therm_sweeps = parse(Int64, retrieve(conf, "simulation", "num_therm_sweeps"))
    meas_interval    = parse(Int64, retrieve(conf, "simulation", "meas_interval"))
    ex_interval      = parse(Int64, retrieve(conf, "simulation", "ex_interval"))
    seed             = parse(Int64, retrieve(conf, "simulation", "seed"))
    opt_temps_dist   = get_param(Bool,       conf, "simulation", "opt_temps_dist", true)

    # For loop updates
    loop_num_trial  = parse(Int64, retrieve(conf, "loop_update", "num_trial"))
    max_loop_length  = parse(Int64, retrieve(conf, "loop_update", "max_loop_length"))
    loop_interval  = get_param(Int64, conf, "loop_update", "interval", 10)

    # Read a list of temperatures
    temps = read_temps(temperature_file)
    num_temps = length(temps)
    if num_temps < num_proc
        error("Number of processes > num_temps")
    end
    if rank == 0
       println("num of temperatures = ", num_temps)
    end

    # Decide which temperatures are computed on this process
    start_idx, end_idx = distribute_temps(rank, num_temps, num_proc)
    num_temps_local = end_idx - start_idx + 1

    # Read non-zero elements in the right-upper triangle part of Jij
    Jij = read_Jij(Jij_file, num_spins)

    # Create single-spin flip updater
    model = JModel(num_spins, Jij)
    updater = SingleSpinFlipUpdater(model)

    # Init random number generator
    Random.seed!(seed + rank)

    # Create accumulator
    acc = Accumulator(num_temps_local)

    # Create accumulator for collecting stat for every process
    acc_proc = Accumulator(1)


    # Init spins
    spins_local = [fill((1.,0.,0.),num_spins) for i in 1:num_temps_local]
    
    # Optional init spin configuration.
    is_read_spin_config     = get_param(Bool, conf, "simulation", "read_spin_config", false)
    if is_read_spin_config
        spin_config_file = retrieve(conf, "model", "spin_config")
        spins_local = fill(read_spin_config(spin_config_file,num_spins),num_temps_local)
    end

    energy_local = [compute_energy(model, spins_local[it]) for it in 1:num_temps_local]

    # Replica exchange
    rex = ReplicaExchange{HeisenbergSpin}(temps, start_idx, end_idx, num_spins)
    temps = 0

    # Perform MC
    last_output_time = time_ns()
    if rank == 0
        println("Starting simulation...")
    end

  
    # Create LoopUpdater 
    loop_updater = LoopUpdater{HeisenbergSpin}(num_spins, max_loop_length)

    # For measuring acceptance rates
    single_spin_flip_acc = zeros(Float64, num_temps_local)


    for sweep in 1:num_sweeps
        # Output roughtly every 10 sececonds
        if rank == 0 && time_ns() - last_output_time > 1e+10
            println("Done $sweep sweeps")
            last_output_time = time_ns()
            flush(stdout)
        end

        elpsCPUtime = []
 
        # Single spin flips
        ts_start = CPUtime_us()
        
        for it in 1:num_temps_local
            dE, acc_rate = gaussian_move(updater, 1/rex.temps[it+start_idx-1], model, spins_local[it], is_xy)
            energy_local[it] += dE
            single_spin_flip_acc[it] = acc_rate
        end
        #println("one_sweep", " ",ts_end - ts_start)
        push!(elpsCPUtime, CPUtime_us() - ts_start)

        # Check if energy is correct
        if mod(sweep, 100) == 0
            energy_local = [compute_energy(model, spins_local[it]) for it in 1:num_temps_local]
            for it in 1:num_temps_local
                @assert abs(energy_local[it] - compute_energy(model, spins_local[it])) < 1e-5
            end
        end
         
        # Replica exchange
        ts_start = CPUtime_us()
        if mod(sweep, ex_interval) == 0
            perform!(rex, spins_local, energy_local, comm)
        end
        if opt_temps_dist
            if sweep <= Int(num_therm_sweeps/2) && rex.num_attemps >= 100
                update_temps_dist!(rex,comm)
            end
        end

        push!(elpsCPUtime, CPUtime_us() - ts_start)

        # Loop update
        ts_start = CPUtime_us()
        loop_found_rate = zeros(Float64, num_temps_local)
        loop_acc_rate = zeros(Float64, num_temps_local)
        if loop_num_trial > 0 && mod(sweep, loop_interval) == 0
            for it in 1:num_temps_local
                dE, loop_found_rate[it], loop_acc_rate[it] = multi_loop_update!(loop_updater, loop_num_trial,
                    updater,
                    1/rex.temps[it+start_idx-1],
                    max_loop_length, spins_local[it], rank==0)
                energy_local[it] += dE
            end
        end
        push!(elpsCPUtime, CPUtime_us() - ts_start)

        # Measurement
        ts_start = CPUtime_us()
        if sweep >= num_therm_sweeps && mod(sweep, meas_interval) == 0
            add!(acc, "E", energy_local)
            add!(acc, "E2", energy_local.^2)
            add!(acc, "single_spin_flip_acc", single_spin_flip_acc)
            
            add!(acc, "loop_found_rate" , loop_found_rate)
            add!(acc, "loop_accept_rate", loop_acc_rate)
  
            # order parameters
            T2_op = zeros(Float64,num_temps_local)
            for it in 1:num_temps_local
                T2_op[it] = compute_T2_op(spins_local[it],num_spins)
            end
            add!(acc, "T2_op", T2_op)
  
        end
        push!(elpsCPUtime, CPUtime_us() - ts_start)

        if sweep > num_therm_sweeps
            add!(acc_proc, "CPUtime", [Vector{Float64}(elpsCPUtime)])
        end
    end        

    #for it in 1:num_temps_local
        #dE,num_accept = multi_loop_update!(loop_num_trial,loop_num_reference_sites,updater,1/rex.temps[it+start_idx-1],triangles,max_loop_length,spins_local[it],rank==0)
    #end
    

    # Output results
    E = mean_gather(acc, "E", comm)
    E2 = mean_gather(acc, "E2", comm)
    single_spin_flip_acc = mean_gather(acc, "single_spin_flip_acc", comm)
    loop_found_rate = mean_gather(acc,"loop_found_rate", comm)
    loop_accept_rate = mean_gather(acc,"loop_accept_rate", comm)
    CPUtime = mean_gather_array(acc_proc, "CPUtime", comm)
    T2_op = mean_gather(acc, "T2_op", comm)
    flush(stdout)
    MPI.Barrier(comm)
    if rank == 0
        println()
        println("<E> <E^2> <C>")
        for i in 1:num_temps
            println(rex.temps[i], "  ", E[i], " ", E2[i], " ", ((E2[i]  - E[i]^2) / (rex.temps[i]^2)) / num_spins)
        end
        println()
        
        println("single_spin_flip_acc: ", single_spin_flip_acc)
      
        println("Acceptant rate of loop update: ")
        for i in 1:num_temps
            println(rex.temps[i], " ", loop_found_rate[i], " ", loop_accept_rate[i])
        end
        
        temp_idx = rand(1:num_spins)
       

        println("<CPUtime> ")
        for (i, t) in enumerate(CPUtime)
            println(" rank=", i-1, " : $t")
        end
   
        #write_spin_config("spin_config.txt",spins_local[1])
       
        # overwrite initial temperature distribution.        
        open("temperatures.txt","w") do fp
             println(fp,num_temps)
             for i in 1:num_temps
                 println(fp,rex.temps[i])
             end
        end

        println("T2_op")
        for i in 1:num_temps
            println(rex.temps[i]," ",T2_op[i]) 
        end


    end
    flush(stdout)
    MPI.Barrier(comm)

    # Stat of Replica Exchange MC
    print_stat(rex, comm)
end

end
