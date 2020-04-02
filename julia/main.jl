using Random
using ConfParser
using ArgParse
using MPI
using Test
using CPUTime
using HDF5
using Profile

include("mcmc.jl")
include("accumulator.jl")
include("replica_exchange.jl")
include("loop_update.jl")

# Read a list of temperatures
function read_temps(temperature_file::String)
    temps = Array{Float64}(undef, 0)
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
function read_Jij(Jij_file::String, num_spins)
    Jij = Array{Tuple{SpinIndex,SpinIndex,Float64,Float64,Float64,Int64}}(undef, 0)
    open(Jij_file, "r" ) do fp
        num_Jij_elems = parse(Int64, readline(fp))
        for i in 1:num_Jij_elems
            str = split(readline(fp))
            i         = parse(SpinIndex, str[1])
            j         = parse(SpinIndex, str[2])
            val_x     = parse(Float64,   str[3])
            val_y     = parse(Float64,   str[4])
            val_z     = parse(Float64,   str[5])
            typeofJij = parse(Int64,     str[6])
            if i >= j
                error("Only right-upper triangle part must be given.")
            end
            if i > num_spins || i < 0  || j > num_spins || j < 0
                error("i or j is out of the range [1, num_spins].")
            end
            push!(Jij, (i, j, val_x, val_y, val_z,typeofJij))
        end
    end
    return Jij
end

function compute_magnetization(acc::Accumulator,num_spins::Int64,spins::Array{Array{Tuple{Float64,Float64,Float64},1},1},num_temps::Int64)
    
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

function latest_spin_direction(acc::Accumulator,spins::Array{Array{Tuple{Float64,Float64,Float64},1},1},num_spins::Int64,num_temps::Int64)
    
    spins_x = [[spins[i][j][1] for j in 1:num_spins] for i in 1:num_temps]
    spins_y = [[spins[i][j][2] for j in 1:num_spins] for i in 1:num_temps]
    spins_z = [[spins[i][j][3] for j in 1:num_spins] for i in 1:num_temps]
    """    
    @assert begin
        for it in 1:num_temps
            for spin in 1:num_spins
                temp = spins_x[it][spins]^2 + spins_y[it][spins]^2 + spins_z[it][spins]^2
                abs(1-temp) < 1e-2
            end
        end
    end         
    """ 
    add!(acc, "sx", spins_x) 
    add!(acc, "sy", spins_y)
    add!(acc, "sz", spins_z)
    
end

function make_kagome(num_spins)
    a1 = 2 .* (1.,0.)
    a2 = 2 .* (cos(pi/3),sin(pi/3))

    L = Int(sqrt(num_spins/3))
    unit_cell_position = fill((0.,0.),L^2)
    index = 1

    # complexity O(L^2~N)
    for (i,j) in Iterators.product(1:L,1:L)

        unit_cell_position[index] = (i-1) .* a1 .+ (j-1) .* a2 
        index += 1
    end 
    
    return unit_cell_position
end

function order_parameter(spins::Array{Array{Tuple{Float64,Float64,Float64},1},1},
                         num_spins::Int64, num_temps::Int64, q::Tuple{Float64,Float64}, unit_cell::Array{Tuple{Float64,Float64}})
    
    #unit_cell = make_kagome(num_spins) 

    M2_AF = zeros(num_temps) 

    for temp in 1:num_temps
        gp1 = spins[temp][1:3:num_spins] 
        gp2 = spins[temp][2:3:num_spins] 
        gp3 = spins[temp][3:3:num_spins] 
        

        temp_vec1 = (0.,0.,0.)
        temp_vec2 = (0.,0.,0.)
        temp_vec3 = (0.,0.,0.)

        for i in 1:Int(num_spins/3)

            phase_iqr = -im*dot(q,unit_cell[i])
   
            temp_vec1 = temp_vec1 .+ gp1[i] .* exp(phase_iqr) 
            temp_vec2 = temp_vec2 .+ gp2[i] .* exp(phase_iqr) 
            temp_vec3 = temp_vec3 .+ gp2[i] .* exp(phase_iqr) 
        end

        M2_AF[temp] = norm(temp_vec1)^2 + norm(temp_vec2)^2 + norm(temp_vec3)^2
    end
    
    return M2_AF
end

# Kronecker delta
function delta(i::Int64,j::Int64)
    if i == j
        return 1
    else
        return 0
    end

end

function octopolar_v2(spins,num_spins::Int64,num_temps::Int64)
    
    T = zeros(num_temps)
    
    for i in 1:num_temps
        temp = zeros(3^3)

        for j in 1:num_spins
 
            spin = spins[i][j]
            index = 1

            for (a,b,c) in Iterators.product(1:3,1:3,1:3)
                temp[index] +=  spin[a]*spin[b]*spin[c] - (spin[a]*delta(b,c) + spin[b]*delta(c,a) + spin[c]*delta(a,b))/5
                index += 1 
            end
        end
        T[i] = sum(temp.^2)
 
    end
    return T
end

function octopolar_orderparameter(spins::Array{Array{Tuple{Float64,Float64,Float64},1},1},num_spins::Int64,num_temps::Int64)
    
    op = zeros(num_temps)
     
    for temp in 1:num_temps
        
        for (i,j) in Iterators.product(1:num_spins,1:num_spins)

            op[temp] += dot(spins[temp][i],spins[temp][j])^3 - (3/5)*dot(spins[temp][i],spins[temp][j])
        end

    end
    
    return op
end

function solve(input_file::String, comm)
    if !isfile(input_file)
        error("$input_file does not exists!")
    end
    conf = ConfParse(input_file)
    parse_conf!(conf)

    rank = MPI.Comm_rank(comm)
    num_proc = MPI.Comm_size(comm)

    num_spins = parse(Int64, retrieve(conf, "model", "num_spins"))
    Jij_file = retrieve(conf, "model", "Jij")
    temperature_file = retrieve(conf, "model", "temperatures")
    is_xy = parse(Bool, retrieve(conf, "model", "xy_spins"))
    if rank == 0 && is_xy
        println("Using XY spins")
    end

    num_sweeps       = parse(Int64, retrieve(conf, "simulation", "num_sweeps"))
    num_therm_sweeps = parse(Int64, retrieve(conf, "simulation", "num_therm_sweeps"))
    meas_interval    = parse(Int64, retrieve(conf, "simulation", "meas_interval"))
    ex_interval      = parse(Int64, retrieve(conf, "simulation", "ex_interval"))
    seed             = parse(Int64, retrieve(conf, "simulation", "seed"))

    # For loop updates
    loop_num_trial  = parse(Int64, retrieve(conf, "loop_update", "num_trial"))
    loop_num_reference_sites  = parse(Int64, retrieve(conf, "loop_update", "num_reference_sites"))
    max_loop_length  = parse(Int64, retrieve(conf, "loop_update", "max_loop_length"))

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
    energy_local = [compute_energy(model, spins_local[it]) for it in 1:num_temps_local]

    # Replica exchange
    rex = ReplicaExchange{HeisenbergSpin}(temps, start_idx, end_idx, num_spins)
    temps = 0

    # Perform MC
    last_output_time = time_ns()
    if rank == 0
        println("Starting simulation...")
    end

    # Find all triangles for loop updates
    triangles = find_triangles(model, updater)
  
    # Create LoopUpdater 
    loop_updater = LoopUpdater(num_spins)

    # For measuring acceptance rates
    single_spin_flip_acc = zeros(Float64, num_temps_local)
    
    # For function order_parameter()
    unit_cell = make_kagome(num_spins) 

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
        if sweep <= Int(num_therm_sweeps/2) && rex.num_attemps >= 100
            update_temps_dist!(rex,comm)
        end
        push!(elpsCPUtime, CPUtime_us() - ts_start)

        # Loop update
        ts_start = CPUtime_us()
        accept_rate = zeros(Float64, num_temps_local)
        if loop_num_trial > 0
            for it in 1:num_temps_local
                dE, num_accept = multi_loop_update!(loop_updater, loop_num_trial,
                    loop_num_reference_sites,updater,
                    1/rex.temps[it+start_idx-1],
                    triangles, max_loop_length, spins_local[it], rank==0)
                energy_local[it] += dE
                accept_rate[it]   = num_accept 
            end
        end
        push!(elpsCPUtime, CPUtime_us() - ts_start)

        # Measurement
        ts_start = CPUtime_us()
        if sweep >= num_therm_sweeps && mod(sweep, meas_interval) == 0
            add!(acc, "E", energy_local)
            add!(acc, "E2", energy_local.^2)
            add!(acc, "single_spin_flip_acc", single_spin_flip_acc)
            
            add!(acc, "accept_rate", accept_rate)

            #compute_magnetization(acc, num_spins, spins_local, num_temps_local)
            add!(acc,"M2_AF",order_parameter(spins_local,num_spins,num_temps_local,(0.,0.),unit_cell))
            add!(acc,"op",octopolar_v2(spins_local,num_spins,num_temps_local))
  
        end
        push!(elpsCPUtime, CPUtime_us() - ts_start)

        if sweep > num_therm_sweeps
            add!(acc_proc, "CPUtime", [Array{Float64}(elpsCPUtime)])
        end
    end        

    #for it in 1:num_temps_local
        #dE,num_accept = multi_loop_update!(loop_num_trial,loop_num_reference_sites,updater,1/rex.temps[it+start_idx-1],triangles,max_loop_length,spins_local[it],rank==0)
    #end
    

    # Output results
    E = mean_gather(acc, "E", comm)
    E2 = mean_gather(acc, "E2", comm)
    single_spin_flip_acc = mean_gather(acc, "single_spin_flip_acc", comm)
    accept_rate = mean_gather(acc,"accept_rate", comm)
    #ss = mean_gather_array(acc, "ss", comm)
    #Mz2 = mean_gather(acc, "Mz2", comm)
    #M2 = mean_gather(acc, "M2", comm)
    #M4 = mean_gather(acc, "M4", comm)
    CPUtime = mean_gather_array(acc_proc, "CPUtime", comm)
    #M1 = mean_gather_array(acc, "M1", comm)
    #M2 = mean_gather_array(acc, "M2", comm)
    #M3 = mean_gather_array(acc, "M3", comm)
    op = mean_gather(acc, "op", comm)
    M2_AF = mean_gather(acc,"M2_AF",comm)
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
            println(rex.temps[i], " ", accept_rate[i])
        end
        
        """  
        println("octopolar: ",op/num_spins^2)
        
        open("g.dat", "w") do fp
            for i in 1:num_temps
                g = (3 - (M4[i]/(M2[i]^2))) / 2
                println(fp, rex.temps[i], " ", g)
            end
        end
        """ 
        
        println("temperatures less than 0.001: ",length(rex.temps[rex.temps .< 0.001]))

       
        println("AFTER CLOSE H5FILE")
        
        println("<CPUtime> ")
        for (i, t) in enumerate(CPUtime)
            println(" rank=", i-1, " : $t")
        end

    end
    flush(stdout)
    MPI.Barrier(comm)

    # Stat of Replica Exchange MC
    print_stat(rex)

end

s = ArgParseSettings()
@add_arg_table s begin
    "input"
        help = "input file"
        required = true
end
args = parse_args(ARGS, s)

# Initialize MPI environment
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

if rank == 0
    @time solve(args["input"], comm)
else
    solve(args["input"], comm)
end

MPI.Finalize()
