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


function read_upward_triangles(file_name::String,num_spins::Int64)
   
    L = Int(sqrt(num_spins/3))
    utris = fill((0,0,0),L^2)

    open(file_name,"r") do fp

        @assert L^2 == parse(Int64, readline(fp)) "!match num_upward_triangles. See 2d.ini and head of utriangles.txt "

        for i in 1:L^2
            str = split(readline(fp))
            s1 = parse(Int64, str[1])
            s2 = parse(Int64, str[2])
            s3 = parse(Int64, str[3])
            utris[i] = (s1,s2,s3)
        end
    end

    return utris
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
            
            @assert isapprox(sx^2+sy^2+sz^2, 1.0)

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
    min_attemps_update_temps_dist  = get_param(Int64,  conf, "simulation", "min_attemps_update_temps_dist", 100)

    if opt_temps_dist && 5*min_attemps_update_temps_dist*ex_interval > num_therm_sweeps/2
        error("num_therm_sweeps is too small for optimizing temps_dist!")
    end


    # Read non-zero elements in the right-upper triangle part of Jij
    Jij = read_Jij(Jij_file, num_spins)

    # Create single-spin flip updater
    model = JModel(num_spins, Jij)


end

end
