module ClassicalMC

using ConfParser
using Profile
using ProfileView
using ProfileSVG
include("mcmc.jl")
include("loop_update.jl")

export profile_loop_update

# To use @code_warntype
using InteractiveUtils

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


# Read non-zero elements in the right-upper triangle part of Jij
function read_Jij(Jij_file::String,num_spins::Int64)
    Jij = Array{Tuple{SpinIndex,SpinIndex,Float64,Float64,Float64,Int64}}(undef, 0)
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

function get_param(type, conf, block, key, default_value)
    if haskey(conf, block, key)
        return parse(type, retrieve(conf, block, key))
    else
        return convert(type, default_value)
    end
end

function profile_loop_update(input_file::String)
    Profile.init(delay = 0.0001)

    conf = ConfParse(input_file)
    parse_conf!(conf)

    seed = parse(Int64, retrieve(conf, "simulation", "seed"))

    num_spins = get_param(Int64,conf,"model","num_spins",0)
    Jij_file  = retrieve(conf,"model","Jij")
    Jij       = read_Jij(Jij_file,num_spins)
 
    # for profiling single process is enough.
    temperature_file = retrieve(conf,"model","temperatures")
    temps            = read_temps(temperature_file) 
    
    # laod spin configration thermaly converged.
    spin_config_file = retrieve(conf,"model","spin_config")
    spins            = read_spin_config(spin_config_file,num_spins)
  
    # Create LoopUpdater
    max_loop_length  = parse(Int64,retrieve(conf,"loop_update","max_loop_length"))
    loop_updater     = LoopUpdater{HeisenbergSpin}(num_spins,max_loop_length)
   
    # Create SingleSpinFlipUpdater 
    model   = JModel(num_spins,Jij)
    updater = SingleSpinFlipUpdater(model)

    loop_num_trial = parse(Int64,retrieve(conf,"loop_update","num_trial"))

    # Force compilation
    multi_loop_update!(loop_updater,loop_num_trial,updater,1/temps[1],
                                max_loop_length,spins,true)

    @code_warntype multi_loop_update!(loop_updater,loop_num_trial,updater,1/temps[1],
                                max_loop_length,spins,true)

    @profile multi_loop_update!(loop_updater,loop_num_trial,updater,1/temps[1],
                                max_loop_length,spins,true)
    

    Profile.print(IOContext(stdout, :displaysize => (24, 500)))
    ProfileSVG.save("profile.svg")
    Profile.init()
end

end #module ClassicalMC

using ArgParse
using .ClassicalMC

s = ArgParseSettings()
@add_arg_table! s begin
    "input" 
         help = "input_file"
         required = true
end
args = parse_args(ARGS,s)

profile_loop_update(args["input"])
