using ConfParser
using Profile
using ArgParse
include("mcmc.jl")
include("loop_update.jl")
include("classical_mc.jl")
using .ClassicalMC


function profile_loop_update(input_file::String)
    Profile.init()

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
    @profile multi_loop_update!(loop_updater,loop_num_trial,updater,1/temps[1],
                                max_loop_length,spins,true)

    Profile.print()
    Profile.init()
end

s = ArgParseSettings()
@add_arg_table s begin
    "input" 
         help = "input_file"
         required = true
end
args = parse_args(ARGS,s)

profile_loop_update(args["input"])
