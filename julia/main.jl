using Random
using ConfParser
using ArgParse

include("mcmc.jl")
include("accumulator.jl")

# Read a list of temperatures
function read_temps(temperature_file::String)
    temps = Array{Float64}(undef, 0)
    num_temps = 0
    open(temperature_file) do file
        num_temps = parse(Int64, readline(file))
        for ln in eachline(file)
            temp = parse(Float64, ln)
            push!(temps, temp)
        end
    end
    return temps
end

# Read non-zero elements in the right-upper triangle part of Jij
function read_Jij(Jij_file::String, num_spins)
    Jij = Array{Tuple{SpinIndex,SpinIndex,Float64}}(undef, 0)
    open(Jij_file, "r" ) do fp
        num_Jij_elems = parse(Int64, readline(fp))
        for i in 1:num_Jij_elems
            str = split(readline(fp))
            i = parse(SpinIndex, str[1])
            j = parse(SpinIndex, str[2])
            val = parse(Float64, str[3])
            if i >= j
                error("Only right-upper triangle part must be given.")
            end
            if i > num_spins || i < 0  || j > num_spins || j < 0
                error("i or j is out of the range [1, num_spins].")
            end
            push!(Jij, (i, j, val))
        end
    end
    return Jij
end

function solve(input_file::String)
    if !isfile(input_file)
        error("$input_file does not exists!")
    end
    conf = ConfParse(input_file)
    parse_conf!(conf)

    num_spins = parse(Int64, retrieve(conf, "model", "num_spins"))
    Jij_file = retrieve(conf, "model", "Jij")
    temperature_file = retrieve(conf, "model", "temperatures")

    num_sweeps = parse(Int64, retrieve(conf, "simulation", "num_sweeps"))
    num_therm_sweeps = parse(Int64, retrieve(conf, "simulation", "num_therm_sweeps"))
    meas_interval = parse(Int64, retrieve(conf, "simulation", "meas_interval"))
    seed = parse(Int64, retrieve(conf, "simulation", "seed"))

    # Read a list of temperatures
    temps = read_temps(temperature_file)
    num_temps = length(temps)
    println("numb of temperatures = ", num_temps)

    # Read non-zero elements in the right-upper triangle part of Jij
    Jij = read_Jij(Jij_file, num_spins)

    # Create single-spin flip updater
    model = IsingModel(num_spins, Jij)
    updater = SingleSpinFlipUpdater(model)

    # Init random number generator
    Random.seed!(seed)

    # Create accumulator
    acc = Accumulator()

    # Perform num_sweeps sweeps
    for it in 1:num_temps
        println("T ", temps[it])
        spins = ones(Int, num_spins)
        energy = compute_energy(model, spins)
        for sweep in 1:num_sweeps
            energy += one_sweep!(updater, 1/temps[it], model, spins)
            if sweep > num_therm_sweeps && mod(sweep, meas_interval) == 0
                add!(acc, "E", energy)
                add!(acc, "E2", energy^2)
            end
        end
        println("<E> ", mean(acc, "E"))
        println("<E^2> ", mean(acc, "E2"))
    end
end


s = ArgParseSettings()
@add_arg_table s begin
    "input"
        help = "input file"
        required = true
end
args = parse_args(ARGS, s)
println(args)
solve(args["input"])
