using Random

include("mcmc.jl")

function solve()
    num_spins = 100
    J = 1.0
    beta = 100.0
    num_sweeps = 10

    # Create model with nearest neighbor J
    Jij = Array{Tuple{SpinIndex,SpinIndex,Float64}}(undef, num_spins-1)
    for ispin = 1:(num_spins-1)
        Jij[ispin] = (ispin, ispin+1, J)
    end
    model = IsingModel(num_spins, Jij)

    # Create single-spin flip updater
    updater = SingleSpinFlipUpdater(model)
    
    # Init random number generator
    seed = 1234
    Random.seed!(seed)
    
    # Perform num_sweeps sweeps
    spins = ones(Int, num_spins)
    energy = compute_energy(model, spins)
    println(energy)
    for sweep in 1:num_sweeps
        energy += one_sweep(updater, beta, model, spins)
    end
    println(energy, "?=", compute_energy(model, spins))
end

solve()
