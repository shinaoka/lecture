using Random

include("single_spin_flip.jl")

function solve()
    num_spins = 100
    J = 1.0
    beta = 100.0
    num_sweeps = 10

    Jij = Array{Tuple{SpinIndex,SpinIndex,Float64}}(undef, num_spins-1)

    # Nearest neighbor J
    for ispin = 1:(num_spins-1)
        Jij[ispin] = (ispin, ispin+1, J)
    end
    
    model = IsingModel(num_spins, Jij)
    
    # Create random number generator
    seed = 1234
    Random.seed!(seed)
    
    updater = SingleSpinFlipUpdater(model)
    
    spins = ones(Int, num_spins)
    energy = compute_energy(model, spins)
    println(energy)
    for sweep in 1:10
        energy += one_sweep(updater, beta, model, spins)
    end
    println(energy, "?=", compute_energy(model, spins))
end

solve()
