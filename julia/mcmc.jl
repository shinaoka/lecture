SpinIndex = Int

struct IsingModel
    # List of non-zero entries of Jij
    num_spins::Int
    Jij::Array{Tuple{SpinIndex,SpinIndex,Float64}}
end

function compute_energy(model::IsingModel, spins::Array{Int})
    return -sum([intr[3] * spins[intr[1]] * spins[intr[2]] for intr in model.Jij])
end

struct SingleSpinFlipUpdater
    num_spins::Int
    coord_num::Array{Int}
    connection::Array{Tuple{SpinIndex,Float64}}

    function SingleSpinFlipUpdater(model::IsingModel)
        num_spins = model.num_spins
        Jij = model.Jij

        coord_num = zeros(Int, num_spins)

        # Figure out which spins each spin is connected to
        connection_tmp = [Set{Tuple{SpinIndex,Float64}}() for _ in 1:num_spins]
        num_Jij = size(Jij, 1)
        for i_pair = 1:num_Jij
            i, j = Jij[i_pair][1:2]
            push!(connection_tmp[i], (j, Jij[i_pair][3]))
            push!(connection_tmp[j], (i, Jij[i_pair][3]))
        end
        max_coord_num = maximum([length(connection_tmp[ispin]) for ispin in 1:num_spins])

        connection = Array{Tuple{SpinIndex,Float64}}(undef, max_coord_num, num_spins)
        coord_num = Array{Int}(undef, num_spins)
        for ispin = 1:num_spins
            coord_num[ispin] = length(connection_tmp[ispin])
            connection[1:coord_num[ispin], ispin] = collect(connection_tmp[ispin])
        end

        new(num_spins, coord_num, connection)
    end
end

function one_sweep(updater::SingleSpinFlipUpdater, beta::Float64, model::IsingModel, spins::Array{Int})
    dE::Float64 = 0
    for ispin in 1:model.num_spins
        si_old = spins[ispin]

        # Compute effective field from the rest of spins
        eff_h::Float64 = 0.0
        for ic in 1:updater.coord_num[ispin]
            c = updater.connection[ic, ispin]
            eff_h += c[2] * spins[c[1]]
        end

        # Flip spin
        dE_ud = -2 * eff_h
        prob_up = 1/(1 + exp(beta * dE_ud))
        if rand() < prob_up
            spins[ispin] = 1
        else
            spins[ispin] = -1
        end

        # Compute energy change
        dE += - eff_h * (spins[ispin] - si_old)
    end
    return dE
end
