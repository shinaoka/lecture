using LinearAlgebra
using StaticArrays
using CPUTime
using Random

const SpinIndex = UInt32

const IsingSpin = Int8
const HeisenbergSpin = Tuple{Float64,Float64,Float64}

export JModel, SingleSpinFlipUpdater

struct JModel
    # List of non-zero entries of Jij
    num_spins::Int
    Jij::Vector{Tuple{SpinIndex,SpinIndex,Float64,Float64,Float64,Int64}}
end

function compute_energy(model::JModel, spins::AbstractVector{IsingSpin})
    return -sum([intr[3] * spins[intr[1]] * spins[intr[2]] for intr in model.Jij])
end

function compute_energy(model::JModel, spins::AbstractVector{HeisenbergSpin})

    energy = 0.0

    """
    for intr in model.Jij
        for j in 1:3
            energy += intr[j+2] * spins[intr[1]][j] * spins[intr[2]][j]
        end
    end
    """
    for intr in model.Jij
        spin1 = spins[intr[1]]
        spin2 = spins[intr[2]]
        for j in 1:3
            energy += intr[j+2] * spin1[j] * spin2[j]
        end
    end
    
    return -energy
end

function propose_unifo()
   work = MVector(0.0, 0.0)
   i = 0
   s = 0.0
   while true
       for i=1:2
           work[i] = 2 * rand() - 1
       end
       s  = work[1]^2 + work[2]^2
       if s < 1
           break
       end
   end
   sqrt_tmp = sqrt(1-s)
   return 2*work[1]*sqrt_tmp, 2*work[2]*sqrt_tmp, 1-(2*s)

end

struct SingleSpinFlipUpdater
    num_spins::Int
    coord_num::Vector{UInt16}
    connection::Matrix{Tuple{SpinIndex,Float64,Float64,Float64,Int64}}

    connected_sites::Matrix{SpinIndex}
    is_NN::Matrix{Bool}

    nn_coord_num::Vector{UInt16}
    nn_sites::Matrix{SpinIndex}
end

function SingleSpinFlipUpdater(model::JModel)
    num_spins = model.num_spins
    Jij = model.Jij

    coord_num = zeros(Int, num_spins)

    # Figure out which spins each spin is connected to
    connection_tmp = [Set{Tuple{SpinIndex,Float64,Float64,Float64,Int64}}() for _ in 1:num_spins]
    nn_connection_tmp = [Set{Tuple{SpinIndex,Float64,Float64,Float64}}() for _ in 1:num_spins]
    num_Jij = size(Jij, 1)
    for i_pair = 1:num_Jij
        i, j = Jij[i_pair][1:2]
        is_nn = Jij[i_pair][6] != 0
        push!(connection_tmp[i], (j, Jij[i_pair][3], Jij[i_pair][4], Jij[i_pair][5],Jij[i_pair][6]))
        push!(connection_tmp[j], (i, Jij[i_pair][3], Jij[i_pair][4], Jij[i_pair][5],Jij[i_pair][6]))
        if is_nn
            push!(nn_connection_tmp[i], (j, Jij[i_pair][3], Jij[i_pair][4], Jij[i_pair][5]))
            push!(nn_connection_tmp[j], (i, Jij[i_pair][3], Jij[i_pair][4], Jij[i_pair][5]))
        end
    end
    max_coord_num = maximum([length(connection_tmp[ispin]) for ispin in 1:num_spins])

    connection = Matrix{Tuple{SpinIndex,Float64,Float64,Float64,Int64}}(undef, max_coord_num, num_spins)
    coord_num = Vector{UInt16}(undef, num_spins)
    connected_sites = Matrix{SpinIndex}(undef, (max_coord_num, num_spins))
    is_NN = Matrix{Bool}(undef, (max_coord_num, num_spins))
    for ispin = 1:num_spins
        coord_num[ispin] = length(connection_tmp[ispin])
        connection[1:coord_num[ispin], ispin] = collect(connection_tmp[ispin])
        for (ic, val) in enumerate(connection_tmp[ispin])
            connected_sites[ic, ispin] = val[1]
            is_NN[ic, ispin] = (val[5] != 0)
        end
    end

    nn_coord_num = zeros(UInt16, num_spins)
    max_nn_coord_num = maximum([length(nn_connection_tmp[ispin]) for ispin in 1:num_spins])
    nn_sites = Matrix{SpinIndex}(undef, (max_nn_coord_num, num_spins))
    for ispin = 1:num_spins
        nn_coord_num[ispin] = length(nn_connection_tmp[ispin])
        for (ic, val) in enumerate(nn_connection_tmp[ispin])
            nn_sites[ic, ispin] = val[1]
        end
    end

    return SingleSpinFlipUpdater(num_spins, coord_num, connection, connected_sites, is_NN,
        nn_coord_num, nn_sites)
end


function one_sweep(updater::SingleSpinFlipUpdater, beta::Float64, model::JModel, spins::AbstractVector{IsingSpin})
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

function one_sweep(updater::SingleSpinFlipUpdater, beta::Float64, model::JModel, spins::AbstractVector{HeisenbergSpin})
    dE::Float64 = 0
    num_acc = 0
    for ispin in 1:model.num_spins

        # Compute effective field from the rest of spins
        eff_h::HeisenbergSpin = (0.0, 0.0, 0.0)
        for ic in 1:updater.coord_num[ispin]
            c = updater.connection[ic, ispin]
            eff_h = eff_h .+ (c[2]*spins[c[1]][1], c[3]*spins[c[1]][2], c[4]*spins[c[1]][3])
        end


        # Propose a new spin state
        si_new = propose_unifo()

         
        # Flip spin
        dE_prop = -dot(si_new .- spins[ispin], eff_h)
        if rand() < exp(-beta*dE_prop)
            spins[ispin] = si_new
            dE += dE_prop
            num_acc += 1
        end

        
        # Over relaxation.
        spins[ispin] = (2*dot(eff_h, spins[ispin])/(norm(eff_h)^2)) .* eff_h .- spins[ispin]
                
    end

    return dE, num_acc/model.num_spins
end



function gaussian_move(updater::SingleSpinFlipUpdater, beta::Float64, model::JModel, spins::AbstractVector{HeisenbergSpin}, xy::Bool=false)
    dE::Float64 = 0
    sigma_g     = sqrt(beta^-1)
    num_acc     = 0
    coeff_z     = xy ? 0.0 : 1.0
    for ispin in 1:model.num_spins
        # Compute effective field from the rest of spins
        eff_h::HeisenbergSpin = (0.0, 0.0, 0.0)
        for ic in 1:updater.coord_num[ispin]
            c = updater.connection[ic, ispin]
            eff_h = eff_h .+ (c[2]*spins[c[1]][1], c[3]*spins[c[1]][2], c[4]*spins[c[1]][3])
        end

        # Propose a new spin direction : Gaussian trial move 
        si_new = spins[ispin] .+ sigma_g .* (randn(Random.GLOBAL_RNG), randn(Random.GLOBAL_RNG), coeff_z*randn(Random.GLOBAL_RNG))
        si_new = si_new ./ norm(si_new)
 
        # Flip spin
        dE_prop = -dot(si_new .- spins[ispin], eff_h)
        if rand(Random.GLOBAL_RNG) < exp(-beta*dE_prop)
            spins[ispin] = si_new
            dE += dE_prop
            num_acc += 1
        end
        
        # Over relaxation.
        spins[ispin] = (2*dot(eff_h, spins[ispin])/(norm(eff_h)^2)) .* eff_h .- spins[ispin]
    end

    return dE, num_acc/model.num_spins
end

function potts_update(updater::SingleSpinFlipUpdater, beta::Float64, model::JModel, spins::AbstractVector{HeisenbergSpin})
    dE::Float64 = 0
    num_acc     = 0
    potts_dof   = [0.,2pi/3,4pi/3]
    for ispin in 1:model.num_spins
        # Compute effective field from the rest of spins
        eff_h::HeisenbergSpin = (0.0, 0.0, 0.0)
        for ic in 1:updater.coord_num[ispin]
            c = updater.connection[ic, ispin]
            eff_h = eff_h .+ (c[2]*spins[c[1]][1], c[3]*spins[c[1]][2], c[4]*spins[c[1]][3])
        end

        # Propose a new spin direction :
        theta  = rand(potts_dof)
        si_new = (cos(theta),sin(theta),0.)

        # Flip spin
        dE_prop = -dot(si_new .- spins[ispin], eff_h)
        if rand(Random.GLOBAL_RNG) < exp(-beta*dE_prop)
            spins[ispin] = si_new
            dE += dE_prop
            num_acc += 1
        end

        # Over relaxation.
        spins[ispin] = (2*dot(eff_h, spins[ispin])/(norm(eff_h)^2)) .* eff_h .- spins[ispin]
    end

    return dE, num_acc/model.num_spins
end
