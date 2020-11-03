using LinearAlgebra
using StaticArrays
using CPUTime
using Random

const SpinIndex = UInt32

const IsingSpin = Int8
const HeisenbergSpin = Tuple{Float64,Float64,Float64}
const UniqueJijIdx = UInt8

export JModel, SingleSpinFlipUpdater

struct UniqueJij
    Jxyz::SVector{3,Float64}
    flag_nn::UInt8
end

struct JModel
    # List of non-zero entries of Jij
    num_spins::Int
    unique_Jij::Vector{UniqueJij}
    Jij::Vector{Tuple{SpinIndex,SpinIndex,UInt8}}
end


# Read non-zero elements in the right-upper triangle part of Jij
function JModel(Jij_file::String, num_spins::Int64)
    unique_Jij = UniqueJij[]
    Jij = Tuple{SpinIndex,SpinIndex,UInt8}[]
    open(Jij_file, "r" ) do fp
        @assert num_spins == parse(Int64, readline(fp)) "!match num_spins. See 2d.ini and head of Jij.txt"

        num_unique_Jij = parse(Int, readline(fp))
        for i in 1:num_unique_Jij
            str     = split(readline(fp))
            i_read  = parse(Float64,  str[1])
            val_x   = parse(Float64,  str[2])
            val_y   = parse(Float64,  str[3])
            val_z   = parse(Float64,  str[4])
            flag_nn = parse(UInt8,    str[5])
            if i_read != i
                error("Invalid input for unique Jij!")
            end
            Jxy = SVector{3,Float64}(val_x, val_y, val_z)
            push!(unique_Jij, UniqueJij(Jxy, flag_nn))
        end

        num_Jij_elems = parse(Int64, readline(fp))
        resize!(Jij, num_Jij_elems)
        for ielem in 1:num_Jij_elems
            str      = split(readline(fp))
            i        = parse(SpinIndex, str[1])
            j        = parse(SpinIndex, str[2])
            uJij_idx = parse(UInt8,     str[3])
            if i >= j
                error("Only right-upper triangle part must be given.")
            end
            if i > num_spins || i <= 0  || j > num_spins || j <= 0
                error("i or j is out of the range [1, num_spins].")
            end
            if uJij_idx > num_unique_Jij
                error("Invalid index for unique Jij!")
            end
            Jij[ielem] = (i, j, uJij_idx)
        end
    end

    return JModel(num_spins, unique_Jij, Jij)
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
        uj_idx = intr[3]
        for j in 1:3
            energy += model.unique_Jij[uj_idx].Jxyz[j] * spin1[j] * spin2[j]
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
    unique_Jij::Vector{UniqueJij}
    coord_num::Vector{UInt16}
    connection::Matrix{Tuple{SpinIndex,UInt8}}
    nn_coord_num::Vector{UInt16}
    nn_sites::Matrix{SpinIndex}
end

function get_J(updater, ispin, icoord)
    jspin, uidx = updater.connection[icoord, ispin]
    return jspin,
     updater.unique_Jij[uidx].Jxyz[1],
     updater.unique_Jij[uidx].Jxyz[2],
     updater.unique_Jij[uidx].Jxyz[3]
end

function SingleSpinFlipUpdater(model::JModel)
    num_spins = model.num_spins
    Jij = model.Jij

    coord_num = zeros(Int, num_spins)

    # Compute coord_num
    coord_num = zeros(UInt16, num_spins)
    nn_coord_num = zeros(UInt16, num_spins)
    for i_pair in eachindex(model.Jij)
        Jij_pair = Jij[i_pair]
        i, j = Jij_pair[1:2]
        coord_num[i] += 1
        coord_num[j] += 1
        if model.unique_Jij[Jij_pair[3]].flag_nn == 1
            nn_coord_num[i] += 1
            nn_coord_num[j] += 1
        end
    end
    max_coord_num = maximum(coord_num)
    max_nn_coord_num = maximum(nn_coord_num)

    # (idx of connected site, index of unique Jij, is_nn)
    connection = Matrix{Tuple{SpinIndex,UInt8}}(undef, max_coord_num, num_spins)
    coord_idx = ones(UInt16, num_spins)
    nn_coord_idx = ones(UInt16, num_spins)
    nn_sites = fill(typemax(SpinIndex), (max_nn_coord_num, num_spins))
    for i_pair in eachindex(model.Jij)
        Jij_pair = Jij[i_pair]
        i, j = Jij_pair[1:2]
        connection[coord_idx[i], i] = (j, Jij_pair[3])
        connection[coord_idx[j], j] = (i, Jij_pair[3])
        coord_idx[i] += 1
        coord_idx[j] += 1
        if model.unique_Jij[Jij_pair[3]].flag_nn == 1
            nn_sites[nn_coord_idx[i], i] = j
            nn_sites[nn_coord_idx[j], j] = i
            nn_coord_idx[i] += 1
            nn_coord_idx[j] += 1
        end
    end

    return SingleSpinFlipUpdater(num_spins, model.unique_Jij, coord_num, connection, nn_coord_num, nn_sites)
end


function one_sweep(updater::SingleSpinFlipUpdater, beta::Float64, model::JModel, spins::AbstractVector{IsingSpin})
    dE::Float64 = 0
    for ispin in 1:model.num_spins
        si_old = spins[ispin]
        # Compute effective field from the rest of spins
        eff_h::Float64 = 0.0
        for ic in 1:updater.coord_num[ispin]
            jsite, Jx, Jy, Jz = get_J(updater, ispin, ic)
            eff_h += Jz * spins[jsite]
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

function effective_field(spins::AbstractVector{HeisenbergSpin}, updater, ispin)
    eff_h::HeisenbergSpin = (0.0, 0.0, 0.0)
    for ic in 1:updater.coord_num[ispin]
        c = updater.connection[ic, ispin]
        jsite, Jx, Jy, Jz = get_J(updater, ispin, ic)
        jspin = spins[jsite]
        eff_h = eff_h .+ (Jx*jspin[1], Jy*jspin[2], Jz*jspin[3])
    end
    eff_h
end

function one_sweep(updater::SingleSpinFlipUpdater, beta::Float64, model::JModel, spins::AbstractVector{HeisenbergSpin})
    dE::Float64 = 0
    num_acc = 0
    for ispin in 1:model.num_spins
        # Compute effective field from the rest of spins
        eff_h = effective_field(spins, updater, ispin)

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
        eff_h = effective_field(spins, updater, ispin)

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


