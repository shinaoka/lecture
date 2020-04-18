include("mcmc.jl")

using LinearAlgebra
using StaticArrays
using Random

struct LoopUpdater{T}
    num_spins::Int
    work::Array{Int}
    spins_on_loop::Array{UInt}
    new_spins::Array{T}
end

function LoopUpdater{T}(num_spins::Int64, max_loop_length::Int64) where T
    work   = zeros(Int,num_spins)
    spins_on_loop  = zeros(Int, max_loop_length)
    new_spins  = Array{T}(undef, max_loop_length)
    return LoopUpdater{T}(num_spins,work, spins_on_loop, new_spins)
end


# rewrite based on Sec.3,B of Stefan Schnabel and David P. Landau(2012)
function find_loop(spins,
                   spins_on_loop,
                   updater::SingleSpinFlipUpdater,
                   first_spin_idx,
                   second_spin_idx,
                   max_length::Int, 
                   work::Array{Int}, verbose::Bool=false)
    #=
    All elements of work must be initialized to zero.
    =#
    
    num_spins = updater.num_spins

    @assert length(work) >= num_spins

    work[first_spin_idx]  = 1
    spins_on_loop[1]      = first_spin_idx
    work[second_spin_idx] = 2
    spins_on_loop[2]      = second_spin_idx
    loop_length::Int64    = 2
    spin_before::Int64    = first_spin_idx
    current_spin_idx::Int64 = second_spin_idx

    max_coord_num = maximum(updater.coord_num)
    candidate_spins = zeros(UInt, max_coord_num)
    
    #if verbose
       #println("colors_on_loop $(colors_on_loop)")
    #end

    success = false
    status = -1

  # for speed up we need to add sum_boundary_spins to constructor LoopUpdater and then MVectorize.
    sum_boundary_spins::HeisenbergSpin = (0.,0.,0.)

    while loop_length < max_length
        #if verbose
           #println("current_spin_idx $(current_spin_idx)")
        #end
        # Search connected spins
        n_candidate = 0
        for ins in 1:updater.nn_coord_num[current_spin_idx]
            # candidate must be either the first spin on the loop (work[ns]==1) or unvisited (work[ns]==0)
            ns::SpinIndex = updater.nn_sites[ins,current_spin_idx]
            if  work[ns] <= 1 && ns != spin_before
                n_candidate += 1
                candidate_spins[n_candidate] = ns
            end
        end
        
        """  
        temp_n_candidate = 0
        for i in candidate_spins
            if i !== 0
                temp_n_candidate += 1
            end
        end
        println("DEBUG A:",temp_n_candidate)
        #println("DEBUG B:",candidate_spins)
        """
        println("DEBUG C:",n_candidate)
        
        # we need only sum of boundary spins.
        if n_candidate !== 2
            sum_boundary_spins = sum_boundary_spins .+ spins[current_spin_idx]
            continue
        end
        #println("DEBUG A:", sum_boundary_spins)
        

        if n_candidate == 0
            status = 1
            break
        end

        # next spin index must be determined by value of inner product between one before spin.
        inner_prod = zeros(Float64,length(candidate_spins))
        
        for idx in 1:length(candidate_spins)
            """
            #println("DEBUG A:",spins[spin_before])
            #println("DEBUG B:",spins[candidate_spins[idx]])
            #println("DEBUG C:",candidate_spins)
            println("DEBUG D:",candidate_spins[idx])
            println("DEBUG E:",idx)
            """
            # candidate_spins has some 0 elements
            # so spins[candidate[idx]] is acces 0th element of spins.
            if candidate_spins[idx] == 0
                continue
            end

            inner_prod[idx] = dot(spins[spin_before],spins[candidate_spins[idx]])
        end
        
        #next_spin_idx = candidate_spins[rand(1:n_candidate)]
       
        next_spin_idx = candidate_spins[findmax(inner_prod)[2]] # findmax() returns (max element,its index)
        #println("DEBUG ?:", next_spin_idx)

        if work[next_spin_idx] == 1
            # OK, we've returned to the starting point.
            status = 0
            success = true
            break
        end

        spin_before = current_spin_idx
        current_spin_idx = next_spin_idx
        loop_length += 1
        work[current_spin_idx] = loop_length
        spins_on_loop[loop_length] = current_spin_idx
    end
    #if verbose
       #println("status $(status) $(loop_length)")
    #end

    # Reset all elements of work to 0
    for l=1:loop_length
        work[spins_on_loop[l]] = 0
    end

    if success
        return loop_length, sum_boundary_spins
    else
        # Return the 0-length array of the same for type-stability
        return 0, (0.,0.,0.)
    end
end

function reflect_spins_on_loop!(loop_length::Int64,
                                spins::Array{HeisenbergSpin},
                                new_spins_on_loop::Array{HeisenbergSpin},
                                spins_on_loop::Array{UInt},
                                updater::SingleSpinFlipUpdater,
                                sum_boundary_spins::HeisenbergSpin)
         
         #implement Equ.(9)
         @assert mod(loop_length,2) == 0 "loop_length must to be even."
         perpendicular_vec = (0.,0.,0.)
         for i in 1:Int(loop_length/2)
             temp = spins[spins_on_loop[2i]] - spins[spins_on_loop[2i-1]]
             perpendicular_vec = perpendicular_vec .+ temp
         end

         #impliment Equ.(10)
         normal_vec = normalize(cross(sum_boundary_spins,cross(sum_boundary_spins,perpendicular_vec)))

         #impliment Equ.(11)
         for i in 1:loop_length
             spin_old = spins[spins_on_loop[i]]
             new_spins_on_loop[i] = normalize(spin_old .- 2 * dot(spin_old,noemal_vec) .* normal_vec)
         end
end

function compute_dE_loop(updater::SingleSpinFlipUpdater,
                          loop_length::Int,
                          spin_idx_on_loop::Array{UInt},
                          spins::Array{HeisenbergSpin},
                          new_spins_on_loop::Array{HeisenbergSpin},
                          work::Array{Int},
                          verbose::Bool=false)
    #=
    Compute change in energy
    All elements of work must be initialized to zero.
    =#

    num_spins = updater.num_spins

    for isp_loop in 1:loop_length
        work[spin_idx_on_loop[isp_loop]] = isp_loop
    end

    dE = 0.0
    for isp_loop in 1:loop_length
        ispin = spin_idx_on_loop[isp_loop]
        si_old = spins[ispin]
        for ic in 1:updater.coord_num[ispin]
            c = updater.connection[ic, ispin]
            jspin, Jx, Jy, Jz = c

            si_old = spins[ispin]
            sj_old = spins[jspin]
            si_new = new_spins_on_loop[isp_loop]
            # If the connected site is on the loop
            if work[jspin] != 0
                sj_new = new_spins_on_loop[work[jspin]]
                dE_spin = (Jx * sj_old[1] * si_old[1] + Jy * sj_old[2] * si_old[2] + Jz * sj_old[3] * si_old[3])
                dE_spin -= (Jx * sj_new[1] * si_new[1] + Jy * sj_new[2] * si_new[2] + Jz * sj_new[3] * si_new[3])
                dE += 0.5 * dE_spin
            else
                d_si = si_new .- si_old
                dE -= (Jx * sj_old[1] * d_si[1] + Jy * sj_old[2] * d_si[2] + Jz * sj_old[3] * d_si[3])
            end
        end
    end

    for isp_loop in 1:loop_length
        work[spin_idx_on_loop[isp_loop]] = 0
    end

    return dE
end


function metropolis_method!(beta::Float64,dE::Float64,
                            spins::AbstractArray{HeisenbergSpin},
                            loop_length::Int,
                            spin_idx_on_loop::Array{UInt},
                            new_spins_on_loop::Array{HeisenbergSpin},
                            num_accept::Int64)::Float64
    
    if rand(Random.GLOBAL_RNG) < exp(-beta*dE)
        spins[spin_idx_on_loop[1:loop_length]] = new_spins_on_loop[1:loop_length]
        num_accept += 1
        
        return dE
    else
        #p = exp(-beta*dE)
        #println("update failed: $(dE) $(p) $(dE*beta)")
        return 0.0
    end
end

function multi_loop_update!(loop_updater::LoopUpdater, num_trial::Int64,
                      updater::SingleSpinFlipUpdater,beta::Float64,
                      max_length::Int,
                      spins::AbstractArray{HeisenbergSpin},
                      verbose::Bool=false)
    
    # No copy
    work = loop_updater.work
    spins_on_loop = loop_updater.spins_on_loop
    new_spins = loop_updater.new_spins

    num_spins = updater.num_spins
    max_coord_num = maximum(updater.coord_num)
    
    dE   = 0.
    num_accept = 0
    num_loop_found = 0
     
    for i=1:num_trial
        
        first_spin_idx = rand(1:num_spins) 
        candidate_second_spin_idx = zeros(UInt,max_coord_num)
        for ins in 1:updater.nn_coord_num[first_spin_idx]
            candidate_second_spin_idx[ins] = updater.nn_sites[ins,first_spin_idx]
        end
        second_spin_idx = rand(candidate_second_spin_idx)

        loop_length,sum_boundary_spins = find_loop(spins,spins_on_loop,updater,first_spin_idx,
                                                   second_spin_idx,max_length,work,verbose)

        # for check detailed balance condition satisfied,test if find_loop() could find inverse loop.
        cp_spins_on_loop = copy(spins_on_loop)
        first_spin_idx_inv  = spins_on_loop[1:loop_length][end]
        second_spin_idx_inv = spins_on_loop[1:loop_length][end-1]
        loop_length_inv,sum_boundary_spins_inv = find_loop(spins,spins_on_loop,updater,first_spin_idx_inv,
                                                           second_spin_idx_inv,max_length,work,verbose)
   
        if !all(reverse(spins_on_loop),cp_spins_on_loop)
            continue
        end
        num_loop_found += 1
             
        reflect_spins_on_loop!(loop_length,spins,new_spins,spins_on_loop,updater,sum_boundary_spins)
        dE_loop = compute_dE_loop(updater,loop_length,spins_on_loop,spins,new_spins,work,verbose)
        dE += metropolis_method!(beta,dE_loop,spins,loop_length,spins_on_loop,new_spins,num_accept)

        # in metropolis_method!(),num_accept += 1 when update is accepted.
        """
        if dE_tmp != 0
            num_accept += 1
        end
        """
    end

    #if verbose
       #println("multi_loop: $(1/beta) $(timings)")
    #end

    return dE, num_loop_found/num_trial, num_accept/num_trial
end

