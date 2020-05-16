include("mcmc.jl")

using LinearAlgebra
using StaticArrays
using Random

struct LoopUpdater{T}
    num_spins::Int
    work::Vector{Int}
    spins_on_loop::Vector{UInt}
    new_spins::Vector{T}
end

function LoopUpdater{T}(num_spins::Int64, max_loop_length::Int64) where T
    work          = zeros(Int,num_spins)
    spins_on_loop = zeros(Int, max_loop_length)
    new_spins     = Vector{T}(undef, max_loop_length)
    return LoopUpdater{T}(num_spins,work, spins_on_loop, new_spins)
end


# rewrite based on Sec.3,B of Stefan Schnabel and David P. Landau(2012)
function find_loop(spins,
                   spins_on_loop,
                   updater::SingleSpinFlipUpdater,
                   first_spin_idx,
                   second_spin_idx,
                   max_length::Int, 
                   work::Vector{Int}, verbose::Bool=false,check_n_candidate::Bool=false)
    #=
    All elements of work must be initialized to zero.
    =#
    

    num_spins = updater.num_spins

    @assert length(work) >= num_spins

    work[first_spin_idx]    = 1
    work[second_spin_idx]   = 2
    spins_on_loop[1]        = first_spin_idx
    spins_on_loop[2]        = second_spin_idx
    loop_length::Int64      = 2
    spin_before::Int64      = first_spin_idx
    current_spin_idx::Int64 = second_spin_idx

    max_coord_num   = maximum(updater.nn_coord_num)
    candidate_spins = zeros(UInt, max_coord_num)

    success = false
    status = -1

    sum_boundary_spins::HeisenbergSpin = (0.,0.,0.)
    inner_prod = zeros(Float64, max_coord_num)
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

        if n_candidate == 0
            status = 1
            break
        end
        

        # next spin index must be determined by value of inner product between one before spin.
        for idx in 1:n_candidate
            inner_prod[idx] = dot(spins[spin_before],spins[candidate_spins[idx]])
        end

        max_inner_prod = findmax(inner_prod[1:n_candidate])[2] # findmax() returns (max element,its index)
        next_spin_idx  = candidate_spins[1:n_candidate][max_inner_prod]

        for idx in 1:n_candidate
            if candidate_spins[idx] == next_spin_idx
                continue
            end
            sum_boundary_spins = sum_boundary_spins .+ spins[candidate_spins[idx]]
        end


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
                                spins::Vector{HeisenbergSpin},
                                new_spins_on_loop::Vector{HeisenbergSpin},
                                spins_on_loop::Vector{UInt},
                                updater::SingleSpinFlipUpdater,
                                sum_boundary_spins::HeisenbergSpin)
         
    #implement Equ.(9)
    @assert mod(loop_length,2) == 0 "loop_length must to be even."
    perpendicular_vec = zeros(Float64,3)
    for i in 1:loop_length
        perpendicular_vec .+= (-1)^i * collect(spins[spins_on_loop[i]])
    end

    #implement Equ.(10)
    sum_boundary_spins = collect(sum_boundary_spins)
    normal_vec = normalize(cross(sum_boundary_spins,cross(sum_boundary_spins,perpendicular_vec)))
         
    #println("DEBUG A: ",normal_vec[1]," ",normal_vec[2]," ",normal_vec[3])
    #println("DEBUG B: ",dot(normal_vec,sum_boundary_spins))

    #implement Equ.(11)
    for i in 1:loop_length
        spin_old = collect(spins[spins_on_loop[i]])
        new_spins_on_loop[i] = Tuple(normalize(spin_old - 2 * dot(spin_old,normal_vec) * normal_vec))
    end

end


function compute_dE_loop(updater::SingleSpinFlipUpdater,
                          loop_length::Int,
                          spin_idx_on_loop::Vector{UInt},
                          spins::Vector{HeisenbergSpin},
                          new_spins_on_loop::Vector{HeisenbergSpin},
                          work::Vector{Int},
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
                dE_spin  = (Jx * sj_old[1] * si_old[1] + Jy * sj_old[2] * si_old[2] + Jz * sj_old[3] * si_old[3])
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


function multi_loop_update!(loop_updater::LoopUpdater, num_trial::Int64,
                            updater::SingleSpinFlipUpdater,beta::Float64,
                            max_length::Int,
                            spins::AbstractVector{HeisenbergSpin},
                            verbose::Bool=false)
    
    # No copy
    work = loop_updater.work
    spins_idx_on_loop = loop_updater.spins_on_loop
    new_spins_on_loop = loop_updater.new_spins

    num_spins = updater.num_spins
    max_coord_num = maximum(updater.coord_num)
    
    dE             = 0.
    num_accept     = 0
    num_loop_found = 0
 
    # temp vairable for check detailed balance condition(dbc).
    num_accept_dbc     = 0
    num_loop_found_dbc = 0
 
    for i=1:num_trial
      
        t1_s = time_ns()
        first_spin_idx = rand(1:num_spins) 
        candidate_second_spin_idx = zeros(UInt,max_coord_num)
        nn_coord_num = updater.nn_coord_num[first_spin_idx]
        for ins in 1:nn_coord_num
            candidate_second_spin_idx[ins] = updater.nn_sites[ins,first_spin_idx]
        end
        second_spin_idx = rand(candidate_second_spin_idx[1:nn_coord_num])
        
        loop_length,sum_boundary_spins = find_loop(spins,spins_idx_on_loop,updater,first_spin_idx,
                                                   second_spin_idx,max_length,work,verbose)

        if loop_length == 0 || mod(loop_length,2) !== 0
            continue
        end
        num_loop_found += 1
        num_loop_found_dbc += 1

        t1_e = time_ns()
        #if verbose
            #println("loop_find: $(1/beta) $(t1_e - t1_s)")
        #end
  
        t2_s = time_ns()
        before_flipped_spins = copy(spins[spins_idx_on_loop[1:loop_length]])
        reflect_spins_on_loop!(loop_length,spins,new_spins_on_loop,spins_idx_on_loop,updater,sum_boundary_spins)
        dE_loop = compute_dE_loop(updater,loop_length,spins_idx_on_loop,spins,new_spins_on_loop,work,verbose)

        # implement metropolis method.
        temp_r = rand(Random.GLOBAL_RNG)
        if temp_r < exp(-beta*dE_loop)
            spins[spins_idx_on_loop[1:loop_length]] = new_spins_on_loop[1:loop_length]
            num_accept += 1
            num_accept_dbc += 1
        else
            continue
        end
        
        t2_e = time_ns() 
        #if verbose
            #println("metropolis_method: $(1/beta) $(t2_e - t2_s)")
        #end

        #println("DEBUG A': ",temp_r," ",beta," ",dE_loop)
    
        t3_s = time_ns()
        # for check detailed balance condition satisfied,test if find_loop() could find inverse loop.
        cp_spins_idx_on_loop = copy(spins_idx_on_loop[1:loop_length])
        first_spin_idx_inv   = spins_idx_on_loop[loop_length]
        second_spin_idx_inv  = spins_idx_on_loop[loop_length-1]
        loop_length_inv,sum_boundary_spins_inv = find_loop(spins,spins_idx_on_loop,updater,first_spin_idx_inv,
                                                           second_spin_idx_inv,max_length,work,verbose)

        if reverse(spins_idx_on_loop[1:loop_length]) != cp_spins_idx_on_loop || loop_length !== loop_length_inv
            spins[cp_spins_idx_on_loop] = before_flipped_spins    
            num_loop_found -= 1
            num_accept     -= 1
            continue
        end
        
        t3_e = time_ns()
        #if verbose
            #println("find_loop_inv: $(1/beta) $(t3_e - t3_s)")
        #end

        dE += dE_loop
        
    end
    
    if verbose
        println("DEBUG A : $(num_loop_found_dbc/num_trial)")
        println("DEBUG A': $(num_loop_found/num_trial)")
        println("DEBUG B : $(num_accept_dbc/num_trial)")
        println("DEBUG B': $(num_accept/num_trial)")
    end
    
    #if num_accept != 0
        #println("DEBUG F: ", num_accept)
    #end

    return dE, num_loop_found/num_trial, num_accept/num_trial
end

