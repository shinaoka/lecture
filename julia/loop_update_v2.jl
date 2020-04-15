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
function find_loop(spins::Array{HeisenbergSpin}, spins_on_loop, updater::SingleSpinFlipUpdater,first_spin_idx::Int, second_spin_idx::Int,max_length::Int, work::Array{Int}, verbose::Bool=false)::Int64
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
    spin_before::Int64    = -1
    current_spin_idx::Int64 = second_spin_idx

    max_coord_num = maximum(updater.coord_num)
    candidate_spins = zeros(UInt, max_coord_num)

    #if verbose
       #println("colors_on_loop $(colors_on_loop)")
    #end
    success = false
    status = -1
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
        inner_prod = zeros(Flaot64,length(candidate_spins))
        for idx in length(candidate_spins)
            inner_prod[idx] = dot(spins[current_spin_idx],spins[candidate_spins[idx]])
        end
        
        next_spin_idx = candidate_spins[findmax(inner_prod)[2]] # findmax() returns (max element,its index)
        
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
        return loop_length
    else
        # Return the 0-length array of the same for type-stability
        return 0
    end
end


function find_boundary_spins(spins_on_loop, updater::SingleSpinFlipUpdater, verbose::Bool=false)::Int64
    #=
    All elements of work must be initialized to zero.
    =#
    num_spins = updater.num_spins

    @assert length(work) >= num_spins

    max_coord_num = maximum(updater.coord_num)
    candidate_spins = zeros(UInt, max_coord_num)

    #if verbose
       #println("colors_on_loop $(colors_on_loop)")
    #end
   
    boundary_spin_idx = zeros(UInt,num_spins)
    counter = 1
    for idx in spins_on_loop
        for ins in 1:updater.nn_coord_num[idx]
            ns::SpinIndex = updater.nn_sites[ins,idx]
            
            if in(ns,spins_on_loop) || in(ns,boundary_spins_idx)
                continue
            end

            boundary_spins_idx[counter] = ns
            counter += 1
        end
    end

    return boundary_spins_idx
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


function one_loop_update!(beta::Float64,x_axis,y_axis,z_axis,
                          num_accept::Int64,
                          spins::AbstractArray{HeisenbergSpin},
                          new_spins::AbstractArray{HeisenbergSpin},
                          updater::SingleSpinFlipUpdater,
                          colors::Array{Color},
                          colors_on_loop::Tuple{Color,Color},
                          first_spin_idx::Int,
                          max_length::Int,
                          work::Array{Int},
                          spin_idx_on_loop::Array{UInt},
                          verbose::Bool=false)::Tuple{Bool,Float64, Float64, Float64, Float64, Float64}
    t1 = time_ns()
    #spin_idx_on_loop = find_loop(spin_idx_on_loop, updater,colors,colors_on_loop,first_spin_idx,max_length,work,verbose)
    loop_length = find_loop(spin_idx_on_loop, updater,colors,colors_on_loop,first_spin_idx,max_length,work,verbose)
    t2 = time_ns()
    if loop_length == 0
        return false, 0.0,  t2-t1, 0.0, 0.0, 0.0
    end
    spins_on_loop = spins[spin_idx_on_loop[1:loop_length]]

    mk_new_spins_on_loop(loop_length,spins_on_loop,new_spins,colors_on_loop,x_axis,y_axis,z_axis)
    t3 = time_ns()
    dE = compute_dE_loop(updater,loop_length,spin_idx_on_loop,spins,new_spins,work,verbose)
    t4 = time_ns()
    
    #if verbose 
        #println("DEBUG: ", "loop_length", " ", "dE")
        #println("DEBUG: ",  loop_length , " ",  dE )
    #end
            
    t5 = time_ns() 
    #println("one_loop_update: find_loop spin_flip compute_dE metropolis")
    #println("one_loop_update: $(t5-t4) $(t4-t3) $(t3-t2) $(t2-t1)")
    dE = metropolis_method!(beta,dE,spins,colors,colors_on_loop,loop_length,spin_idx_on_loop,new_spins,num_accept)
    t6 = time_ns() 

    return true, dE, t2-t1, t3-t2, t4-t3, t6-t5
end
  
function multi_loop_update!(loop_updater::LoopUpdater, num_trial::Int64,num_reference::Int64,
                            updater::SingleSpinFlipUpdater,beta::Float64,
                            triangles::Array{Tuple{Int,Int,Int}},
                            max_length::Int,
                            spins::AbstractArray{HeisenbergSpin},
                            verbose::Bool=false)
    indices,x_axis,y_axis,normal_vec = estimate_loc_coord(spins,num_reference)
    if length(indices) == 0
        return 0.0, 0.0
    end

    # No copy
    colors = loop_updater.colors
    work = loop_updater.work
    spins_on_loop = loop_updater.spins_on_loop
    new_spins = loop_updater.new_spins

    timings = zeros(Float64, 3)
    dtheta  = min(50*sqrt(beta^-1),pi/3)

    t1_s = time_ns()
    mk_init_colors!(updater,spins,x_axis,y_axis,normal_vec,indices,triangles,colors,dtheta)
    t1_e = time_ns()
    timings[1] += t1_e - t1_s
    #if verbose
       #println("mk_init_color: $(1/beta) $(t1_e-t1_s)")
    #end

    dE   = 0.
    counter    = 0
    num_accept = 0 
    num_loop_found = 0

    
    for i=1:num_trial
        counter += 1
        t2_s = time_ns()
        first_spin_idx,colors_on_loop = mk_init_condition(length(spins),colors) 
        t2_e = time_ns()
        timings[2] += t2_e - t2_s
        #if verbose
           #println("mk_init_condition: $(1/beta) $(t2_e-t2_s)")
        #end
        if first_spin_idx == -1
            continue
        end

        t3_s = time_ns()
        loop_found, dE_tmp, tt1, tt2, tt3, tt4 = one_loop_update!(beta,x_axis,y_axis,normal_vec,num_accept,
               spins, new_spins, updater,colors,colors_on_loop,first_spin_idx,max_length,work,spins_on_loop,verbose) 
        t3_e = time_ns()
        #if verbose
           #println("one_loop_update $(1/beta) $(tt1) $(tt2) $(tt3) $(tt4)")
        #end
        timings[3] += t3_e - t3_s
        dE += dE_tmp
        if loop_found
            num_loop_found += 1
        end
     
        if dE_tmp != 0
            num_accept += 1
        end
        #println("execution: init_condition one_update")
        #println("execution: $(t6-t5) $(t7-t6)")
    end

    #println("multi_loop: coloring work execution")
    #println("multi_loop: $(t2-t1) $(t3-t2) $(t4-t3)")
    #if verbose
       #println("multi_loop: $(1/beta) $(timings)")
    #end
   
    return dE, num_loop_found/num_trial, num_accept/num_trial
end
