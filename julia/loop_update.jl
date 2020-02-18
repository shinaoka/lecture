include("mcmc.jl")

using LinearAlgebra

function estimate_plane(spins::AbstractArray{HeisenbergSpin})
    #=
    Estimate a coplanar plane and return its (normalized) normal vector
    =#

    num_spins = length(spins)

    # Normal vector perpendicular to the plane
    normal_vec = zeros(Float64, 3)

    for i=1:num_spins
        for j=i+1:num_spins
           cv = cross(collect(spins[i]), collect(spins[j]))
           if dot(cv, normal_vec) < 0
               cv *=  -1
           end
           normal_vec += cv
        end
    end
    normal_vec /= norm(normal_vec)

    return normal_vec
end

function estimate_axes(spins::AbstractArray{HeisenbergSpin})
    #=
    argument spins is random sampled array from original spins
    =#
    
    norm_vec   = estimate_plane(spins)
    first_spin = collect(spins[1])
    
    x_axis  = first_spin - dot(first_spin,norm_vec)*norm_vec
    x_axis /= norm(x_axis)
    y_axis  = cross(norm_vec,x_axis)

    return x_axis,y_axis,norm_vec
end

function mk_reference_system(spins::AbstractArray{HeisenbergSpin},num_rand_spins)
 
    num_spins  = length(spins)
    
    # random spin sampling and paint them black
    rand_spins = fill((0.,0.,0.),num_rand_spins) 
    
    for i=1:num_rand_spins
        rand_idx = rand(1:num_spins)
        rand_spins[i] = spins[rand_idx]
    end

    return rand_spins
end
 
@enum Color red blue green black

# paint red blue green differently on kagome site.
# For simple test,pass this function reference and color array as a argument.
function paint_rbg_differently(spins::AbstractArray{HeisenbergSpin},reference_system)
    
    x_axis,y_axis,norm_vec = estimate_axes(reference_system)
    
    num_spins = length(spins)
    colors = Color.([0 for i=1:num_spins])

    # assign three colors to the site 
    for i=1:num_spins
        
        spin = collect(spins[i])
        proj_unit_vec  = spin - (dot(spin,norm_vec))*norm_vec
        proj_unit_vec /= norm(proj_unit_vec)  
        
        if dot(proj_unit_vec,x_axis) >= 1/2
                 colors[i] = red
        else
            if dot(proj_unit_vec,y_axis) >= 0
                 colors[i] = blue
            else
                 colors[i] = green
            end
        end

    end
 
    return colors 
end    
    

function find_loop(updater::SingleSpinFlipUpdater, colors::Array{Color}, colors_on_loop::Tuple{Color,Color}, first_spin_idx::Int, max_length::Int, work::Array{Int}, check::Bool=false)
    #=
    All elements of work must be initialized to zero.
    =#
    num_spins = updater.num_spins

    # Optional expensive check
    if check
        @assert all(work .== 0)
    end

    @assert length(work) >= num_spins
    @assert colors[first_spin_idx] == colors_on_loop[1] || colors[first_spin_idx] == colors_on_loop[2]

    work[first_spin_idx] = 1
    spins_on_loop = zeros(UInt, max_length)
    spins_on_loop[1] = first_spin_idx
    loop_length = 1
    spin_before = -1
    current_spin_idx = first_spin_idx

    max_coord_num = maximum(updater.coord_num)
    candidate_spins = zeros(UInt, max_coord_num)

    success = false
    while loop_length <= max_length
        # Search connected spins
        n_candidate = 0
        next_color = colors[current_spin_idx]==colors_on_loop[1] ? colors_on_loop[2] : colors_on_loop[1]
        for ins=1:updater.coord_num[current_spin_idx]
            # candidate must be either the first spin on the loop (work[ns]==1) or unvisited (work[ns]==0)
            ns = updater.connection[ins,current_spin_idx][1]
            # Is ns a nearest neighborb site from current site?
            isnn = updater.connection[ins,current_spin_idx][5] == 1
            if isnn && colors[ns]==next_color && work[ns] <= 1 && ns != spin_before
                n_candidate += 1
                candidate_spins[n_candidate] = ns
            end
        end
        if n_candidate == 0
            break
        end

        next_spin_idx = candidate_spins[rand(1:n_candidate)]
        if work[next_spin_idx] == 1
            # OK, we've returned to the starting point.
            success = true
            break
        end

        spin_before = current_spin_idx
        current_spin_idx = next_spin_idx
        loop_length += 1
        work[current_spin_idx] = loop_length
        spins_on_loop[loop_length] = current_spin_idx
    end

    # Reset all elements of work to 0
    for l=1:loop_length
        work[spins_on_loop[l]] = 0
    end
    # Optional expensive check
    if check
        @assert all(work .== 0)
    end

    if success
        return spins_on_loop[1:loop_length]
    else
        # Return the 0-length array of the same for type-stability
        return zeros(UInt,0)
    end
end


function compute_dE_loop(updater::SingleSpinFlipUpdater,
                          spin_idx_on_loop::Array{UInt},
                          spins::Array{HeisenbergSpin},
                          new_spins_on_loop::Array{HeisenbergSpin},
                          work::Array{Int},
                          check::Bool=false)
    #=
    Compute change in energy
    All elements of work must be initialized to zero.
    =#

    # Optional expensive check
    if check
        @assert all(work .== 0)
    end

    num_spins = updater.num_spins
    num_spins_loop = length(spin_idx_on_loop)

    for isp_loop in 1:num_spins_loop
        work[spin_idx_on_loop[isp_loop]] = isp_loop
    end

    dE = 0.0
    for isp_loop in 1:num_spins_loop
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

    for isp_loop in 1:num_spins_loop
        work[spin_idx_on_loop[isp_loop]] = 0
    end

    # Optional expensive check
    if check
        @assert all(work .== 0)
    end

    return dE
end
