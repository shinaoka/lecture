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

function estimate_axes(spin::HeisenbergSpin,vec::Array{Float64,1})
    #=
    estimate x and y axis in a plane determined by the vec as (normal) vectors.
    =#
    
    x_axis  = collect(spin) - dot(collect(spin),vec)*vec
    x_axis /= norm(x_axis)
    y_axis  = cross(vec,x_axis)

    return x_axis,y_axis
end

@enum Color red blue green black

function mk_reference(spins::AbstractArray{HeisenbergSpin},num_reference::Int64)
    #= 
    Make a reference system for estimate local spin coordination.
    =#

    rand_idx  = [rand(1:length(spins)) for i=1:num_reference] 
    reference = copy(spins[rand_idx])

    return reference,rand_idx
end

function paint_black!(colors::Array{Color},indices)
    
    colors[indices] = [black for i=1:length(indices)]
end

   
function paint_rbg_differently!(spins::AbstractArray{HeisenbergSpin},x_axis,y_axis,normal_vec,colors::Array{Color})
    #=
    rbg = red blue green. paint spins rbg based on its direction in local spin coordination.
    =#    
    
     
    num_spins = length(spins)

    for i=1:num_spins
        
        if colors[i] !== black
                
            spin = collect(spins[i])
            proj_unit_vec  = spin - (dot(spin,normal_vec))*normal_vec
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

    end
 
end    

function find_breaking_triangle!(updater::SingleSpinFlipUpdater,colors::Array{Color})
    #= 
    breaking triangle has more than two site painted the same color.
    =#
    
    #First,find nearest neighbor pair which have same color.
    same_color_pairs = []
    
    for idx=1:length(colors)
        for ins=1:updater.coord_num[idx]
            ns = updater.connection[ins,idx][1]
            isnn = updater.connection[ins,idx][5] == 1
            if isnn && colors[idx]==colors[ins] 
                push!(same_color_pairs,(idx,ns))  
            end
        end
    end
    
    #Second,find nearest neighbor site from each same color pairs. 
    breaking_triangle = []
    for idx=1:length(same_color_pairs)
        i,j = same_color_pairs[idx]

        nnset_i = []
        for ins=1:updater.coord_num[i]
            ns = updater.connection[ins,i][1]
            isnn = updater.connection[ins,i][5] == 1
            if isnn 
                push!(nnset_i,ns)  
            end
        end

        nnset_j = []
        for ins=1:updater.coord_num[j]
            ns = updater.connection[ins,j][1]
            isnn = updater.connection[ins,j][5] == 1
            if isnn 
                push!(nnset_j,ns)  
            end
        end
       
        for temp in nnset_i
            if in(temp,nnset_j)
                push!(breaking_triangle,(i,j,temp)) 
            end
        end

    end
    
    #Finally,assign black to sites on breaking triangles.
    for i=1:length(breaking_triangle)
        for j in collect(breaking_triangle[i])
            colors[j] = black
        end
    end
 
end

function parallel_flip(spin::HeisenbergSpin,color::Color,x_axis,y_axis,z_axis)

    # Color typed variable red=0 blue=1 green=2 black=3.
    theta = Int(color) * 2pi/3
    mirror_vec = cos(theta) * x_axis + sin(theta) * y_axis
 
    spin = collect(spin)
    spin_on_plane = spin - (dot(spin,z_axis)) * z_axis
    
    # parallel flip
    temp = (dot(spin_on_plane,mirror_vec)) * mirror_vec 
    flipped_spin = 2 * temp - spin_on_plane

    return Tuple(flipped_spin + (dot(spin,z_axis)) * z_axis)
end

function mk_new_spins_on_loop(spins_on_loop,colors_on_loop::Tuple{Color,Color},x_axis,y_axis,z_axis)
    
    # find the color there are not on a loop
    color = red
    for ic in [red,blue,green]
        if !in(ic,collect(colors_on_loop))
            color = ic
        end
    end

    loop_length       = length(spins_on_loop)
    new_spins_on_loop = fill((0.,0.,0.),loop_length)

    for i=1:loop_length
        new_spins_on_loop[i] = parallel_flip(spins_on_loop[i],color,x_axis,y_axis,z_axis) 
    end
    return new_spins_on_loop
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

function update_colors!(colors::Array{Color},colors_on_loop::Tuple{Color,Color},spin_idx_on_loop::Array{UInt})

    for idx=spin_idx_on_loop
        colors[idx] == colors_on_loop[1] ? colors[idx] = colors_on_loop[2] : colors[idx] = colors_on_loop[1]
    end
  
end

function metropolis_method!(beta::Float64,dE::Float64,
                            spins::AbstractArray{HeisenbergSpin},
                            colors::Array{Color},
                            colors_on_loop::Tuple{Color,Color},
                            spin_idx_on_loop::Array{UInt},
                            new_spins_on_loop::Array{HeisenbergSpin})

    loop_length = length(spins_on_loop)

    if rand() < exp(-beta*dE)
        
        spins[spin_idx_on_loop] = new_spins_on_loop 
        update_colors!(colors,colors_on_loop,spin_idx_on_loop)
        return dE
    
    else
        return 0
    end
    
end

#################################################################

function estimate_loc_coord(spins,num_reference)
    
    reference,indices = mk_reference(spins,num_reference)
    normal_vec = estimate_plane(reference)
    x_axis,y_axis = estimate_axes(reference[1],normal_vec)
  
    return indices,x_axis,y_axis,normal_vec
end

function mk_init_colors(updater::SingleSpinFlipUpdater,spins::AbstractArray{HeisenbergSpin},x_axis,y_axis,z_axis,indices)

   colors = [red for i=1:length(spins)]
   paint_black!(colors,indices)

   paint_rbg_differently!(spins,x_axis,y_axis,z_axis,colors) 
   find_breaking_triangle!(updater,colors) 
   
   return colors
end

function one_loop_update!(beta::Float64,x_axis,y_axis,z_axis,
                         spins::AbstractArray{HeisenbergSpin},
                         updater::SingleSpinFlipUpdater,
                         colors::Array{Color},
                         colors_on_loop::Tuple{Color,Color},
                         first_spin_idx::Int,
                         max_length::Int,
                         work::Array{Int},
                         check::Bool=false)

    spin_idx_on_loop = find_loop(updater,colors,colors_on_loop,first_spin_idx,max_length,work,check)
    
    spins_on_loop = copy(spins[spin_idx_on_loop])

    new_spins_on_loop = mk_new_spins_on_loop(spins_on_loop,colors_on_loop,x_axis,y_axis,z_axis)
    dE = compute_dE_loop(updater,spin_idx_on_loop,spins,new_spins_on_loop,work,check)
    
    return  metropolis_method!(beta,spins,colors,colors_on_loop,spin_idx_on_loop,new_spins_on_loop)

end


