include("mcmc.jl")

using LinearAlgebra
using StaticArrays
using Random

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

    if norm(normal_vec) > 1e-10
        normal_vec /= norm(normal_vec)
        return normal_vec
    else
        return [0.0, 0.0, 0.0]
    end
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

function estimate_axes(spins::Array{HeisenbergSpin}, vec::Array{Float64})
    #=
    estimate x and y axis in a plane determined by the vec as (normal) vectors.
    =#

    function project(spin)
       tmp = collect(spin) - dot(collect(spin),vec)*vec
       return tmp/norm(tmp)
    end
    
    x_axis0  = project(spins[1])
    x_axis = sum([project(spins[i]) for i in eachindex(spins) if dot(project(spins[i]), x_axis0) > 1/2])
    x_axis /= norm(x_axis)

    y_axis  = cross(vec,x_axis)

    return x_axis,y_axis
end

@enum Color red blue green black

function mk_reference(spins::AbstractArray{HeisenbergSpin},num_reference::Int64)
    #= 
    Make a reference system for estimate local spin coordination.
    =#

    #rand_idx  = [rand(1:length(spins)) for i=1:num_reference] 
    rand_idx  = rand(1:length(spins), num_reference)
    reference = spins[rand_idx]

    return reference,rand_idx
end

function paint_black!(colors::Array{Color},indices)
    
    for i in indices
        colors[i] = black
    end

end

   
function paint_rbg_differently!(spins::AbstractArray{HeisenbergSpin},x_axis,y_axis,normal_vec,colors::Array{Color})
    #=
    rbg = red blue green. paint spins rbg based on its direction in local spin coordination.
    =#    
     
    # Allocate static arrays as work space (much faster than using dynamic arrays)
    xvec = SVector(x_axis[1], x_axis[2], x_axis[3])
    yvec = SVector(y_axis[1], y_axis[2], y_axis[3])
    zvec = SVector(normal_vec[1], normal_vec[2], normal_vec[3])
    work = MVector(0.0, 0.0, 0.0)
    
    for i in eachindex(spins)
        if colors[i] == black
            continue
        end
        for c = 1:3
            work[c] = spins[i][c]
        end
        work[:]  = normalize(work - (dot(work,zvec))*zvec)
        if dot(work, xvec) >= 1/2
            colors[i] = red
        else
            colors[i] = dot(work, yvec) >= 0 ? blue : green
        end
    end
end    

function paint_rbg_differently_v2!(spins::AbstractArray{HeisenbergSpin},x_axis,y_axis,normal_vec,colors::Array{Color},dtheta::Float64)
    #=
    rbg = red blue green. paint spins rbg based on its direction in local spin coordination.
    =#    
     
    # Allocate static arrays as work space (much faster than using dynamic arrays)
    xvec = SVector(x_axis[1], x_axis[2], x_axis[3])
    yvec = SVector(y_axis[1], y_axis[2], y_axis[3])
    zvec = SVector(normal_vec[1], normal_vec[2], normal_vec[3])
    work = MVector(0.0, 0.0, 0.0)

    for i in eachindex(spins)
        if colors[i] == black
            continue
        end
        for c = 1:3
            work[c] = spins[i][c]
        end
        work[:]  = normalize(work - (dot(work,zvec))*zvec)
        if dot(work, xvec) >= cos(dtheta)
            colors[i] = red
        elseif cos(2pi/3-dtheta)<= dot(work,xvec) <= cos(2pi/3+dtheta) 
            if dot(work,yvec) > 0
                colors[i] = blue
            else 
                colors[i] = green
            end
        else
            colors[i] = black
        end
    end
end    

function paint_rbg_differently_v3!(spins::AbstractArray{HeisenbergSpin},x_axis,y_axis,normal_vec,colors::Array{Color},dtheta::Float64)
    #=
    rbg = red blue green. paint spins rbg based on its direction in local spin coordination.
    =#    
     
    # Allocate static arrays as work space (much faster than using dynamic arrays)
    xvec = SVector(x_axis[1], x_axis[2], x_axis[3])
    yvec = SVector(y_axis[1], y_axis[2], y_axis[3])
    zvec = SVector(normal_vec[1], normal_vec[2], normal_vec[3])
    work = MVector(0.0, 0.0, 0.0)
    
     
    for i in eachindex(spins)
        if colors[i] == black
            continue
        end
        for c = 1:3
            work[c] = spins[i][c]
        end
        work[:]  = normalize(work - (dot(work,zvec))*zvec)
        x = dot(work,xvec)
        y = dot(work,yvec)
        temp_angle = atan(y,x)

        if - dtheta <= temp_angle <= dtheta
            colors[i] = red
        elseif  2pi/3 - dtheta <= temp_angle <=  2pi/3 + dtheta
            colors[i] = blue
        elseif -2pi/3 - dtheta <= temp_angle <= -2pi/3 + dtheta
            colors[i] = green 
        else
            colors[i] = black
        end
    end
end    

function find_triangles(model::JModel, updater::SingleSpinFlipUpdater)
    #= 
    Find all triangles
    =#
    num_spins = updater.num_spins

    # Creat a set of nn sites for each site
    nn_sites = [Set{Int}() for i in 1:num_spins]
    for isite=1:num_spins
        for ins=1:updater.coord_num[isite]
            if updater.connection[ins,isite][5] == 1
                push!(nn_sites[isite], updater.connection[ins,isite][1])
            end
        end
    end

    triangles = Set{Tuple{Int,Int,Int}}()
    for (i, j, Jx, Jy, Jz, is_nn) in model.Jij
        if is_nn != 1
            continue
        end
        common_nn = intersect(nn_sites[i], nn_sites[j])
        @assert length(common_nn) <= 1
        if length(common_nn) != 0
            t_sites = sort([i, j, pop!(common_nn)])
            push!(triangles, Tuple(t_sites))
        end
    end
    return collect(triangles)
end

function find_breaking_triangle!(updater::SingleSpinFlipUpdater, triangles::Array{Tuple{Int,Int,Int}}, colors::Array{Color})
    #= 
    breaking triangle has more than two site painted the same color.
    =#
    is_black = zeros(Bool, length(colors))
    rgb::UInt8 = 0
    rgb |= (1 << UInt8(red))
    rgb |= (1 << UInt8(green))
    rgb |= (1 << UInt8(blue))
    for three_sites in triangles
        x::UInt8 = 0
        for isite in three_sites
            x |= (1 << UInt8(colors[isite]))
        end
        if x != rgb
            for isite in three_sites
                is_black[isite] = true
            end
        end
    end
    for isite in eachindex(colors)
        if is_black[isite]
            colors[isite] = black
        end
    end
end    

function parallel_flip(loop_length::Int,
    spins::Array{HeisenbergSpin},
    new_spins::Array{HeisenbergSpin},
    color::Color,x_axis,y_axis,z_axis)

    xvec = SVector(x_axis[1], x_axis[2], x_axis[3])
    yvec = SVector(y_axis[1], y_axis[2], y_axis[3])
    zvec = SVector(z_axis[1], z_axis[2], z_axis[3])
    spin = MVector(0.0, 0.0, 0.0)
    spin_on_plane = MVector(0.0, 0.0, 0.0)

    # Color typed variable red=0 blue=1 green=2 black=3.
    theta = Int(color) * 2pi/3
    mirror_vec = cos(theta) * xvec + sin(theta) * yvec

    for i in 1:loop_length
        for ix in 1:3
           spin[ix] = spins[i][ix]
        end
        spin_on_plane[:] = spin - (dot(spin,zvec)) * zvec

        # parallel flip
        temp = (dot(spin_on_plane, mirror_vec)) * mirror_vec 
        flipped_spin = 2 * temp - spin_on_plane
        new_spin = flipped_spin + (dot(spin, zvec)) * zvec
        new_spins[i] = (new_spin[1], new_spin[2], new_spin[3])
    end
end

function mk_new_spins_on_loop(loop_length, spins_on_loop,
    new_spins_on_loop,
    colors_on_loop::Tuple{Color,Color},x_axis,y_axis,z_axis)
    
    # find the color there are not on a loop
    color = red
    for ic in [red,blue,green]
        if ic != colors_on_loop[1] && ic != colors_on_loop[2]
            color = ic
        end
    end

    #new_spins_on_loop = fill((0.,0.,0.),loop_length)
    parallel_flip(loop_length, spins_on_loop, new_spins_on_loop, color,x_axis,y_axis,z_axis) 
    #return new_spins_on_loop
end

function find_loop(spins_on_loop, updater::SingleSpinFlipUpdater, colors::Array{Color},
    colors_on_loop::Tuple{Color,Color}, first_spin_idx::Int, max_length::Int, work::Array{Int}, check::Bool=false)::Int64
    #=
    All elements of work must be initialized to zero.
    =#
    num_spins = updater.num_spins

    # Optional expensive check
    @assert length(work) >= num_spins
    @assert colors[first_spin_idx] == colors_on_loop[1] || colors[first_spin_idx] == colors_on_loop[2]

    work[first_spin_idx] = 1
    spins_on_loop[1] = first_spin_idx
    loop_length::Int64 = 1
    spin_before::Int64 = -1
    current_spin_idx::Int64 = first_spin_idx

    max_coord_num = maximum(updater.coord_num)
    candidate_spins = zeros(UInt, max_coord_num)

    #if check
       #println("colors_on_loop $(colors_on_loop)")
    #end
    success = false
    status = -1
    while loop_length < max_length
        #if check
           #println("current_spin_idx $(current_spin_idx)")
        #end
        # Search connected spins
        n_candidate = 0
        next_color = colors[current_spin_idx]==colors_on_loop[1] ? colors_on_loop[2] : colors_on_loop[1]
        for ins in 1:updater.coord_num[current_spin_idx]
            # candidate must be either the first spin on the loop (work[ns]==1) or unvisited (work[ns]==0)
            ns::SpinIndex = updater.connected_sites[ins,current_spin_idx]
            # Is ns a nearest neighborb site from current site?
            if updater.is_NN[ins,current_spin_idx] && colors[ns]==next_color && work[ns] <= 1 && ns != spin_before
                n_candidate += 1
                candidate_spins[n_candidate] = ns
            end
        end
        if n_candidate == 0
            status = 1
            break
        end

        next_spin_idx = candidate_spins[rand(1:n_candidate)]
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
    #if check
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

   

function compute_dE_loop(updater::SingleSpinFlipUpdater,
                          loop_length::Int,
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
    #if check
        #@assert all(work .== 0)
    #end

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

    # Optional expensive check
    #if check
        #@assert all(work .== 0)
    #end

    return dE
end

function update_colors!(colors::Array{Color},colors_on_loop::Tuple{Color,Color},loop_length::Int,spin_idx_on_loop::Array{UInt})

    for idx = spin_idx_on_loop[1:loop_length]
        colors[idx] == colors_on_loop[1] ? colors[idx] = colors_on_loop[2] : colors[idx] = colors_on_loop[1]
    end
  
end

function metropolis_method!(beta::Float64,dE::Float64,
                            spins::AbstractArray{HeisenbergSpin},
                            colors::Array{Color},
                            colors_on_loop::Tuple{Color,Color},
                            loop_length::Int,
                            spin_idx_on_loop::Array{UInt},
                            new_spins_on_loop::Array{HeisenbergSpin},
                            num_accept::Int64)::Float64
    
    if rand(Random.GLOBAL_RNG) < exp(-beta*dE)
        spins[spin_idx_on_loop[1:loop_length]] = new_spins_on_loop[1:loop_length]
        update_colors!(colors,colors_on_loop,loop_length,spin_idx_on_loop)
        num_accept += 1
        
        return dE
    else
        #p = exp(-beta*dE)
        #println("update failed: $(dE) $(p) $(dE*beta)")
        return 0.0
    end
end


function estimate_loc_coord(spins,num_reference)
    reference,indices = mk_reference(spins,num_reference)
    normal_vec = estimate_plane(reference)
    if norm(normal_vec) > 1e-10
        x_axis,y_axis = estimate_axes(reference, normal_vec)
        return indices, x_axis, y_axis, normal_vec
    else
        return Array{Int}[], Array{Float64}[], Array{Float64}[], Array{Float64}[]
    end
end

function mk_init_colors!(updater::SingleSpinFlipUpdater,spins::AbstractArray{HeisenbergSpin},x_axis,y_axis,z_axis,indices,triangles::Array{Tuple{Int,Int,Int}},colors::Array{Color},dtheta::Float64)
   for i in eachindex(colors)
       colors[i] = red
   end
   t1 = time_ns()
   paint_black!(colors,indices)
   t2 = time_ns()
   paint_rbg_differently_v3!(spins,x_axis,y_axis,z_axis,colors,dtheta) 
   t3 = time_ns()
   find_breaking_triangle!(updater,triangles,colors) 
   t4 = time_ns()
   #println("mk_init_c ", "paint black", " ", "paint rbg", " ", "find b-triangles")
   #println("mk_init_c ", t2-t1, " ", t3-t2, " ", t4-t3)
   
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
                          check::Bool=false)::Tuple{Bool,Float64, Float64, Float64, Float64, Float64}
    t1 = time_ns()
    #spin_idx_on_loop = find_loop(spin_idx_on_loop, updater,colors,colors_on_loop,first_spin_idx,max_length,work,check)
    loop_length = find_loop(spin_idx_on_loop, updater,colors,colors_on_loop,first_spin_idx,max_length,work,check)
    t2 = time_ns()
    if loop_length == 0
        return false, 0.0,  t2-t1, 0.0, 0.0, 0.0
    end
    spins_on_loop = spins[spin_idx_on_loop[1:loop_length]]

    mk_new_spins_on_loop(loop_length,spins_on_loop,new_spins,colors_on_loop,x_axis,y_axis,z_axis)
    t3 = time_ns()
    dE = compute_dE_loop(updater,loop_length,spin_idx_on_loop,spins,new_spins,work,check)
    t4 = time_ns()
    
    #if check 
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

function mk_init_condition(num_spins::Int64,colors::Array{Color})
    
    first_spin_idx = -1
    color1 = black

    first_spin_idx = rand(1:num_spins)
    color1 = colors[first_spin_idx]

    if color1 != black
        temp_colors = Set([red, blue, green])
        delete!(temp_colors, color1)
        return first_spin_idx, (color1, rand(temp_colors))
    end

    if color1 == black
        return -1, (black, black)
    end
end

struct LoopUpdater{T}
    num_spins::Int
    colors::Array{Color}
    work::Array{Int}
    spins_on_loop::Array{UInt}
    new_spins::Array{T}
end

function LoopUpdater{T}(num_spins::Int64, max_loop_length::Int64) where T
    colors = fill(red,num_spins)
    work   = zeros(Int,num_spins)
    spins_on_loop  = zeros(Int, max_loop_length)
    new_spins  = Array{T}(undef, max_loop_length)
    return LoopUpdater{T}(num_spins,colors,work, spins_on_loop, new_spins)
end

  
function multi_loop_update!(loop_updater::LoopUpdater, num_trial::Int64,num_reference::Int64,
                            updater::SingleSpinFlipUpdater,beta::Float64,
                            triangles::Array{Tuple{Int,Int,Int}},
                            max_length::Int,
                            spins::AbstractArray{HeisenbergSpin},
                            check::Bool=false)
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
    #if check
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
        #if check
           #println("mk_init_condition: $(1/beta) $(t2_e-t2_s)")
        #end
        if first_spin_idx == -1
            continue
        end

        t3_s = time_ns()
        loop_found, dE_tmp, tt1, tt2, tt3, tt4 = one_loop_update!(beta,x_axis,y_axis,normal_vec,num_accept,
               spins, new_spins, updater,colors,colors_on_loop,first_spin_idx,max_length,work,spins_on_loop,check) 
        t3_e = time_ns()
        #if check
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
    #if check
       #println("multi_loop: $(1/beta) $(timings)")
    #end
   
    return dE, num_loop_found/num_trial, num_accept/num_trial
end
