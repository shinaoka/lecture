# For unit test measurement function moved from classical_mc.jl.

include("mcmc.jl")

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


# It is important to keep computational complexity O(num_spins)
function compute_m2_af(spins::Vector{HeisenbergSpin},num_spins::Int64,
                       triangles::Array{Tuple{Int64,Int64,Int64},1})
    
    m_af = fill((0.,0.,0.),3)
    for itri in triangles
        
        for j in 1:3
            m_af[j] = m_af[j] .+ spins[itri[j]]
        end

    end
    
    m2_af = 0.
    for i in 1:3
        m2_af += sum(m_af[i].^2)
    end

    return m2_af/(3*length(triangles)^2)
end


delta(a,b) = ifelse(a==b,1,0) 

function compute_T2_op(spins::Vector{HeisenbergSpin},num_spins::Int64)

    T_op = zeros(Float64,27)
    
    for ispin in 1:num_spins
        spin = spins[ispin]
        idx  = 1
        for (a,b,c) in Iterators.product(1:3,1:3,1:3)
            
            T_op[idx] += spin[a]*spin[b]*spin[c] - (spin[a]*delta(b,c)+spin[b]*delta(c,a)+spin[c]*delta(a,b))/5
          
            idx += 1
        end

    end

    return sum(T_op.^2) / (num_spins^2)
end

