# For unit test measurement function moved from classical_mc.jl.
include("mcmc.jl")
include("loop_update.jl")

# It is important to keep computational complexity O(num_spins)
function compute_m2_af(spins::Vector{HeisenbergSpin},
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


function compute_loop_length(spins::Vector{HeisenbergSpin},
                             updater::SingleSpinFlipUpdater,
                             loop_updater::LoopUpdater,
                             max_loop_length,verbose)
     
    work = loop_updater.work
    spins_idx_on_loop = loop_updater.spins_on_loop
    new_spins_on_loop = loop_updater.new_spins

    num_spins = updater.num_spins
    max_coord_num = maximum(updater.coord_num)

    first_spin_idx = rand(1:num_spins)
    candidate_second_spin_idx = zeros(UInt,max_coord_num)
    nn_coord_num = updater.nn_coord_num[first_spin_idx]
    for ins in 1:nn_coord_num
        candidate_second_spin_idx[ins] = updater.nn_sites[ins,first_spin_idx]
    end
    second_spin_idx = rand(candidate_second_spin_idx[1:nn_coord_num])

    loop_length = 0
    loop_length,sum_boundary_spins = find_loop(spins,spins_idx_on_loop,updater,first_spin_idx,
                                                   second_spin_idx,max_loop_length,work,verbose)  

    if mod(loop_length,2) == 0
        return loop_length
    else
        return 0
    end

end

