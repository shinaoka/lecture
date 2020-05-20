# For unit test measurement function moved from classical_mc.jl.

include("mcmc.jl")

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

    return 6*m2_af/(num_spins^2)
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

