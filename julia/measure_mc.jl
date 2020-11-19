# For unit test measurement function moved from classical_mc.jl.
using FFTW
include("mcmc.jl")
include("loop_update.jl")

#For more efficient memory access, convert Vector to two-dimensional Array.
function convert_spins_to_array(spins::Vector{HeisenbergSpin})
    N = length(spins)
    spins_array = Array{Float64,2}(undef, 3, N)
    convert_spins_to_array!(spins, spins_array)
    spins_array
end

function convert_spins_to_array!(spins::Vector{HeisenbergSpin},spins_array::AbstractMatrix{Float64})
    N = length(spins)
    for i in 1:N, j=1:3
        spins_array[j,i] = spins[i][j]
    end
end

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


function compute_mq(q,kagome,spins,triangles)
    @assert length(keys(kagome)) == length(spins)

    sq = fill((0.0+0im,0.0+0im,0.0+0im),3)
    for i in triangles
        for j in 1:3
            idx = i[j]
            R = kagome[idx]
            sq[j] = sq[j] .+ spins[idx].*exp(2π*(q⋅R)*im)
        end
    end

    ss = 0.

    for i in 1:3
        ss += sq[i] ⋅ sq[i]
    end

    return real(ss)/(3*length(triangles)^2)
end

function compute_m_120degrees(spins)
    cos_sum = sum((cos(3*atan(x[2], x[1])) for x in spins))
    sin_sum = sum((sin(3*atan(x[2], x[1])) for x in spins))
    (cos_sum^2 + sin_sum^2)/length(spins)^2
end


function compute_sisj(num_triangles, spins, triangles)
    num_spins = length(spins)
    sisj = Array{Float64,3}(undef, num_spins, 3, num_triangles)
    for it in 1:num_triangles
        for isublatt in 1:3
            for ispin in 1:num_spins
               sisj[ispin,isublatt,it] = dot(spins[ispin], spins[triangles[it][isublatt]])
            end
        end
    end
    sisj
end

# vector spin chirality defined as sum of outer product of spins in unit triangular with counterclockwise rotation.
#=
function compute_vector_chirality(spins::Vector{HeisenbergSpin},
                                  triangles::Array{Tuple{Int64,Int64,Int64},1})

    vc = 0.0
    num_spins = length(spins)
    num_triangles = length(triangles)
    @assert num_spins == 3num_triangles
    spins = collect.(spins)
    for i in triangles
        for j in 1:3
            s1 = spins[i[j]]
            s2 = spins[i[ifelse(j==3,1,j+1)]]
            vc += s1[1]s2[2] - s1[2]s2[1]
            #vc += cross(s[j],s[ifelse(j==3,1,j+1)])[3]
        end
    end

    return vc / num_spins
end
=#

mycross(s1, s2) = s1[1]s2[2] - s1[2]s2[1]

function compute_vector_chirality(spins::AbstractArray{Float64,2},
                                  triangles::Vector{Tuple{Int64,Int64,Int64}})

    vc = 0.0
    num_spins = size(spins)[2]
    num_triangles = length(triangles)
    @assert num_spins == 3num_triangles
    for i in triangles
        for j in 1:3
            s1 = view(spins, :, i[j])
            s2 = view(spins, :, i[ifelse(j==3,1,j+1)])
            vc += mycross(s1, s2)
        end
    end
    return vc / num_spins
end

function compute_all_vector_chiralities(spins::AbstractArray{Float64,2},
                                  triangles::Vector{Tuple{Int64,Int64,Int64}})
    num_spins = size(spins)[2]
    num_triangles = length(triangles)
    vc = zeros(Float64, num_triangles)
    @assert num_spins == 3num_triangles
    for it in eachindex(triangles)
        i = triangles[it]
        for j in 1:3
            s1 = view(spins, :, i[j])
            s2 = view(spins, :, i[ifelse(j==3,1,j+1)])
            vc[it] += mycross(s1, s2)
        end
    end
    vc
end

function compute_vector_chiralities(spins::AbstractArray{Float64,2}, utriangles, dtriangles)
    @assert length(utriangles) == length(dtriangles)
    num_spins = size(spins)[2]
    uc_all = compute_all_vector_chiralities(spins, utriangles)
    dc_all = compute_all_vector_chiralities(spins, dtriangles)
    uc, dc = sum(uc_all)/num_spins, sum(dc_all)/num_spins
    (uc+dc)^2/3, (uc-dc)^2/3, uc_all[1] * [uc_all; dc_all]
end

function compute_ferro_vector_chirality(spins::AbstractArray{Float64,2}, utriangles, dtriangles)
    @assert length(utriangles) == length(dtriangles)
    fvc = compute_vector_chirality(spins,utriangles) + compute_vector_chirality(spins,dtriangles)
    fvc ^= 2
    return fvc / 3 # vector spin chirality of q=0 state is 3,larger than that of all other states.
end

function compute_af_vector_chirality(spins::AbstractArray{Float64,2}, utriangles, dtriangles)
    @assert length(utriangles) == length(dtriangles)
    fvc = compute_vector_chirality(spins,utriangles) - compute_vector_chirality(spins,dtriangles)
    fvc ^= 2
    return fvc / 3 # vector spin chirality of √3×√3 state is 3,larger than that of all other states.
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

    num_spins = updater.num_spins
    max_coord_num = maximum(updater.coord_num)

    first_spin_idx = rand(1:num_spins)
    candidate_second_spin_idx = zeros(UInt,max_coord_num)
    nn_coord_num = updater.nn_coord_num[first_spin_idx]
    for ins in 1:nn_coord_num
        candidate_second_spin_idx[ins] = updater.nn_sites[ins,first_spin_idx]
    end
    second_spin_idx = rand(candidate_second_spin_idx[1:nn_coord_num])

    loop_length,sum_boundary_spins = find_loop(spins,spins_idx_on_loop,updater,first_spin_idx,
                                                   second_spin_idx,max_loop_length,work,verbose)  

    if mod(loop_length,2) == 0
        return loop_length
    else
        return 0
    end

end


function compute_1dcorr(spins::Vector{IsingSpin})
    
    num_spins = length(spins)
    spins_q = fft(spins)
    corr_q  = spins_q .* conj(spins_q)
    
    return real(ifft(corr_q)) / num_spins
end


function compute_corr(spins::Vector{HeisenbergSpin})


end


function compute_corr_length(corr,rdata,param_init)
  
    model(t,p) = p[1]*exp.(-t/p[2])
    fit = curve_fit(model,rdata,corr,param_init)
    
    return fit.param[2]
end

