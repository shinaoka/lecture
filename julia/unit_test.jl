using Test

include("measure_mc.jl")
include("loop_update_tests.jl")

println("unit test results")

function x_model(num_spins)

    return [fill((1.,0.,0.),num_spins)]
end

function random_model(num_spins)    
    spins = fill((0.,0.,0.),num_spins)
    for i in 1:num_spins
        spins[i] = Tuple(normalize(rand(3)))
    end

    return [spins]
end

# octopolar_v2 and order_parameter are former implementations measuring order parameters.
function octopolar_v2(spins,num_spins::Int64,num_temps::Int64)

    T = zeros(num_temps)

    for i in 1:num_temps
        temp = zeros(3^3)

        for j in 1:num_spins

            spin = spins[i][j]
            index = 1

            for (a,b,c) in Iterators.product(1:3,1:3,1:3)
                temp[index] +=  spin[a]*spin[b]*spin[c] - (spin[a]*delta(b,c) + spin[b]*delta(c,a) + spin[c]*delta(a,b))/5
                index += 1
            end
        end
        T[i] = sum(temp.^2)

    end
    return T/(num_spins^2)
end


# It looks like wrong,need to be fixed.
function order_parameter(spins::Vector{Vector{Tuple{Float64,Float64,Float64}}},
                         num_spins::Int64, num_temps::Int64, q::Tuple{Float64,Float64}, unit_cell::Vector{Tuple{Float64,Float64}})
    
    #unit_cell = make_kagome(num_spins) 

    M2_AF = zeros(num_temps) 

    for temp in 1:num_temps
        gp1 = spins[temp][1:3:num_spins] 
        gp2 = spins[temp][2:3:num_spins] 
        gp3 = spins[temp][3:3:num_spins] 
        
        temp_vec1 = (0.,0.,0.)
        temp_vec2 = (0.,0.,0.)
        temp_vec3 = (0.,0.,0.)

        for i in 1:Int(num_spins/3)

            phase_iqr = -im*dot(q,unit_cell[i])
   
            temp_vec1 = temp_vec1 .+ gp1[i] .* exp(phase_iqr) 
            temp_vec2 = temp_vec2 .+ gp2[i] .* exp(phase_iqr) 
            temp_vec3 = temp_vec3 .+ gp2[i] .* exp(phase_iqr) 
        end

        M2_AF[temp] = norm(temp_vec1)^2 + norm(temp_vec2)^2 + norm(temp_vec3)^2
    end
    
    return M2_AF
end


function test_compute_m2_af()
    
    num_spins = 6
    spins = x_model(num_spins) 
   
    triangles = [(1,2,3),(4,5,6)]
    m2_af = compute_m2_af(spins[1],triangles)
    @test isapprox(m2_af,1.0)

end


function test_compute_T2_op(spins::Vector{Vector{HeisenbergSpin}})
    
    num_spins = length(spins[1])
    T2_op = compute_T2_op(spins[1],num_spins)
    @test isapprox(T2_op,octopolar_v2(spins,num_spins,1)[1])

end


function test_compute_1dcorr(spins::Vector{IsingSpin})
   
    num_spins = length(spins)

    corr_rough = zeros(Float64,num_spins) 
    for i in 1:num_spins
        for j in 1:num_spins
            if num_spins < i+j < 2*num_spins
                corr_rough[j] += spins[i]*spins[mod(i+j,num_spins)]
                continue
            elseif i+j == 2*num_spins
                corr_rough[num_spins] += spins[num_spins]^2
                continue
            end
        corr_rough[j] += spins[i]*spins[i+j]
        end
    end
    prepend!(corr_rough,corr_rough[end])
    pop!(corr_rough)
    corr_rough /= num_spins

    corr_wfft  = compute_1dcorr(spins) 
  
    #println("rough: ",corr_rough)
    #println("wfft: ",corr_wfft)
  
    @test all(isapprox.(corr_rough,corr_wfft))
 
    rdata = [(i-1) for i in 1:num_spins]
    param_init = [0.5,0.5]
    corr_length_wfft  = compute_corr_length(corr_wfft ,rdata,param_init)
    corr_length_rough = compute_corr_length(corr_rough,rdata,param_init)

    println(corr_length_wfft)
    println(corr_length_rough)
  
    @test isapprox(corr_length_wfft,corr_length_rough)

end


function test_compute_corr(spins::Vector{Vector{HeisenbergSpin}})
  
end 
    

test_compute_m2_af()
test_compute_T2_op(x_model(6))

num_spins = 100
test_ising_spins = [rand(Int8.([-1,1])) for i in 1:num_spins]
test_compute_1dcorr(test_ising_spins)
