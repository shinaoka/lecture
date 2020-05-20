using Test

include("measure_mc.jl")

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


function test_m2_af()
    
    num_spins = 6
    spins = x_model(num_spins) 
   
    triangles = [(1,2,3),(4,5,6)]
    m2_af = compute_m2_af(spins[1],num_spins,triangles)
    @test isapprox(m2_af,1.0)

end


function test_T2_op(spins::Vector{Vector{HeisenbergSpin}})
    
    num_spins = length(spins[1])
    T2_op = compute_T2_op(spins[1],num_spins)
    @test isapprox(T2_op,octopolar_v2(spins,num_spins,1)[1])

end

test_m2_af()
test_T2_op(x_model(6))

