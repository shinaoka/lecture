# test
using BenchmarkTools
using Test

include("accumulator.jl") 

num_spins = 3
num_temps = 4

test1 = [fill((1.,0.,0.),num_spins) for i in 1:num_temps] 

function compute_magnetization(acc::Accumulator,num_spins::Int64,spins::Array{Array{Tuple{Float64,Float64,Float64},1},1},num_temps::Int64)

    mx = zeros(Float64, num_temps)
    my = zeros(Float64, num_temps)
    mz = zeros(Float64, num_temps)
    m2 = zeros(Float64, num_temps)

    for i in 1:num_temps
      for j in 1:num_spins
        mx[i] += spins[i][j][1]
        my[i] += spins[i][j][2]
        mz[i] += spins[i][j][3]
      end
      m2[i] = mx[i]^2 + my[i]^2 + mz[i]^2
    end
    add!(acc,"mx",mx)    
end

acc = Accumulator(num_temps)

compute_magnetization(acc, num_spins, test1, num_temps)

println(acc)

compute_magnetization(acc, num_spins, test1, num_temps)

println(acc)

println(mean(acc, "mx"))
