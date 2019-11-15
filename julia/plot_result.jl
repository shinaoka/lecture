
using Plots
Plots.gr()

#L = 2

# make array having kagome's lattice point.
function mk_kagome(L::Int64)
    
    a1 = (2.,0.)
    a2 = (2*cos(pi/3),2*sin(pi/3))
    
    index = 1
    temp = fill((0.,0.),3*L^2)
    for (i,j) in Iterators.product(1:L,1:L)
        A = (i-1).* a1 .+ (j-1).*a2
        B = A .+ a1./2
        C = A .+ a2./2
        temp[index] = A
        temp[index+1] = B
        temp[index+2] = C
        index += 3
    end
    
    return temp
end

# make function plotting x-y spin on any lattice.
function plot_spin_direction(lattice::Array{Tuple{Float64,Float64}},spins_x::Array{Float64,1},spins_y::Array{Float64,1},num_spins)
    
    lattice_x = [lattice[i][1] for i in 1:length(lattice)]
    lattice_y = [lattice[i][2] for i in 1:length(lattice)]
    
    temp_sx = zeros(num_spins)
    temp_sy = zeros(num_spins)
    for it in 1:num_spins
        temp_norm = norm((spins_x[it],spins_y[it]))
        temp_sx[it] = spins_x[it] / temp_norm    
        temp_sy[it] = spins_y[it] / temp_norm    
    end

    quiver(lattice_x,lattice_y, quiver=(temp_sx,temp_sy),title="spins' direction",framestyle=:box)

end

"""
# test code.
test = [((cos(pi/2 + (i-1)*2pi/3), sin(pi/2 + (i-1)*2pi/3),0.)./3) for i in 1:3*L^2] 
test_x = [test[i][1] for i in 1:length(test)]
test_y = [test[i][2] for i in 1:length(test)]

plot_spin_direction(mk_kagome(L),test_x,test_y,1)
"""
