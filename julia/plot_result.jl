using LinearAlgebra
using Plots
Plots.gr()

"""
# prameters for system.
L = 16
l = 1.
num_stack = 1
num_spins = (3*L^2)*num_stack
"""

# lattice vectors
Type3dVector = Tuple{Float64,Float64,Float64}
lat_vec1 = l .* (2.,0.,0.)
lat_vec2 = l .* (2*cos(π/3),2*sin(π/3),0.)
lat_vec3 = l .* (0.,0.,1.)

function mk_stacked_structure(L::Int64,num_stack::Int64,a1::Type3dVector,a2::Type3dVector,a3::Type3dVector)
    index = 1
    temp = fill((0.,0.,0.), num_spins)
    for layer in 1:num_stack
        for (i,j) in Iterators.product(1:L,1:L)
            A = (i-1).* a1 .+ (j-1).*a2 .+ (layer-1).*a3
            B = A .+ a1./2
            C = A .+ a2./2
            temp[index  ] = A
            temp[index+1] = B
            temp[index+2] = C
            index += 3
        end
    end
    return temp
end

# make function plotting x-y spin on any lattice.
function plot_spin_direction(lattice::Array{Tuple{Float64,Float64,Float64},1},spins_x::Array{Float64,1},spins_y::Array{Float64,1},num_spins)

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
test = [((cos(pi/2 + (i-1)*2pi/3), sin(pi/2 + (i-1)*2pi/3),0.)./3) for i in 1:num_spins] 
test_x = [test[i][1] for i in 1:length(test)]
test_y = [test[i][2] for i in 1:length(test)]

plot_spin_direction(mk_stacked_structure(L,num_stack,lat_vec1,lat_vec2,lat_vec3),test_x,test_y,num_spins)
savefig("test.png")

"""
