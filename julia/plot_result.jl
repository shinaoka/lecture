using LinearAlgebra
using HDF5
using PyPlot

# prameters for system.
num_stack = 1
num_spins = h5open("12_0.1.h5","r") do fp
                read(fp,"num_spins")
end   

L = Int(sqrt(num_spins/3))

# lattice vectors
Type3dVector = Tuple{Float64,Float64,Float64}
l = 1.
lat_vec1 = l .* (2.,0.,0.)
lat_vec2 = l .* (2*cos(pi/3),2*sin(pi/3),0.)
lat_vec3 = l .* (0.,0.,1.)



function mk_stacked_structure(L::Int64,num_stack::Int64,a1::Type3dVector,a2::Type3dVector,a3::Type3dVector,num_spins::Int64)
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

function plot_spin_direction(lattice::Array{Tuple{Float64,Float64,Float64},1},spins_x::Array{Float64,1},spins_y::Array{Float64,1},num_spins::Int64)
    c = "red"
    lw = 0.5
    ls = :dash
    
    lattx = [lattice[i][1] for i in 1:length(lattice)]
    latty = [lattice[i][2] for i in 1:length(lattice)]

    temp_sx = zeros(num_spins)
    temp_sy = zeros(num_spins)
    for it in 1:num_spins
        temp_norm = norm((spins_x[it],spins_y[it]))
        temp_sx[it] = spins_x[it] / temp_norm
        temp_sy[it] = spins_y[it] / temp_norm
    end
    
    PyPlot.quiver(lattx,latty,temp_sx,temp_sy,pivot=:middle)
    
    for i in 1:length(lattx)
        for j in 1:length(latty)
            if i < j && abs(1-norm(lattice[i].-lattice[j])) < 1e-5
            
                PyPlot.plot([lattx[i],lattx[j]],[latty[i],latty[j]],color=c)
                
            end
        end
    end
    
end

# pull latest spin data from .h5 file.
sx = h5open("12_0.1.h5","r") do fp
         read(fp,"sx")
end

sy = h5open("12_0.1.h5","r") do fp
         read(fp,"sy")
end

sz = h5open("12_0.1.h5","r") do fp
         read(fp,"sz")
end

"""
println("sx: ",sx)
println("sy: ",sy)
println("sz: ",sz)
"""

plot_spin_direction(mk_stacked_structure(L,num_stack,lat_vec1,lat_vec2,lat_vec3,num_spins),sx,sy,num_spins)
savefig("12_0.1.png")





