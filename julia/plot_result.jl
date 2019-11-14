
using Plots
Plots.gr()

L = 2
const HeisenbergSpin = Tuple{Float64,Float64,Float64}

# make array having kagome's lattice point.
function mk_kagome(L::Int64)
    
    a1 = (2.,0.)
    a2 = (2*cos(π/3),2*sin(π/3))
    
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
function plot_spin(lattice::Array{Tuple{Float64,Float64}}, spins::AbstractArray{HeisenbergSpin})
    
    lattice_x = [lattice[i][1] for i in 1:length(lattice)]
    lattice_y = [lattice[i][2] for i in 1:length(lattice)]
    
    spins_x   = [spins[i][1] for i in 1:length(spins)]
    spins_y   = [spins[i][2] for i in 1:length(spins)]
    
    quiver(lattice_x,lattice_y, quiver=(spins_x,spins_y),title="spins' direction",framestyle=:box)
    
    savefig("SpinOnKagome")
end

test = [((cos(pi/2 + (i-1)*2pi/3), sin(pi/2 + (i-1)*2pi/3),0.)./3) for i in 1:3*L^2] 
plot_spin(mk_kagome(L), test)
