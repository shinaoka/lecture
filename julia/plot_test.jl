using Plots

L = 1

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
lattice = mk_kagome(L)
println(lattice)
plot(lattice, seriestype =:scatter,title="Kagome")
