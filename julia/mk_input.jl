using LinearAlgebra

#parameters for Jij
L = 4
const Jx= 1.0
const Jy= 1.0
const Jz= 1.0

#parameters for temepratures.
num_temps = 48
min_T = 0.1
max_T = 1.5

function input_temperatures(num_temps::Int64,min_T::Float64,max_T::Float64)
    
    open("temperatures.txt", "w") do fp
       println(fp, num_temps)
       temp = range(min_T, stop=max_T, length=num_temps)
       for i in temp
          println(fp, i)
       end
     end
end

# Output temperatures.
input_temperatures(num_temps,min_T,max_T)

function input_Jij(couplings::Array{Any,1}, Jx::Float64, Jy::Float64, Jz::Float64)
    
    open("Jij.txt", "w") do fp
      println(fp, length(couplings))
      for pair in couplings
          println(fp,pair[1]," ",pair[2]," ",Jx," ",Jy," ", Jz)
      end
    end
end

# some following functions generate input_Jij's argument array of tuple.

# make dictionary has as keys kagome' site index corresponding to coordination.
function mk_kagome(L::Int64)
    
    a1 = (2.,0.)
    a2 = (2*cos(pi/3),2*sin(pi/3))
    
    index = 1
    temp = Dict()
    for (i,j) in Iterators.product(1:L,1:L)
        A = (i-1).* a1 .+ (j-1).*a2
        B = A .+ a1./2
        C = A .+ a2./2
        get!(temp,index  ,A)
        get!(temp,index+1,B)
        get!(temp,index+2,C)
        index += 3
    end
    return temp
end

function compute_distance(L::Int64,p1::Tuple{Float64,Float64},p2::Tuple{Float64,Float64})
    
    # base lattice vector 
    a1 = (2.,0.)
    a2 = (2*cos(pi/3),2*sin(pi/3))
    
    temp = []
    for (i,j) in Iterators.product(1:3, 1:3)
        tempi = L*(i-2)
        tempj = L*(j-2)
    
        lattice_vec = tempi.*a1 .+ tempj.*a2
        cp_point = p2 .+ lattice_vec
        push!(temp, norm(p1.-cp_point))
    end

    return minimum(temp)
end

function mk_coupling(L::Int64)
    nn_coupling = []
    kagome = mk_kagome(L)
    
    for key1 in keys(kagome)
        for key2 in keys(kagome)
            
            distance = compute_distance(L,kagome[key1],kagome[key2])
            
            if key1 < key2 && abs(1-distance) < 1e-5
                push!(nn_coupling, (key1,key2))
            end
        end
    end
    return sort!(nn_coupling)
end

# output Jij
input_Jij(mk_coupling(L),Jx,Jy,Jz)
