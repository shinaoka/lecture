using LinearAlgebra

# prameters for system.
L = 9
l = 1.
num_stack = 1
num_spins = (3*L^2)*num_stack

# lattice vectors
Type3dVector = Tuple{Float64,Float64,Float64}
lat_vec1 = l .* (2.,0.,0.)
lat_vec2 = l .* (2*cos(pi/3),2*sin(pi/3),0.)
lat_vec3 = l .* (0.,0.,1.)

# length betweeen 1st&2nd nearest neightbor sites.
len1 = 2*cos(pi/3)
len2 = 2*sin(pi/3)

# parameter for Interaction.
SSInteraction = Type3dVector
J_1stnn  = (-1.,-1.,-2.)
J_2ndnn  = (-0.005,-0.005,-0.005)
J_intlay = (1.,1.,1.)

#parameters for temepratures.
num_temps = 24
min_T = 0.2
max_T = 0.25

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


# some following functions generate input_Jij's argument array of tuple.

#make stacked structure for any given lattice vector.
function mk_stacked_structure(L::Int64,num_stack::Int64,a1::Type3dVector,a2::Type3dVector,a3::Type3dVector)
    index = 1
    temp = Dict()
    for layer in 1:num_stack
        for (i,j) in Iterators.product(1:L,1:L)
            A = (i-1).* a1 .+ (j-1).*a2 .+ (layer-1).*a3
            B = A .+ a1./2
            C = A .+ a2./2
            get!(temp,index  ,A)
            get!(temp,index+1,B)
            get!(temp,index+2,C)
            index += 3
        end
    end
    return temp
end

function compute_distance(L::Int64,n::Int64,p1::Type3dVector,p2::Type3dVector,a1::Type3dVector,a2::Type3dVector,a3::Type3dVector)

    temp = []
    for (i,j,k) in Iterators.product(1:3,1:3,1:3)
        tempi = L*(i-2)
        tempj = L*(j-2)
        tempk = n*(k-2)

        lattice_vec = tempi.*a1 .+ tempj.*a2 .+ tempk.*a3
        cp_point = p2 .+ lattice_vec
        push!(temp, norm(p1.-cp_point))
    end

    return minimum(temp)
    
end

function mk_interaction(L::Int64,num_stack::Int64,a1::Type3dVector,a2::Type3dVector,a3::Type3dVector,J1::SSInteraction,J2::SSInteraction,J3::SSInteraction,len1::Float64,len2::Float64)
    interaction = []
    crystal = mk_stacked_structure(L,num_stack,a1,a2,a3)
    
    
    for key1 in keys(crystal)
        for key2 in keys(crystal)
            if key1 < key2
                
                distance = compute_distance(L,num_stack,crystal[key1],crystal[key2],a1,a2,a3)
                z1 = crystal[key1][3]
                z2 = crystal[key2][3]

                # make x-y plane interaction.
                if z1 == z2
                    
                    # make 1st nearest neighbor interaction.
                    if abs(len1-distance) < 1e-10
                        push!(interaction, (key1,key2,J1[1],J1[2],J1[3]))
                       
                    
                    # make 2nd nearest neighbor interaction.    
                    elseif abs(len2-distance) < 1e-10 
                        push!(interaction, (key1,key2,J2[1],J2[2],J2[3]))
                    

                    end     
                    
                # make inter layer interaction.
                else
                    if abs(len1-distance) < 1e-10
                        push!(interaction, (key1,key2,J3[1],J3[2],J3[3]))
                    end
                    
                end
            end
        end
    end
    
    return sort!(interaction)
end

function input_Jij(interaction::Array{Any,1})

    open("Jij.txt", "w") do fp
        
      println(fp,length(interaction))
        
      for intr in interaction
          println(fp,intr[1]," ",intr[2]," ",intr[3]," ",intr[4]," ", intr[5])
      end
        
    end
    
end

# output Jij
println("num_spins: ",num_spins)
input_Jij(mk_interaction(L,num_stack,lat_vec1,lat_vec2,lat_vec3,J_1stnn,J_2ndnn,J_intlay,len1,len2))
