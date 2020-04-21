using LinearAlgebra

# prameters for system.
L = 4
num_spins = 3*L^2

# lattice vectors
Type3dVector = Tuple{Float64,Float64,Float64}

# length betweeen 1st&2nd nearest neightbor sites.
len1 = 2*cos(pi/3)
len2 = 2*sin(pi/3)

# parameter for Interaction.
SSInteraction = Type3dVector
tempJ1 = -1.
tempJ2 = -0.0
J1 = (tempJ1,tempJ1,tempJ1)
J2 = (tempJ2,tempJ2,tempJ2)

#parameters for temepratures.
num_temps = 1
min_T = 0.00001
max_T = min_T

function input_temperatures(num_temps::Int64,min_T::Float64,max_T::Float64)
    logT = LinRange(log(min_T), log(max_T), num_temps)
    open("temperatures.txt", "w") do fp
       println(fp, num_temps)
       for i in logT
          println(fp, exp(i))
       end
     end
end

# Output temperatures.
input_temperatures(num_temps,min_T,max_T)


# some following functions generate input_Jij's argument array of tuple.

function mk_kagome(L::Int64)
    
    a1 = (2.,0.)
    a2 = (2cos(pi/3),2sin(pi/3))
    index = 1
    kagome = Dict()
    for (i,j) in Iterators.product(1:L,1:L)
        A = (i-1) .* a1 .+ (j-1) .* a2
        B = A .+ a1./2
        C = A .+ a2./2
        get!(kagome,index  ,A)
        get!(kagome,index+1,B)
        get!(kagome,index+2,C)
        index += 3
    end

    return kagome

end

function compute_distance(L::Int64,p1,p2)
    
    a1 = (2.,0.)
    a2 = (2cos(pi/3),2sin(pi/3))

    candidate_distance = []
    for (i,j) in Iterators.product(1:3,1:3)
        tempi = L*(i-2)
        tempj = L*(j-2)

        lattice_vec = tempi .* a1 .+ tempj .* a2
        cp_p2 = p2 .+ lattice_vec
        push!(candidate_distance, norm(p1.-cp_p2))
    end

    return minimum(candidate_distance)
    
end

function mk_interaction(L::Int64,J1::SSInteraction,J2::SSInteraction,len1::Float64,len2::Float64)
    interaction = []
    kagome = mk_kagome(L)
    
    
    for (site1,site2) in Iterators.product(keys(kagome),keys(kagome))
        
        if site1 < site2
                
            distance = compute_distance(L,kagome[site1],kagome[site2])
    
            # make 1st nearest neighbor interaction.    
            # last element 1 is sign that this interaction is J1  
            if isapprox(len1,distance)
                push!(interaction, (site1,site2,J1[1],J1[2],J1[3],1))
            
            # make 2nd nearest neighbor interaction.    
            # last element 0 is sign that this interaction is J2  
            elseif isapprox(len2,distance)
                push!(interaction, (site1,site2,J2[1],J2[2],J2[3],0))
            
            end
        end
    end
    
    return sort!(interaction)

end

function input_Jij(interaction::Array{Any,1})

    open("Jij.txt", "w") do fp
        
        println(fp,length(interaction))
        
        for intr in interaction
            println(fp,intr[1]," ",intr[2]," ",intr[3]," ",intr[4]," ",intr[5]," ",intr[6])
        end
        
    end
    
end

# output Jij
println("num_spins: ",num_spins)
input_Jij(mk_interaction(L,J1,J2,len1,len2))
