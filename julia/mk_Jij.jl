using LinearAlgebra

# prameters for system.
L = 192
num_spins = 3*L^2

# lattice vectors
Type3dVector = Tuple{Float64,Float64,Float64}

# length betweeen 1st&2nd nearest neightbor sites.
len1 = 2*cos(pi/3)
len2 = 2*sin(pi/3)

# parameter for Interaction.
SSInteraction = Type3dVector
tempJ1 = -1.
tempJ2 = -0.02
J1 = (tempJ1,tempJ1,tempJ1)
J2 = (tempJ2,tempJ2,tempJ2)


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
        lattice_vec = L * (i-2) .* a1 .+ L * (j-2) .* a2
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

function input_Jij(num_spins::Int64,interaction::Array{Any,1})

    open("Jij.txt", "w") do fp

        println(fp,num_spins)
        println(fp,length(interaction))
        
        for intr in interaction
            println(fp,intr[1]," ",intr[2]," ",intr[3]," ",intr[4]," ",intr[5]," ",intr[6])
        end
        
    end
    
end

# output Jij
println("num_spins: ",num_spins)
ts = time_ns()
input_Jij(num_spins,mk_interaction(L,J1,J2,len1,len2))
println("elapsed time $(time_ns()- ts)")