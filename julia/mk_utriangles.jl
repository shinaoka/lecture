using LinearAlgebra

# prameters for system.
L = 500

function mk_upward_triangles(L,file_name)

    kagome = Dict()
    idx = 1
    L = 2L

    for (i,j) in Iterators.product(1:L,1:L)
        mod_site = mod.((i,j),L)
        if mod(mod_site[1],2) == 0 && mod(mod_site[2],2) == 1
            continue
        end
        get!(kagome,mod_site,idx)
        idx += 1
    end

    upward_triangles = []
    
    for key in keys(kagome)
        
        site1 = mod.(key .+ (1,0),L)
        site2 = mod.(key .+ (1,1),L)
        if in(site1,keys(kagome)) && in(site2,keys(kagome))
            push!(upward_triangles,(kagome[key],kagome[site1],kagome[site2]))
        end
    
    end

    open(file_name,"w") do fp

        num_triangles = length(upward_triangles)
        println(fp,num_triangles)

        for i in 1:num_triangles
            i_triangle = upward_triangles[i]
            println(fp,i_triangle[1]," ",i_triangle[2]," ",i_triangle[3])
        end
    end

end


@time mk_upward_triangles(L,"utriangles.txt")
