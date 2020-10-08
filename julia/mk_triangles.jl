using LinearAlgebra

# prameters for system.
L = 3

function mk_triangles(L,output_file1,output_file2)

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

    upward_triangles   = []
    downward_triangles = []
    
    for key in keys(kagome)
        
        usite1 = mod.(key .+ (1,0),L)
        usite2 = mod.(key .+ (1,1),L)
        dsite1 = mod.(key .+ (-1,0),L)
        dsite2 = mod.(key .+ (-1,-1),L)
        
        if in(usite1,keys(kagome)) && in(usite2,keys(kagome))
            push!(upward_triangles,(kagome[key],kagome[usite1],kagome[usite2]))
        end

        if in(dsite1,keys(kagome)) && in(dsite2,keys(kagome))
            push!(downward_triangles,(kagome[key],kagome[dsite1],kagome[dsite2]))
        end
    
    end

    @assert length(upward_triangles) == length(downward_triangles)

    open(output_file1,"w") do fp

        num_triangles = length(upward_triangles)
        println(fp,num_triangles)

        for i in upward_triangles
            println(fp,i[1]," ",i[2]," ",i[3])
        end
    end

    open(output_file2,"w") do fp

        num_triangles = length(downward_triangles)
        println(fp,num_triangles)

        for i in downward_triangles
            println(fp,i[1]," ",i[2]," ",i[3])
        end
    end


end


@time mk_triangles(L,"utriangles.txt","dtriangles.txt")
