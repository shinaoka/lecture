#parameters for system
L = 3
J1 = -1.
J2 = 0.05
J1 = (J1,J1,J1)
J2 = (J2,J2,J2)
flag1 = 1
flag2 = 0


function mk_kagome(L)
   
    kagome = Dict()
    idx = 1
    L = 2L    

    for (i,j) in Iterators.product(1:L,1:L)
        mod_site = mod.((i,j),L)
        if mod(mod_site[1],2) == 0 && mod(mod_site[2],2) == 1
            continue
        end
        kagome[mod_site] = copy(idx) 
        idx += 1
    end
    
    return kagome
end


function mk_kagome2(L,output_file)

    kagome = Dict()
    idx = 1
    L = 2L

    for (i,j) in Iterators.product(1:L,1:L)
        mod_site = mod.((i,j,0),L)
        if mod(mod_site[1],2) == 0 && mod(mod_site[2],2) == 1
            continue
        end

        kagome[idx] = copy.(mod_site)
        idx += 1
    end

    open(output_file,"w") do fp 
        num_sites = length(keys(kagome))
        println(fp,num_sites)
        for i in 1:num_sites
            site = kagome[i]
            println(fp,i," ",site[1]," ",site[2]," ",site[3])
        end
    end

    return kagome
end

mk_kagome2(L,"kagome.txt")


function mk_triangles(L)
    kagome = mk_kagome(L)
    L = 2L
    upward_triangles   = []
    downward_triangles = []
    
    for key in keys(kagome)
        
        usite1 = mod.(key .+ (1,0),L)
        usite2 = mod.(key .+ (1,1),L)
        dsite1 = mod.(key .+ (-1,0),L)
        dsite2 = mod.(key .+ (-1,-1),L)
        
        if in(usite1,keys(kagome)) && in(usite2,keys(kagome))
            push!(upward_triangles,[kagome[key],kagome[usite1],kagome[usite2]])
        end

        if in(dsite1,keys(kagome)) && in(dsite2,keys(kagome))
            push!(downward_triangles,[kagome[key],kagome[dsite1],kagome[dsite2]])
        end
    
    end

    @assert length(upward_triangles) == length(downward_triangles)

    return upward_triangles, downward_triangles

end


function write_triangles(triangles,output_file)

    num_triangles = length(triangles)
    open(output_file,"w") do fp
        num_triangles = length(triangles)
        println(fp,num_triangles)

        for i in triangles
            println(fp,i[1]," ",i[2]," ",i[3])
        end
    end

end

ud = mk_triangles(L)
u  = ud[1]
d  = ud[2]
@assert length(u) == length(d)
write_triangles(u,"utriangles.txt")
write_triangles(d,"dtriangles.txt")


function mk_q0(L)

    num_spins = 3L^2
    spins = fill((0.,0.,0),num_spins)
    
    utriangles = mk_triangles(L)[1]
    for i in utriangles
        for j in 1:3
            spins[i[j]] = (cos(2j*pi/3),sin(2j*pi/3),0.0)
        end
    end

    return spins

end


function assin_number!(triangles)

    num_triangles = length(triangles)
    smallest_site_indices = zeros(Int64,num_triangles)
    for i in 1:num_triangles
        smallest_site_indices[i] = minimum(triangles[i])
    end

    sort!(smallest_site_indices)
    #println("debug b:",smallest_site_indices)

    temps = Dict()
    idx = 1
    for isite in smallest_site_indices
        temps[isite] = copy(idx)
        idx += 1
    end

    #println("debug c:",temps)

    for i in 1:num_triangles
        for j in 1:3
            site_ij = triangles[i][j]
            if in(site_ij,keys(temps))
                push!(triangles[i],temps[site_ij])
            end
        end
    end

    #println("debug d:",triangles)
end

function classify(L,triangles)
    num_triangles = length(triangles)
    #crassify each numberd upward triangle row to three groups.
    g1 = []
    g2 = []
    g3 = []
    for i in 1:num_triangles
        row = div(triangles[i][4],L) + 1
        if mod(triangles[i][4],L) == 0
            row -= 1
        end
        #println("debug f:",triangles[i]," ",column)

        mod_row = mod(row,3)
        if     mod_row == 0
            push!(g1,triangles[i])
        elseif mod_row == 1
            push!(g2,triangles[i])
        else   mod_row == 2
            push!(g3,triangles[i])
        end
    end

    return g1,g2,g3
end


function sort_triangle_in_group!(g)
    @assert mod(length(g),3) == 0
    reverse!.(g)
    sort!(g)
    reverse!.(g)
end

function put_spin_config!(spins,g,theta)
    count = 1
    for i in 1:length(g)
        for j in 1:3
            arg = theta + 2j*pi/3 + mod(count,3)*4pi/3
            spins[g[i][j]] = (cos(arg),sin(arg),0)
        end
        count += 1
    end

end


function mk_sqrt3(L)

    triangles = mk_triangles(L)[1]
    assin_number!(triangles)
    gs = classify(L,triangles)
    sort_triangle_in_group!.(gs)
    num_spins = 3L^2
    spins = fill((0.,0.,0.),num_spins)
    for i in 1:3
        theta = 4pi*i/3 + pi/2
        put_spin_config!(spins,gs[i],theta)
    end

    return spins
end


function write_spins(spins,file_name)
    open(file_name,"w") do fp
        num_spins = length(spins)
        println(fp,num_spins)
        for s in spins
            println(fp,s[1]," ",s[2]," ",s[3])
        end
    end
end

# write initial spin configs q0 and sqrt3.
q0 = mk_q0(L)
write_spins(q0,"q0.txt")
sqrt3 = mk_sqrt3(L)
write_spins(sqrt3,"sqrt3.txt")


function mk_Jij_kagome_nn(L,J1_x,J1_y,J1_z,flag1)

    kagome = mk_kagome(L)
    Jij = []
   
    L = 2L
    for key in keys(kagome) 

        for (dx,dy) in Iterators.product(0:1,0:1)
            if dx == dy == 0
                continue
            end
            
            nn_mod_site = mod.((key .+ (dx,dy)),L)

            if in(nn_mod_site, keys(kagome))
                site_i = kagome[key]
                site_j = kagome[nn_mod_site]
                
                if site_i > site_j
                    site_i,site_j = site_j,site_i 
                end

                push!(Jij,(site_i,site_j,J1_x,J1_y,J1_z,flag1))
            end
        end
    end

    return sort(Jij)
end
    

function mk_Jij_kagome_nnn(L,J1,flag1,J2,flag2)

    kagome = mk_kagome(L)
    Jij = mk_Jij_kagome_nn(L,J1[1],J1[2],J1[3],flag1)

    L = 2L
    
    for key in keys(kagome)

        for (dx,dy) in Iterators.product(-2:2,-2:2)
      
            if (dx,dy) == (1,2) || (dx,dy) == (2,1) || (dx,dy) == (1,-1)

                nn_mod_site = mod.((key .+ (dx,dy)),L)

                if in(nn_mod_site, keys(kagome))
                    site_i = kagome[key]
                    site_j = kagome[nn_mod_site]

                    if site_i > site_j
                        site_i,site_j = site_j,site_i
                    end
                    #println(site_i," ",site_j)

                    push!(Jij,(site_i,site_j,J2[1],J2[2],J2[3],flag2))
                end
            end
        end
    end
    return sort(Jij)
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

Jij_kagome_nnn = mk_Jij_kagome_nnn(L,J1,flag1,J2,flag2)
@time input_Jij(3L^2,Jij_kagome_nnn)


