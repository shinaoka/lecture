L = 5000
J1 = -1.
J1 = (J1,J1,J1)
J2 = -0.02
J2 = (J2,J2,J2)
flag1 = 1
flag2 = 0

function mk_kagome_v2(L)
   
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
    
    return kagome
end


function mk_Jij_kagome_nn(L,J1_x,J1_y,J1_z,flag1)

    kagome = mk_kagome_v2(L)
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

    kagome = mk_kagome_v2(L)
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


