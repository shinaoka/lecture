
L = 3

function write_init_non_eq_state(L,u_triangle_file,output_file)

    num_spins = 3L^2
    spins = fill((0.,0.,0),num_spins)
    
    open(u_triangle_file,"r") do fp
        @assert L^2 == parse(Int64,readline(fp))
       
        for l in 1:L^2
            str = split(readline(fp))
            for i in 1:3
                site = parse(Int64, str[i])
                spins[site] = (cos(2i*pi/3),sin(2i*pi/3),0.0)
            end
            
        end
    end
 
    open(output_file,"w") do fp
        
        println(fp,num_spins)
        for i in 1:num_spins
            spin = spins[i]
            println(fp,spin[1]," ",spin[2]," ",spin[3])    
        end
    end

end


write_init_non_eq_state(L,"utriangles.txt","init_non_eq_state.txt")
