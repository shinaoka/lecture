
function write_init_non_eq_state(L)
     
    num_spins = 3L^2
    open("init_non_eq_state.txt","w") do fp 

        println(fp,num_spins)
        for i in 1:num_spins
        
            """
            theta = 10*rand()
            phi   = 10*rand()
            site  = (sin(theta)cos(phi),sin(theta)sin(phi),cos(theta))
            """
            site = (1.,0.,0.)
            println(fp,i," ",site[1]," ",site[2]," ",site[3])
        end
    end
end

write_init_non_eq_state(16)
