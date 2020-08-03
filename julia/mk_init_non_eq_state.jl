
function write_init_non_eq_state(L)
     
    num_spins = 3L^2
    open("init_non_eq_state.txt","w") do fp 

        println(fp,num_spins)
        for i in 1:num_spins
            println(fp, i, " ", cos(2i*pi/3), " ", sin(2i*pi/3), " " , 0.0)
        end
    end
end

write_init_non_eq_state(48)
