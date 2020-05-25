using LinearAlgebra

# prameters for system.
L = 192

function mk_upward_triangles(file_name::String,L::Int64)
    
    open(file_name,"w") do fp 
        println(fp,L^2)
        idx = 1
        for i in 1:L^2
            println(fp,idx," ",idx+1," ",idx+2)
            idx += 3
        end
    end
  
end

ts = time_ns()
mk_upward_triangles("utriangles.txt",L)
println("elapsed time $(time_ns()- ts))")
