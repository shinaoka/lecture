#parameters for Jij
const d = 3
L = 10
num_spins = L^d
const Jx= 1.0
const Jy= 1.0
const Jz= 1.0

#parameters for temepratures.
num_temps = 48
min_T = 0.001
max_T = 0.005

function input_temperatures(num_temps::Int64,min_T::Float64,max_T::Float64)
    
    open("temperatures.txt", "w") do fp
       println(fp, num_temps)
       temp = range(min_T, stop=max_T, length=num_temps)
       for i in temp
          println(fp, i)
       end
     end
end

function input_Jij(couplings::Array{Any,1}, Jx::Float64, Jy::Float64, Jz::Float64)
    
    open("Jij.txt", "w") do fp
      println(fp, length(couplings))
      for pair in couplings
          println(fp,pair[1]," ",pair[2]," ",Jx," ",Jy," ", Jz)
      end
    end
end

# some following functions generate input_Jij's argument array of tuple.


