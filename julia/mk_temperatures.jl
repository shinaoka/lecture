#parameters for temepratures.
num_temps = 10
min_T = 1e-2
max_T = 1e-1

function input_temperatures(num_temps::Int64,min_T::Float64,max_T::Float64)
    logT = LinRange(log(min_T), log(max_T), num_temps)
    open("temperatures.txt", "w") do fp
       println(fp, num_temps)
       for i in logT
          println(fp, exp(i))
       end
     end
end

# Output temperatures.
input_temperatures(num_temps,min_T,max_T)
