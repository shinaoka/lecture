num_insert  = 10
start_temps = 2e-2
end_temps   = 3e-2
input_file  = "temperatures.txt"
output_file = "temperatures.txt"

# Read a list of temperatures
function read_temps(temperature_file::String)
    temps = Vector{Float64}(undef, 0)
    num_temps = 0
    open(temperature_file) do file
        num_temps = parse(Int64, readline(file))
        for l in 1:num_temps
            temp = parse(Float64, readline(file))
            push!(temps, temp)
        end
    end

    # Check if temperatures are in ascending order or descending order
    if !all(temps[1:num_temps-1] .< temps[2:num_temps]) && !all(temps[1:num_temps-1] .> temps[2:num_temps])
        error("Temperatures must be given either in ascending order or in descending order!")
    end

    return temps
end


function insert_temps(num_insert,start_temps,end_temps,input_file,output_file)
    
    insert_temps = LinRange(start_temps,end_temps,num_insert)
    original_temps = read_temps(input_file)
    
    target_temps = union(original_temps, insert_temps)
    sort!(target_temps)

    open(output_file,"w") do fp
        num_temps    = length(target_temps)
        println(fp,num_temps)
        for idx in 1:num_temps
            println(fp,target_temps[idx])
        end
    end

end


insert_temps(num_insert,start_temps,end_temps,input_file,output_file)
