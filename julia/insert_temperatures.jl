num_insert = 100
begin_idx  = 1
end_idx    = 2
file_name  = "temperatures.txt"

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


function insert_temps(num_insert,begin_insert_idx,end_insert_idx,file_name)
    
    temps = read_temps(file_name)
    num_temps = length(temps)

    @assert 1 <= begin_insert_idx < end_insert_idx <= num_temps

    inserted_temps= Vector{Float64}(undef, 0)

    for idx in 1:begin_insert_idx
        push!(inserted_temps,temps[idx])
    end
    
    delta_T = (temps[end_insert_idx] - temps[begin_insert_idx]) / num_insert

    for insert_idx in 1:num_insert
        temp = temps[begin_insert_idx] + insert_idx*delta_T
        push!(inserted_temps,temp)
    end

    for idx in end_insert_idx+1:num_temps
        push!(inserted_temps,temps[idx])
    end

    num_temps = length(inserted_temps)
    open(file_name,"w") do fp
        println(fp,num_temps)
        for idx in 1:num_temps
            println(fp,inserted_temps[idx])
        end
    end

end

insert_temps(num_insert,begin_idx,end_idx,file_name)
