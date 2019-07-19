struct Accumulator
    count::Dict{String,UInt64}
    data::Dict{String,Any}
    num_temps::Integer
    function Accumulator(num_temps::Integer)
        new(Dict{String,UInt64}(), Dict{String,Any}(), num_temps)
    end
end

function add!(acc::Accumulator, name::String, data)
    @assert length(data) == acc.num_temps
    if !haskey(acc.count, name)
        acc.count[name] = 1
        acc.data[name] = copy(data)
    else
        acc.count[name] += 1
        acc.data[name] += data
    end
end

function mean(acc::Accumulator, name::String)
    return acc.data[name]/acc.count[name]
end
