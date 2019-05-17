struct Accumulator
    count::Dict{String,UInt64}
    data::Dict{String,Any}
    function Accumulator()
        new(Dict{String,UInt64}(), Dict{String,Any}())
    end
end

function add!(acc::Accumulator, name::String, data)
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
