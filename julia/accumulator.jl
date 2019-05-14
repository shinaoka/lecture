struct Accumulator
    count::Dict{String,UInt64}
    data::Dict{String,Any}
end

function measure(acc:Accumulator, String::name, data)
    if !haskey(acc.count, name)
        acc.count[name] = 1
        acc.data[name] = copy(data)
    else
        acc.count[name] += 1
        acc.data[name] += data
    end
end

function mean(acc:Accumulator, String::name)
    return acc.data[name]/acc.count[name]
end

