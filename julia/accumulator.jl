using MPI
using Test
using EllipsisNotation

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

function mean_gather(acc::Accumulator, name::String, comm)
    counts = convert(Array{Cint}, MPI.Allgather(acc.num_temps, comm))
    results_local = mean(acc, name)
    return MPI.Gatherv(results_local, counts, 0, comm)
end

function mean_gather_array(acc::Accumulator, name::String, comm)
    rank = MPI.Comm_rank(comm)
    num_proc = MPI.Comm_size(comm)

    results_local = mean(acc, name)
    num_temps_local = length(results_local)
    num_temps = MPI.Allreduce(num_temps_local, MPI.SUM, comm)
    data_size = size(results_local[1])
    data_type = typeof(results_local[1])

    for r in results_local
        @test size(r) == data_size
        @test typeof(r) == data_type
    end

    # Flatten data
    data_local_flatten = collect(Iterators.flatten(results_local))

    counts = convert(Array{Cint}, MPI.Allgather(length(data_local_flatten), comm))
    data_flatten = MPI.Gatherv(data_local_flatten, counts, 0, comm)
    if rank > 0
        return
    end

    right_dim = (num_temps,)
    data = reshape(data_flatten,  (data_size..., right_dim...))

    result = Array{data_type}(undef, num_temps)
    for it in 1:num_temps
        result[it] = data[.., it]
    end
    return result
end

