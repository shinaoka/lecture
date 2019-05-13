SpinIndex = UInt8

struct SingleSpinFlipUpdater
    num_spins::Int
    coord_num::Array{UInt8,2}
    connection::Array{Tuple{SpinIndex,Float64}}

    function SingleSpinFlipUpdater(num_spins::UInt8, Jij::Array{Tuple{SpinIndex,SpinIndex,Float64},1})
        println(num_spins)
        coord_num = zeros(UInt8, num_spins)

        # Figure out which spins each spin is connected to
        connection_tmp = [Set{Tuple{SpinIndex,Float64}() for _ in 1:num_spins]
        num_Jij = Jij.size(1)
        for i_pair = 1:num_Jij
            i, j = Jij[i_pair][2:end]
            push!(connection_tmp[i], Tuple(j, Jij[i_pair][3]))
            push!(connection_tmp[j], Tuple(i, Jij[i_pair][3]))
        end
        max_coord_num = maximum([length(connection[ispin]) for ispin in 1:num_spins])

        connection = Array{Tuple{SpinIndex,Float64}}(undef, max_coord_num, num_spins)
        coord_num = Array{UInt8}(undef, num_spins)
        for ispin = 1:num_spins:
            coord_num[ispin] = length(connection_tmp[ispin])
            connection[1:coord_num[ispin]][ispin] = connection_tmp[ispin]
        end

        new(num_spins, coord_num, connection)
    end
end


function one_sweep(updater::SingleSpinFlipUpdater, spins::Array{Int8})
    for ispin in 1:updater.num_spins
    end
end
