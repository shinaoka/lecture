#parameters for Jij
const d = 3
L = 10
num_spins = L^d
const J = 1.0

#parameters for temepratures.
num_temps = 24
min_T = 1.0
max_T = 1.5

function input_temperatures(num_temps::Int64,min_T::Float64,max_T::Float64)
    
    open("temperatures.txt", "w") do fp
       println(fp, num_temps)
       temp = range(min_T, stop=max_T, length=num_temps)
       for i in temp
          println(fp, i)
       end
     end
end

function input_Jij(couplings::Array{Any,1}, J::Float64)
    
    open("Jij.txt", "w") do fp
      println(fp, length(couplings))
      for pair in couplings
          println(fp,pair[1]," ",pair[2]," ",J)
      end
    end
end

# some following functions generate input_Jij's argument array of tuple.
# assign cartesian coordinate to site index.
function coord2ind(L::Int64, coord::Tuple{Int64,Int64,Int64})

    @assert (1,1,1) <= coord <= (L,L,L) "site1 is (1,1,1),given coordinate is out of cubic."
    index = 1+(coord[1]-1)+(coord[2]-1)*L+(coord[3]-1)*(L^2)
    return Int(index)
end

# transform cartesian to modulo coordinate of L and store modulo-index pairs.
function storemodL_index(d::Int64,L::Int64)

    temp_d = Dict()
    for i in Iterators.product(ntuple(i->1:L, d)...)
        modL = i .% L
        index = coord2ind(L, i)
        get!(temp_d, modL, index)
    end
    return temp_d
end

# find nearest-neighbor site couplings and store them. 
function nn_coupling(d::Int64, L::Int64,temp_d::Dict{Any,Any})

    coupling = []
    for i in Iterators.product(ntuple(i->1:L, d)...)
        index = coord2ind(L, i)
        modL = i .% L
        
        nn_x = (modL .+ (1,0,0)) .% L
        nn_y = (modL .+ (0,1,0)) .% L
        nn_z = (modL .+ (0,0,1)) .% L
        
        push!(coupling, (index, temp_d[nn_x]))
        push!(coupling, (index, temp_d[nn_y]))
        push!(coupling, (index, temp_d[nn_z]))
        
    end
    return coupling
end

# swap site index for input rigth upper Jij.
function swapindex(array)

    temp = []
    for itr in 1:length(array)
        i,j = array[itr][1], array[itr][2]
        if i > j
            push!(temp, (j,i))
        else
            push!(temp, (i,j))
        end
    end
    return sort!(temp)
end


#execute!
input_temperatures(num_temps,min_T,max_T)

modL_index   = storemodL_index(d,L)
nn_couplings = nn_coupling(d, L, modL_index)
input_Jij(swapindex(nn_couplings), J)
