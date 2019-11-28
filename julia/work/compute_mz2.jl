
using HDF5

h5open("spin_config.h5","r") do fp
    num_spins = read(fp,"num_spins")
    temp = read(fp,"temp")
    sz = read(fp,"sz")
    mz2 = (sum(sz))^2/ (num_spins^2)
    
    println("mz2: ",mz2) 
end


