using PyPlot
using HDF5

function plot_mz2(file_name::String)
    fp = h5open(file_name,"r")
    num_spins = read(fp,"num_spins") 
    temps     = read(fp,"temps"    ) 
    M2        = read(fp,"Mz2"      ) 

    m2 = (9*M2) / (num_spins^2)
        
    PyPlot.plot(temps,m2)
    close(fp)
end


plot_mz2("N48.h5")
#plot_mz2("N192.h5")
#plot_mz2("N768.h5")

savefig("M-T")
