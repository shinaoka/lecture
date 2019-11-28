using PyPlot
using HDF5

function temp_plot(file_name::String)
    h5open(file_name,"r") do fp
        temps = read(fp,"temps")
        num_spins = read(fp,"num_spins")
        m2 = read(fp,"Mz2") / (num_spins^2)

        println(m2)
 
        PyPlot.plot(temps,m2)
      
        savefig("M-T.png")
    end
end

temp_plot("phys_quant.h5")
