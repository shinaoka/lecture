using HDF5
using PyPlot

function plot_sis1(file_name)

    fid = h5open(file_name,"r")
    sis1_lowestT = fid["sis1/1th_temps"]
    for i in 1:length(sis1_lowestT)
        ss = sis1_lowestT[i]
        println(i," ",ss)
    end

end

plot_sis1("sis1.h5")

