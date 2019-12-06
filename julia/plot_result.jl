using PyPlot
using HDF5
using LaTeXStrings

function plot_mz2(file_name::String,title,xlabel,ylavel,label,color::String)
    fp = h5open(file_name,"r")
    num_spins = read(fp,"num_spins") 
    temps     = read(fp,"temps"    ) 
    M2        = read(fp,"Mz2"      ) 

    m2 = (9*M2) / (num_spins^2)
        
    PyPlot.plot(temps,m2,color=color,label=label)
    PyPlot.title(title)
    PyPlot.xlabel(xlabel)
    PyPlot.ylabel(ylabel)
    PyPlot.legend(loc="upper right")
    close(fp)
end

title  = "temperature and siize dependence of $(L"\langle M_z^2 \rangle")"
xlabel = L"T[\rm K]"
ylabel = L"\langle M_z^2 \rangle /N_{\rm unit}"

#plot_mz2("N48.h5"  ,title,xlabel,ylabel,L"L=4" ,"red"   )
plot_mz2("N192.h5" ,title,xlabel,ylabel,L"L=8" ,"green" )
plot_mz2("N768.h5" ,title,xlabel,ylabel,L"L=16","blue"  )
plot_mz2("N3072.h5",title,xlabel,ylabel,L"L=32","yellow")

savefig("M-T")

function plot_c(file_name::String)

    fp = h5open(file_name,"r")
    num_spins = read(fp,"num_spins")
    temps     = read(fp,"temps"    )
    E         = read(fp,"E"        )
    E2        = read(fp,"E2"       )

    heat_capa = zeros(Float64,length(temps))
  
    for i in 1:length(heat_capa)
        heat_capa[i] = ((E2[i]  - E[i]^2) / (temps[i]^2)) / num_spins
    end
   
    PyPlot.plot(temps,heat_capa)
    close(fp)
end
#plot_c("N48.h5")
#savefig("c-T.png")

function plot_order_parameter(file_name::String)

    fp = h5open(file_name,"r")
    num_spins = read(fp,"num_spins")
    temps     = read(fp,"temps"    )
    M1        = read(fp,"M1"       )
    M2        = read(fp,"M2"       )
    M3        = read(fp,"M3"       )

    temp = num_spins/3
    PyPlot.plot(temps,M1/(temp^2))
    PyPlot.plot(temps,M2/(temp^2))
    PyPlot.plot(temps,M3/(temp^2))
    close(fp)
end

#plot_order_parameter("N48.h5")
#savefig("order_parameter")
