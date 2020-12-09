using HDF5
using PyPlot

function save_plot_corr(h5file,names)

    fid = h5open(h5file,"r")
    temps = fid["temperatures"]
    println("temperatures=$(temps)")
    num_temps = length(temps)
    println("num_temps=$(num_temps)")

    for name in names
        data = fid["$(name)/mean"]
        name = split(name,"_corr")[1]
        plt.figure()
        for it in 1:num_temps
            plt.plot(data[:,it],label="T=$(temps[it])")
            #plt.plot(data[:,it],label="T=$(temps[it])")
        end
        plt.title("correlation function of $(name)")
        plt.xlabel("MC step")
        plt.ylabel("$(name)")
        plt.yscale("log")
        plt.legend(loc="upper right")
        savefig("$(name)")
    end
    close(fid)

end

names = ["Gt","mq_sqrt3_corr", "mq_q0_corr", "afvc_corr", "fvc_corr","m_120degs_corr"]
#save_plot_corr("2d_out.h5",names)

function save_ave_corr(h5file,num_split,names,num_mcsteps,num_temps)

    num_names = length(names)
    aves = zeros(Float64,num_mcsteps,num_temps,num_names)

    for isplit in 1:num_split
        ith_h5file = replace(h5file,"x"=>"$(isplit)")
        fid = h5open(ith_h5file,"r")
        @assert num_temps == length(fid["temperatures"])

        for iname in 1:num_names
            data = fid["$(names[iname])/mean"]
            @assert num_mcsteps == size(data)[1]
            for it in 1:num_temps
                aves[:,it,iname] .+= data[:,it]
            end
        end
    end

    aves ./= num_split
    println("size of datas: $(size(aves))")


end
            

names = ["Gt","mq_sqrt3_corr", "mq_q0_corr", "afvc_corr", "fvc_corr","m_120degs_corr"]
h5file = "2d_splitx_out.h5"
num_split = 2
num_mcsteps = 201
num_temps = 1

save_ave_corr(h5file,num_split,names,num_mcsteps,num_temps)
    

