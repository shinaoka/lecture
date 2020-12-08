using HDF5
using PyPlot

function save_plot_corr(h5file,names)

    fid = h5open(h5file,"r")
    temps = fid["temperatures"]
    num_temps = length(temps)
    for name in names
        data = fid["$(name)/mean"]
        name = split(name,"_corr")[1]
        plt.figure()
        for it in 1:num_temps
            plt.plot(data[:,it],label="T=$(it)")
        end
        plt.title("correlation function of $(name)")
        plt.xlabel("MC step")
        plt.ylabel("$(name)")
        plt.yscale("log")
        savefig("$(name)")

    end

end

names = ["Gt","mq_sqrt3_corr", "mq_q0_corr", "afvc_corr", "fvc_corr","m_120degs_corr"]
save_plot_corr("2d_out.h5",names)

