using HDF5
using PyPlot
using LsqFit


names = ["Gt","mq_sqrt3_corr", "mq_q0_corr", "afvc_corr", "fvc_corr","m_120degs_corr"]


function save_plot_corr(h5file,names)

    fid = h5open(h5file,"r")
    temps = fid["temperatures"]
    num_temps = length(temps)

    for name in names
        data = fid["$(name)/mean"]
        name = split(name,"_corr")[1]
        plt.figure()
        for it in 1:num_temps
            plt.plot(data[:,it],label="T=$(temps[it])")
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

#save_plot_corr("2d_out.h5",names)


function compute_τ(tdatas,ydatas,p0)
    model(t,p) = p[1] .- (1/p[2])*t
    fit = curve_fit(model, tdatas, ydatas, p0)
    fit.param[2]
end


function compute_Tc(temps,τs,p0)
    model(t,p) = p[1] .+ p[2] ./ sqrt.(t .- p[3])
    fit = curve_fit(model,temps,log.(τs),p0)
    fit.param[3]
end


function execute(fid,name)

    temps = fid["temperatures"]
    num_temps = length(temps)
    @assert num_temps > 1
    ydata = fid["$(name)/mean"]
    tdata = [i for i in 1:length(ydata[:,1])]
    p01 = [1.0,1.0]
    τs = Vector{Float64}(undef,num_temps)
    println("τs=$(τs)")
    
    for it in 1:num_temps
        temp_ydata = log.(abs.(ydata[:,it]))
        τs[it] = compute_τ(tdata,temp_ydata,p01)
    end
    
    p02 = [1.0,1.0,1.0]
    Tc = compute_Tc(temps,τs,p02)
    Tc
    
end

h5file = "2d_out.h5"
fid = h5open(h5file,"r")
name = "Gt"
Tc_of_Gt = execute(fid,name)
#println("Tc of Gt = $(Tc_of_Gt)")

function execute_all(h5file,names)

    fid = h5open(h5file,"r")
    Tcs = Dict()
    for name in names
        Tcs["$(name)"] = execute(fid,name)
    end

    Tcs
end

h5file = "2d_out.h5"
names = ["Gt","mq_sqrt3_corr", "mq_q0_corr", "afvc_corr", "fvc_corr","m_120degs_corr"]

#=
Tcs = execute_all(h5file,names)
for key in keys(Tcs)
    println("Tc of $(key) = $(Tcs[key])")
end
=#


#=
function compute_Tc(p01,p02,temps,data::Array{Float64,2})

    num_temps = length(temps)
    @assert num_temps == size(data)[2]
    tdata = [i for i in 1:length(data[:,1])]
    τs = zeros(Float64,num_temps)
    for it in 1:num_temps
        ydata = log.(data[:,it])
        τs[it] = compute_τ(tdata,ydata,p01)
    end
end
=#



function plot_τ_and_Tc(τs,Tc)
    plt.figure()
    plt.plot(temps,log.(τs),label="log(τ)")
    plt.plot(Tc,l)
    plt.title("")
    plt.xlabel("T")
    plt.ylabel("")
    plt.legend(loc="upper right")
end


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
    aves

end
            

names = ["Gt","mq_sqrt3_corr", "mq_q0_corr", "afvc_corr", "fvc_corr","m_120degs_corr"]
h5file = "2d_splitx_out.h5"
num_split = 2
num_mcsteps = 201
num_temps = 1
