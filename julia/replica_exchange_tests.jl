module ClassicalMC

using Test
using MPI
using Traceur

include("mcmc.jl")
include("replica_exchange.jl")

function test_exchange(comm)
    # Number of temps per process
    num_temps_local = 4
    min_temp, max_temp = 1.0, 100.0
    
    rank = MPI.Comm_rank(comm)
    num_proc = MPI.Comm_size(comm)
    num_temps = num_temps_local * num_proc
    num_spins = 100
    
    temps_init = collect(range(max_temp, stop=min_temp, length=num_proc*num_temps_local))
    rex = ReplicaExchange{HeisenbergSpin}(temps_init, rank*num_temps_local+1,
        rank*num_temps_local+num_temps_local, num_spins)

    config_local = [fill((0., 0., 1.), num_spins) for _ in 1:num_temps_local]
    energy_local = fill(0.0, num_temps_local)
    perform!(rex, config_local, energy_local, comm)

    if rank == 0
        @time perform!(rex, config_local, energy_local, comm)
    else
        perform!(rex, config_local, energy_local, comm)
    end
    #@trace perform!(rex, config_local, energy_local, comm)

    @assert all(rex.num_accepted[1:end-1] .== 2)
end
    
function test_update_dist(comm)
    # Number of temps per process
    num_temps_local = 4
    min_temp, max_temp = 1.0, 100.0
    
    rank = MPI.Comm_rank(comm)
    num_proc = MPI.Comm_size(comm)
    num_temps = num_temps_local * num_proc
    num_spins = 100
    
    #if num_proc == 1
        #error("Please run this test with multiple processes!")
    #end
    
    temps_init = collect(range(max_temp, stop=min_temp, length=num_proc*num_temps_local))
    rex = ReplicaExchange{HeisenbergSpin}(temps_init, rank*num_temps_local+1, rank*num_temps_local+num_temps_local, num_spins)
    
    # Our model: acceptance rate is propotional to the inverse of distance in beta.
    big_int = 1000000
    for iupdate = 1:200
        #if rank==0
            #betas_opt = 1 ./ rex.temps
            #d_betas_opt = abs.(betas_opt[2:end] - betas_opt[1:end-1])
            #println(betas_opt)
            #println(d_betas_opt)
        #end
        rex.num_attemps = big_int
        for i in rex.start_idx:min(rex.end_idx,num_temps-1)
            rex.num_accepted[i-rex.start_idx+1] = round(UInt64, big_int/abs(1/rex.temps[i+1]-1/rex.temps[i]))
        end
        update_temps_dist!(rex, comm)
    end
    
    betas_opt = 1 ./ rex.temps
    d_betas_opt = abs.(betas_opt[2:end] - betas_opt[1:end-1])
    #@assert all(isapprox.(d_betas_opt, 0.33, rtol=1e-3))
    @assert all(isapprox.(d_betas_opt, abs(1/max_temp-1/min_temp)/(num_temps-1), rtol=1e-3))
end

MPI.Init()
comm = MPI.COMM_WORLD
    
test_update_dist(comm)
test_exchange(comm)
    
MPI.Finalize()
end
