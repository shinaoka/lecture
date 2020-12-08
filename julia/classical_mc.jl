module ClassicalMC

export solve,get_param,read_temps,read_Jij,read_spin_config,read_triangles

include("mcmc.jl")
include("accumulator.jl")
include("replica_exchange.jl")
include("loop_update.jl")
include("measure_mc.jl")

using Random
using ConfParser
using ArgParse
using MPI
using Test
using CPUTime
using HDF5
using Profile


# Read a list of temperatures
function read_temps(temperature_file::String)
    temps = Vector{Float64}(undef, 0)
    num_temps = 0
    open(temperature_file) do file
        num_temps = parse(Int64, readline(file))
        for l in 1:num_temps
            temp = parse(Float64, readline(file))
            push!(temps, temp)
        end
    end

    # Check if temperatures are in ascending order or descending order
    if !all(temps[1:num_temps-1] .< temps[2:num_temps]) && !all(temps[1:num_temps-1] .> temps[2:num_temps])
        error("Temperatures must be given either in ascending order or in descending order!")
    end

    return temps
end


# Distribute temperature over MPI processes
function distribute_temps(rank, num_temps, num_proc)
    num_temps_local = fill(trunc(Int, num_temps/num_proc), (num_proc,))
    left_over = mod(num_temps, num_proc)
    num_temps_local[1:left_over] .+= 1
    @test sum(num_temps_local) == num_temps
    start_idx = sum(num_temps_local[1:rank]) + 1
    end_idx = start_idx + num_temps_local[rank+1] - 1

    return start_idx, end_idx
end

# Read non-zero elements in the right-upper triangle part of Jij
"""
function read_Jij(Jij_file::String, num_spins::Int64)
    Jij = Vector{Tuple{SpinIndex,SpinIndex,Float64,Float64,Float64,Int64}}(undef, 0)
    open(Jij_file, "r" ) do fp
        
        @assert num_spins == parse(Int64, readline(fp)) "!match num_spins. See 2d.ini and head of Jij.txt"

        num_Jij_elems = parse(Int64, readline(fp))
        for i in 1:num_Jij_elems
            str     = split(readline(fp))
            i       = parse(SpinIndex, str[1])
            j       = parse(SpinIndex, str[2])
            val_x   = parse(Float64,   str[3])
            val_y   = parse(Float64,   str[4])
            val_z   = parse(Float64,   str[5])
            flag_nn = parse(Int64,     str[6])
            if i >= j
                error("Only right-upper triangle part must be given.")
            end
            if i > num_spins || i < 0  || j > num_spins || j < 0
                error("i or j is out of the range [1, num_spins].")
            end
            push!(Jij, (i, j, val_x, val_y, val_z,flag_nn))
        end
    end
    return Jij
end
"""


function read_triangles(file_name::String,num_spins::Int64)
   
    L = Int(sqrt(num_spins/3))
    triangles = fill((0,0,0),L^2)

    open(file_name,"r") do fp

        @assert L^2 == parse(Int64, readline(fp)) "!match number of triangles. See 2d.ini and head of utriangles.txt and dtriangles.txt!"

        for i in 1:L^2
            str = split(readline(fp))
            s1 = parse(Int64, str[1])
            s2 = parse(Int64, str[2])
            s3 = parse(Int64, str[3])
            triangles[i] = (s1,s2,s3)
        end
    end

    return triangles
end


function mk_kagome(L)
    kagome = Dict()
    idx = 1
    L = 2L

    for (i,j) in Iterators.product(1:L,1:L)
        mod_site = mod.((i,j),L)
        if mod(mod_site[1],2) == 0 && mod(mod_site[2],2) == 1
            continue
        end

        kagome[idx] = copy.(mod_site)
        idx += 1
    end
    [kagome[i] for i in 1:length(keys(kagome))]
end



function read_spin_config(file_name::String,num_spins::Int64)

    spins = fill((0.,0.,0.),num_spins)

    open(file_name,"r") do fp

        @assert num_spins == parse(Int64, readline(fp)) "!match num_spins. See 2d.ini and head of spin_config.txt "

        for i in 1:num_spins
            str = split(readline(fp))
            sx = parse(Float64, str[1])
            sy = parse(Float64, str[2])
            sz = parse(Float64, str[3])
            
            @assert isapprox(sx^2+sy^2+sz^2, 1.0)

            spins[i] = (sx,sy,sz)
        end
    end

    return spins
end

function write_spin_config(file_name::String,spins)
    
    num_spins = length(spins)
    open(file_name,"w") do fp
        println(fp,num_spins)
        for i in 1:num_spins 
            sx,sy,sz = spins[i]
            println(fp,sx," ",sy," ",sz)
        end
    end

end
         

function compute_magnetization(acc::Accumulator,num_spins::Int64,spins::Vector{Vector{Tuple{Float64,Float64,Float64}}},num_temps::Int64)
    
    mx = zeros(Float64, num_temps)
    my = zeros(Float64, num_temps)
    mz = zeros(Float64, num_temps)
    m2 = zeros(Float64, num_temps)

    for i in 1:num_temps
      temp_mx = 0.0    
      temp_my = 0.0    
      temp_mz = 0.0    
      for j in 1:num_spins
        temp_mx += spins[i][j][1]
        temp_my += spins[i][j][2]
        temp_mz += spins[i][j][3]
      end
      mx[i] = temp_mx
      my[i] = temp_my
      mz[i] = temp_mz
      m2[i] = float(mx[i]^2 + my[i]^2+ mz[i]^2)
    end
    
    add!(acc, "Mx2",mx.^2)    
    add!(acc, "My2",my.^2)    
    add!(acc, "Mz2",mz.^2)
    #add!(acc, "M2", m2)
    #add!(acc, "M4", m2.^2)
   
 
end


function get_param(type, conf, block, key, default_value)
    if haskey(conf, block, key)
        return parse(type, retrieve(conf, block, key))
    else
        return convert(type, default_value)
    end
end

function save_result(h5file, acc, comm, res_names)
    rank = MPI.Comm_rank(comm)
    for r in res_names
        data = mean_gather(acc, r, comm)
        if rank > 0
            continue
        end
        h5file["obs/"*r] = data
    end
end


function solve(input_file::String, comm, prefix, seed_shift)
    open(prefix*"out", "w") do outf
        solve_(input_file::String, comm, prefix, seed_shift, outf)
    end
end

function solve_(input_file::String, comm, prefix, seed_shift, outf)
    if !isfile(input_file)
        error("$input_file does not exists!")
    end
    conf = ConfParse(input_file)
    parse_conf!(conf)
    
    rank = MPI.Comm_rank(comm)
    num_proc = MPI.Comm_size(comm)

    num_spins = get_param(Int64, conf, "model", "num_spins", 0)
    Jij_file = retrieve(conf, "model", "Jij")
    temperature_file = retrieve(conf, "model", "temperatures")
    is_xy = get_param(Bool, conf, "model", "xy_spins", false)
    if rank == 0 && is_xy
        println(outf, "Using XY spins")
    end

    num_sweeps       = parse(Int64, retrieve(conf, "simulation", "num_sweeps"))
    num_therm_sweeps = parse(Int64, retrieve(conf, "simulation", "num_therm_sweeps"))
    meas_interval    = parse(Int64, retrieve(conf, "simulation", "meas_interval"))
    ex_interval      = parse(Int64, retrieve(conf, "simulation", "ex_interval"))
    seed             = parse(Int64, retrieve(conf, "simulation", "seed"))
    opt_temps_dist   = get_param(Bool,       conf, "simulation", "opt_temps_dist", true)
    min_attemps_update_temps_dist  = get_param(Int64,  conf, "simulation", "min_attemps_update_temps_dist", 100)

    # Si Sj
    num_src_triangles_sisj = get_param(Int64, conf, "simulation", "num_src_triangles_sisj", 5)

    # Non-equilibrium relaxation method.
    use_neq = get_param(Bool, conf, "simulation", "neq", false)

    if opt_temps_dist && 5*min_attemps_update_temps_dist*ex_interval > num_therm_sweeps/2
        error("num_therm_sweeps is too small for optimizing temps_dist!")
    end

    # For loop updates
    loop_num_trial  = parse(Int64, retrieve(conf, "loop_update", "num_trial"))
    max_loop_length  = parse(Int64, retrieve(conf, "loop_update", "max_loop_length"))
    loop_interval  = get_param(Int64, conf, "loop_update", "interval", 10)

    # Read a list of temperatures
    temps = read_temps(temperature_file)
    num_temps = length(temps)
    if num_temps < num_proc
        error("Number of processes > num_temps")
    end
    if rank == 0
       println(outf, "num of temperatures = ", num_temps)
    end

    # Decide which temperatures are computed on this process
    start_idx, end_idx = distribute_temps(rank, num_temps, num_proc)
    num_temps_local = end_idx - start_idx + 1

    # Read non-zero elements in the right-upper triangle part of Jij
    model = JModel(Jij_file, num_spins)

    # Create single-spin flip updater
    updater = SingleSpinFlipUpdater(model)

    # Init random number generator
    Random.seed!(seed + seed_shift)

    # Create accumulator
    acc = Accumulator(num_temps_local)

    # Create accumulator for collecting stat for every process
    acc_proc = Accumulator(1)


    # Init spins
    spins_local = [fill((1.,0.,0.),num_spins) for _ in 1:num_temps_local]

    # Work arrays
    spins_array = [zeros(Float64, 3, num_spins) for _ in 1:num_temps_local]
    
    # Optional init spin configuration.
    is_read_spin_config = get_param(Bool, conf, "simulation", "read_spin_config", false)
    if use_neq && !is_read_spin_config
        error("Set read_spin_config for neq!")
    end
    if is_read_spin_config
        spin_config_file = retrieve(conf, "model", "spin_config")
        spins_local = fill(read_spin_config(spin_config_file,num_spins),num_temps_local)
    end

    # Create upward and downward triangles for computing order parameters.
    utriangles_file    = retrieve(conf, "model", "utriangles")
    upward_triangles   = read_triangles(utriangles_file,num_spins)
    dtriangles_file    = retrieve(conf, "model", "dtriangles")
    downward_triangles = read_triangles(dtriangles_file,num_spins)

    # preoaration for computation of magnetic order parameter mq
    kagome = mk_kagome(Int(sqrt(num_spins/3)))
    site_pos = [kagome[i] ./ 2 for i in 1:num_spins]
    qs   = [(0.0,0.0),(π/3,π/3)]


    energy_local = [compute_energy(model, spins_local[it]) for it in 1:num_temps_local]
    energy_local = [compute_energy(model, spins_local[it]) for it in 1:num_temps_local]
    # Replica exchange
    rex = ReplicaExchange{HeisenbergSpin}(temps, start_idx, end_idx, num_spins)
    temps = 0

    # Perform MC
    last_output_time = time_ns()
    if rank == 0
        println(outf, "Starting simulation...")
    end
  
    # Create LoopUpdater 
    loop_updater = LoopUpdater{HeisenbergSpin}(num_spins, max_loop_length)

    # For measuring acceptance rates
    single_spin_flip_acc = zeros(Float64, num_temps_local)

    # Replica exchange
    ex_rex = get_param(Bool, conf, "simulation", "ex_rex", false)

    # Non-equilibrium relaxation method.
    if use_neq
        if ex_rex
            error("Do not use replica exchange MC with neq!")
        end
        if num_therm_sweeps > 0
            error("Please set num_therm_sweeps to 0 for neq!")
        end
        for it in 1:num_temps_local
            convert_spins_to_array!(spins_local[it], spins_array[it])
        end

        init_spins_local = [copy(s) for s in spins_local]

        # initial value of Ferro and AF vector spin chirality.
        init_fvc        = zeros(Float64,num_temps_local)
        init_afvc       = zeros(Float64,num_temps_local)
        init_mq_q0      = zeros(Float64,num_temps_local)
        init_mq_sqrt3   = zeros(Float64,num_temps_local)
        init_m_120degs  = zeros(Float64,num_temps_local)
        for it in 1:num_temps_local
            init_fvc[it]  = compute_ferro_vector_chirality(spins_array[it],upward_triangles,downward_triangles) 
            init_afvc[it]  = compute_af_vector_chirality(spins_array[it],upward_triangles,downward_triangles) 
            init_mq_q0[it]    = compute_mq((0.,0.),site_pos,spins_local[it],upward_triangles)
            init_mq_sqrt3[it] = compute_mq((1/3,1/3),site_pos,spins_local[it],upward_triangles)
            init_m_120degs[it]= compute_m_120degrees(spins_local[it])
        end
        init_mq_sqrt3 .*= 2 #max value of order parameter for √3×√3 is 0.5

        correlation_func      = [Float64[] for _ in 1:num_temps_local]
        fvc_correlation       = [Float64[] for _ in 1:num_temps_local]
        afvc_correlation      = [Float64[] for _ in 1:num_temps_local]
        mq_q0_correlation     = [Float64[] for _ in 1:num_temps_local]
        mq_sqrt3_correlation  = [Float64[] for _ in 1:num_temps_local]
        m_120degs_correlation = [Float64[] for _ in 1:num_temps_local]
        for it in 1:num_temps_local
            tmp = sum([dot(spins_local[it][i],spins_local[it][i]) for i in 1:num_spins])
            push!(correlation_func[it],tmp / num_spins)
        
            temp_fvc, temp_afvc, vc_corr = compute_vector_chiralities(spins_array[it],upward_triangles,downward_triangles) 
            push!(fvc_correlation[it], init_fvc[it]*temp_fvc)
            push!(afvc_correlation[it], init_afvc[it]*temp_afvc)

            temp_mq_q0     = compute_mq((0.,0.),site_pos,spins_local[it],upward_triangles)
            temp_mq_sqrt3  = compute_mq((1/3,1/3),site_pos,spins_local[it],upward_triangles)
            temp_m_120degs = compute_m_120degrees(spins_local[it])
            push!(mq_q0_correlation[it],init_mq_q0[it]*temp_mq_q0)
            push!(mq_sqrt3_correlation[it],init_mq_sqrt3[it]*2temp_mq_sqrt3)
            push!(m_120degs_correlation[it],init_m_120degs[it]*temp_m_120degs)
        end

        maf_time_evo = [[] for it in 1:num_temps_local]
        for it in 1:num_temps_local
            push!(maf_time_evo[it],compute_m2_af(spins_local[it],upward_triangles))
        end
    end


    for sweep in 1:num_sweeps
        # Output roughtly every 10 sececonds
        if rank == 0 && time_ns() - last_output_time > 1e+10
            println(outf, "Done $sweep sweeps")
            last_output_time = time_ns()
            flush(stdout)
        end

        elpsCPUtime = []
 
        # Single spin flips
        ts_start = CPUtime_us()
        
        for it in 1:num_temps_local
            
            dE, acc_rate = gaussian_move(updater, 1/rex.temps[it+start_idx-1], model, spins_local[it], is_xy)
            
            energy_local[it] += dE
            single_spin_flip_acc[it] = acc_rate
        end
        #println(outf, "one_sweep", " ",ts_end - ts_start)
        push!(elpsCPUtime, CPUtime_us() - ts_start)

        # Check if energy is correct
        if mod(sweep, 100) == 0
            energy_local = [compute_energy(model, spins_local[it]) for it in 1:num_temps_local]
            for it in 1:num_temps_local
                @assert abs(energy_local[it] - compute_energy(model, spins_local[it])) < 1e-5
            end
        end
         
        # Replica exchange
        ts_start = CPUtime_us()
        if mod(sweep, ex_interval) == 0 && ex_rex == true
            perform!(rex, spins_local, energy_local, comm)
        end
        if opt_temps_dist
            if sweep <= Int(num_therm_sweeps/2) && rex.num_attemps >= min_attemps_update_temps_dist
                update_temps_dist!(rex,comm)
            end
        end

        push!(elpsCPUtime, CPUtime_us() - ts_start)

        # Loop update
        ts_start = CPUtime_us()
        loop_found_rate = zeros(Float64, num_temps_local)
        loop_acc_rate = zeros(Float64, num_temps_local)
        if loop_num_trial > 0 && mod(sweep, loop_interval) == 0
            for it in 1:num_temps_local
                dE, loop_found_rate[it], loop_acc_rate[it] = multi_loop_update!(loop_updater, loop_num_trial,
                    updater,
                    1/rex.temps[it+start_idx-1],
                    max_loop_length, spins_local[it], rank==0)
                energy_local[it] += dE
            end
        end
        push!(elpsCPUtime, CPUtime_us() - ts_start)

        # Measurement
        ts_start = CPUtime_us()
        if sweep >= num_therm_sweeps && mod(sweep, meas_interval) == 0
            # Convert spin config to arrays
            for it in 1:num_temps_local
                convert_spins_to_array!(spins_local[it], spins_array[it])
            end

            add!(acc, "E", energy_local)
            add!(acc, "E2", energy_local.^2)
            add!(acc, "single_spin_flip_acc", single_spin_flip_acc)
            
            add!(acc, "loop_found_rate" , loop_found_rate)
            add!(acc, "loop_accept_rate", loop_acc_rate)
          
            # loop length,candidate order parameter
            measured_loop_length = zeros(Int64,num_temps_local)
            """
            for it in 1:num_temps_local 
                measured_loop_length[it] = compute_loop_length(spins_local[it],
                                                               updater,
                                                               loop_updater,
                                                               max_loop_length,
                                                               false)
            end
 
            if !in(0,measured_loop_length)
                add!(acc,"loop_length",measured_loop_length)
            end
            """
            add!(acc,"loop_length",measured_loop_length)

            # square and fourth power of magnetic and chirality order parameters
            m2_af = zeros(Float64,num_temps_local)
            T2_op = zeros(Float64,num_temps_local)
            for it in 1:num_temps_local
                m2_af[it] = compute_m2_af(spins_local[it],upward_triangles)
                T2_op[it] = compute_T2_op(spins_local[it],num_spins)
            end
            add!(acc, "m_af_2", m2_af)
            add!(acc, "T_op_2", T2_op)
            add!(acc, "m_af_4", m2_af.^2)
            add!(acc, "T_op_4", T2_op.^2)

            mq_q0     = zeros(Float64,num_temps_local)
            mq_sqrt3  = zeros(Float64,num_temps_local)
            m_120degs = zeros(Float64,num_temps_local)
            ss        = [zeros(Float64,3,num_spins) for _ in 1:num_temps_local]
            sisj      = Vector{Array{Float64,3}}(undef, num_temps_local)
            for it in 1:num_temps_local
                mq_q0[it]    = compute_mq((0.,0.),site_pos,spins_local[it],upward_triangles)
                mq_sqrt3[it] = compute_mq((1/3,1/3),site_pos,spins_local[it],upward_triangles)
                m_120degs[it]= compute_m_120degrees(spins_local[it])
                sisj[it] = compute_sisj(num_src_triangles_sisj, spins_local[it], upward_triangles)
                #for is in 1:num_spins, j in 1:3
                    #temp_idx = upward_triangles[1][j]
                    #ss[it][j,is] = spins_local[it][is]⋅spins_local[it][temp_idx]
                #end
            end
            add!(acc,"mq0_2",mq_q0)
            add!(acc,"msqrt_2",mq_sqrt3)
            add!(acc,"m120degs_2",m_120degs)
            add!(acc,"mq0_4",mq_q0.^2)
            add!(acc,"msqrt_4",mq_sqrt3.^2)
            add!(acc,"m120degs_4",m_120degs.^2)
            add!(acc,"ss",ss)
            add!(acc,"sisj", sisj)

            # ferro and anti-ferro vector spin chirality
            fvc  = zeros(Float64,num_temps_local)
            afvc = zeros(Float64,num_temps_local)
            vc_corrs = Vector{Float64}[]
            for it in 1:num_temps_local
                fvc[it], afvc[it], vc_corr = compute_vector_chiralities(spins_array[it],upward_triangles,downward_triangles)
                push!(vc_corrs, vc_corr)
            end
            add!(acc,"Ferro_vc_2",fvc)
            add!(acc,"AF_vc_2",afvc)
            add!(acc,"Ferro_vc_4",fvc.^2)
            add!(acc,"AF_vc_4", afvc.^2)
            add!(acc,"vc_corr", vc_corrs)

            # non-equilibrium relaxation method.
            if use_neq
                for it in 1:num_temps_local

                    tmp = sum([dot(init_spins_local[it][i],spins_local[it][i]) for i in 1:num_spins])
                    push!(correlation_func[it],tmp/num_spins)

                    push!(fvc_correlation[it],init_fvc[it]*fvc[it])
                    push!(afvc_correlation[it],init_afvc[it]*afvc[it])

                    push!(mq_q0_correlation[it],init_mq_q0[it]*mq_q0[it])
                    push!(mq_sqrt3_correlation[it],init_mq_sqrt3[it]*2mq_sqrt3[it])
                    push!(m_120degs_correlation[it],init_m_120degs[it]*m_120degs[it])
                    if mod(sweep, 100) == 0
                        #println(outf, sweep, "th: ", tmp/num_spins)
                    end
  
                end

                """ 
                for it in 1:num_temps_local
                    temp_maf = compute_m2_af(spins_local[it],upward_triangles)
                    push!(maf_time_evo[it],maf_time_evo[it][1]*temp_maf)
                end 
                """
            end
        end
        push!(elpsCPUtime, CPUtime_us() - ts_start)

        if sweep > num_therm_sweeps
            add!(acc_proc, "CPUtime", [Vector{Float64}(elpsCPUtime)])
        end
    end        

    #for it in 1:num_temps_local
      #dE,num_accept = multi_loop_update!(loop_num_trial,loop_num_reference_sites,updater,1/rex.temps[it+start_idx-1],triangles,max_loop_length,spins_local[it],rank==0)
    #end
 
    # Output results
    E = mean_gather(acc, "E", comm)
    E2 = mean_gather(acc, "E2", comm)
    single_spin_flip_acc = mean_gather(acc, "single_spin_flip_acc", comm)
    loop_found_rate = mean_gather(acc,"loop_found_rate", comm)
    loop_accept_rate = mean_gather(acc,"loop_accept_rate", comm)
    CPUtime = mean_gather_array(acc_proc, "CPUtime", comm)
    ave_loop_length = mean_gather(acc,"loop_length",comm)
    m2_af = mean_gather(acc, "m_af_2", comm)
    T2_op = mean_gather(acc, "T_op_2", comm)
    m2q_q0 = mean_gather(acc, "mq0_2", comm)
    m2q_sqrt3 = mean_gather(acc, "msqrt_2", comm)
    m120degs = mean_gather(acc, "m120degs_2", comm)
    Ferro_vc2 = mean_gather(acc, "Ferro_vc_2", comm)
    AF_vc2    = mean_gather(acc, "AF_vc_2"   , comm)
    m4_af = mean_gather(acc, "m_af_4", comm)
    T4_op = mean_gather(acc, "T_op_4", comm)
    m4q_q0 = mean_gather(acc, "mq0_4", comm)
    m4q_sqrt3 = mean_gather(acc, "msqrt_4", comm)
    m120degs4 = mean_gather(acc, "m120degs_4", comm)
    Ferro_vc4 = mean_gather(acc, "Ferro_vc_4", comm)
    AF_vc4    = mean_gather(acc, "AF_vc_4"   , comm)
    ss    = mean_gather_array(acc, "ss" , comm)
    vc_corr = mean_gather_array(acc, "vc_corr" , comm)
    flush(stdout)
    MPI.Barrier(comm)

    flush(stdout)
    MPI.Barrier(comm)
  
    for it in 1:num_temps_local
        #write_spin_config("spin_configs/spin_config$(it+start_idx-1).txt",spins_local[it])
    end
  

    # To save to HDF5 file,add correlation functions to accumulator.
    add!(acc,"Gt",correlation_func)
    add!(acc,"fvc_corr",fvc_correlation)
    add!(acc,"afvc_corr",afvc_correlation)
    add!(acc,"mq_q0_corr",mq_q0_correlation)
    add!(acc,"mq_sqrt3_corr",mq_sqrt3_correlation)
    add!(acc,"m_120degs_corr",m_120degs_correlation)

    # Output time evolution of order parameter.
    #=
    if use_neq
        for itemp in 1:num_temps_local

            open(prefix*"Gt_$(itemp+start_idx-1).dat","w") do fp
               for itime in 1:length(correlation_func[itemp])
                   println(fp, itime, " ", correlation_func[itemp][itime])
               end
            end

            open(prefix*"maf_$(itemp+start_idx-1).dat","w") do fp
               #for itime in 1:length(maf_time_evo[itemp])
                   #println(fp, itime, " ", maf_time_evo[itemp][itime])
               #end
            end

            open(prefix*"fvc_$(itemp+start_idx-1).dat","w") do fp
               for itime in 1:length(fvc_correlation[itemp])
                   println(fp, itime, " ", fvc_correlation[itemp][itime])
               end
            end

            open(prefix*"afvc_$(itemp+start_idx-1).dat","w") do fp
               for itime in 1:length(afvc_correlation[itemp])
                   println(fp, itime, " ", afvc_correlation[itemp][itime])
               end
            end

            open(prefix*"mq_q0_$(itemp+start_idx-1).dat","w") do fp
                for itime in 1:length(mq_q0_correlation[itemp])
                   println(fp, itime, " ", mq_q0_correlation[itemp][itime])
                end
            end

            open(prefix*"mq_sqrt3_$(itemp+start_idx-1).dat","w") do fp
                for itime in 1:length(mq_sqrt3_correlation[itemp])
                   println(fp, itime, " ", mq_sqrt3_correlation[itemp][itime])
                end
            end

            open(prefix*"m_120degs_$(itemp+start_idx-1).dat","w") do fp
                for itime in 1:length(m_120degs_correlation[itemp])
                   println(fp, itime, " ", m_120degs_correlation[itemp][itime])
                end
            end

        end
    end
    =#

    if rank == 0
        println(outf)
        open(prefix*"E.txt", "w") do f
            println(f, "#T <E> <E^2> <C>")
            for i in 1:num_temps
                println(f, "$(rex.temps[i]) $(E[i]) $(E2[i]) $(((E2[i]  - E[i]^2) / (rex.temps[i]^2)) / num_spins)")
                println(outf, "$(rex.temps[i]) $(E[i]) $(E2[i]) $(((E2[i]  - E[i]^2) / (rex.temps[i]^2)) / num_spins)")
            end
        end
      
        println(outf, "single_spin_flip_acc: ", single_spin_flip_acc)
        println(outf, "Acceptant rate of loop update: ")
        for i in 1:num_temps
            println(outf, rex.temps[i], " ", loop_found_rate[i], " ", loop_accept_rate[i])
        end

        println(outf, "<CPUtime> ")
        for (i, t) in enumerate(CPUtime)
            println(outf, " rank=", i-1, " : $t")
        end
    
        # update initial temperature distribution.        
        open(prefix*"temperatures.txt","w") do fp
             println(fp,num_temps)
             for i in 1:num_temps
                 println(fp,rex.temps[i])
             end
        end
     
        for i in 1:num_temps
            println(outf, "af2: $(rex.temps[i]) $(m2_af[i])")
            println(outf, "op2: $(rex.temps[i]) $(T2_op[i])")
            println(outf, "Ferro_vc2: $(rex.temps[i]) $(Ferro_vc2[i])")
            println(outf, "AF_vc2: $(rex.temps[i]) $(AF_vc2[i])")
            println(outf, "m2q0: $(rex.temps[i]) $(m2q_q0[i])")
            println(outf, "m2_sqrt3: $(rex.temps[i]) $(m2q_sqrt3[i])")
            println(outf, "m120degs: $(rex.temps[i]) $(m120degs[i])")
            println(outf, "af4: $(rex.temps[i]) $(m4_af[i])")
            println(outf, "op4: $(rex.temps[i]) $(T4_op[i])")
            println(outf, "Ferro_vc4: $(rex.temps[i]) $(Ferro_vc4[i])")
            println(outf, "AF_vc4: $(rex.temps[i]) $(AF_vc4[i])")
            println(outf, "m4q0: $(rex.temps[i]) $(m4q_q0[i])")
            println(outf, "m4_sqrt3: $(rex.temps[i]) $(m4q_sqrt3[i])")
            println(outf, "m120degs4: $(rex.temps[i]) $(m120degs4[i])")
        end

        h5open(prefix*"corr.h5", "w") do fid
            # ss is an array of arrays of shape (3, num_spins).
            # We "reshape" it into an array of shape (3, num_spins, num_temps).
            #tmp = cat(ss[:,:,CartesianIndex()]..., dims=3)
            #println(size(tmp))
            fid["ss"] = cat(ss[:,:,CartesianIndex()]..., dims=3)
            fid["vc"] = cat(vc_corr[:,CartesianIndex()]..., dims=2)
        end
    end
    flush(stdout)
    MPI.Barrier(comm)

    fid = nothing
    if rank == 0
        fid = h5open(prefix*"out.h5" , "w")
        fid["temperatures"] = temps
    end
    save_to_hdf5!(acc, fid, comm)
    if rank == 0
        close(fid)
    end

    """
    obs_names = String[]
    for base_name in ["m_af", "T_op", "mq0", "msqrt", "m120degs", "Ferro_vc", "AF_vc"]
        push!(obs_names, base_name*"_2")
        push!(obs_names, base_name*"_4")
    end
    fid = nothing
    if rank == 0
        fid = h5open(prefix*"out.h5" , "w")
        fid["temperatures"] = temps
    end
    save_result(fid, acc, comm, obs_names)
    """

    # Stat of Replica Exchange MC
    print_stat(rex, comm, outf)
end

end
