using Test

include("measure_mc.jl")
include("loop_update_tests.jl")
include("compute_tau.jl")

println("unit test results")

function x_model(num_spins)

    return [fill((1.,0.,0.),num_spins)]
end

function random_model(num_spins)    
    spins = fill((0.,0.,0.),num_spins)
    for i in 1:num_spins
        spins[i] = Tuple(normalize(rand(3)))
    end

    return [spins]
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


function mk_kagome2(L)
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
    
    return kagome
end


function test_compute_τ()

    test_model(t,p) = p[1] .- (1/p[2])*t
    test_pdata = [5.0,5.0]
    test_tdata = LinRange(0,5,10)
    test_ydata = test_model(test_tdata,test_pdata)

    p0 = [0.1,0.1]
    τ = compute_τ(test_tdata,test_ydata,p0)
    
    @test τ ≈ test_pdata[2]

end
test_compute_τ()


function test_compute_Tc()

    test_model(t,p) = p[1]*exp.(p[2]./sqrt.((t .- p[3])))
    test_pdata = [0.1,0.1,0.1]
    test_tdata = LinRange(1,5,10)
    test_ydata = test_model(test_tdata,test_pdata)
    println("test_ydata: ",test_ydata)

    p0 = [1.0,1.0,1.0]
    Tc = compute_Tc(test_tdata,test_ydata,p0)
    @test Tc ≈ test_pdata[3]

end
test_compute_Tc()


#= function
    test_model(t,p) = p[1]*exp.(-(1/p[2])t)
    test_pdata = [5.0,5.0]
    num_data = 10
    test_tdata = LinRange(0,5,num_data)
    test_ydata = test_model(test_tdata,test_pdata)

    num_temps = 10
    test_data = Array{Float64,2}(undef,num_data,num_temps)
    temps = LinRange(0.11,5,num_temps)
    test_Tc = 0.01
    for it in 1:num_temps
        effective_T = sqrt(temps[it] - test_Tc)
        test_data[:,it] = test_ydata .* exp.(-1/effective_T)
    end

    p01 = [0.1,0.1]
    p02 = [0.1,0.1,0.1]
    Tc = compute_Tc(p01,p02,temps,test_data)

    @test Tc ≈ test_Tc

    #test_model2(t,p) = p[1]*exp.(p[2]./(t .- p[3]))

end
=#


function test_JModel()
    # Make Jij.txt
    # Three spins on a chain with nearest neighbobor (J1) and
    # next nearest neighbobor (J2) interactions with an open boundary condition
    J1, J2 = 1.0, -0.1
    num_spins = 3
    J1_idx, J2_idx = 1, 2
    num_unique_Jij = 2
    fname = "Jij_test_JModel.ini"
    open(fname, "w") do f
        println(f, num_spins)
        println(f, num_unique_Jij)
        println(f, "1  $(J1) $(J1) $(J1) 1")
        println(f, "2  $(J2) $(J2) $(J2) 0")
        println(f, 3)
        println(f, "1 2 $(J1_idx)")
        println(f, "2 3 $(J1_idx)")
        println(f, "1 3 $(J2_idx)")
    end

    model = JModel(fname, num_spins)

    @assert model.num_spins == num_spins
    @assert length(model.unique_Jij) == num_unique_Jij
    @assert length(model.Jij) == 3
    @assert model.unique_Jij[1].Jxyz == SVector(J1, J1, J1)
    @assert model.unique_Jij[1].flag_nn == 1
    @assert model.unique_Jij[2].Jxyz == SVector(J2, J2, J2)
    @assert model.unique_Jij[2].flag_nn == 0

    updater = SingleSpinFlipUpdater(model)
    @assert updater.coord_num == [2, 2, 2]
    @assert updater.nn_coord_num == [1, 2, 1]
    @assert updater.nn_sites[:,1] == [2, typemax(SpinIndex)]
    @assert updater.nn_sites[:,2] == [1, 3]
    @assert updater.nn_sites[:,3] == [2, typemax(SpinIndex)]

    spins = [propose_unifo() for _ in 1:num_spins]
    eff_h = effective_field(spins, updater, 1)
    for s in 1:3
       eff_h_ref = J1*spins[2][s] + J2*spins[3][s]
       @assert eff_h[s] ≈ eff_h_ref
    end
end
test_JModel()


function test_compute_m(L)

    num_spins = 3L^2
    q0    = read_spin_config("q0.txt",num_spins)
    sqrt3 = read_spin_config("sqrt3.txt",num_spins)

    utriangles = read_triangles("utriangles.txt",num_spins)
    num_utriangles = length(utriangles)

    #q = read_q()
    qx = LinRange(0,2π,100)
    qy = LinRange(0,2π,100)
    kagome = mk_kagome2(L)
    #kagome = read_kagome()
    
    mq_q0    = []
    mq_sqrt3 = []
    for i in 1:length(qx), j in 1:length(qy)
        push!(mq_q0   ,compute_mq((qx[i],qy[j]),kagome,q0   ,utriangles))
        push!(mq_sqrt3,compute_mq((qx[i],qy[j]),kagome,sqrt3,utriangles))
    end

    @test maximum(mq_q0) ≈ 1.0
    #@test maximum(mq_sqrt3) ≈ 0.5
    
    num_triangles_sisj = 1
    sisj = compute_sisj(num_triangles_sisj, q0, utriangles)
    for it1 in 1:length(utriangles), isub1 in 1:3
        ispin1 = utriangles[it1][isub1]
        for it2 in 1:num_triangles_sisj, isub2 in 1:3
            ispin2 = utriangles[it2][isub2]
            @assert dot(q0[ispin1], q0[ispin2]) ≈ ifelse(isub1 == isub2, 1.0, -0.5)
            @assert sisj[ispin1, isub2, it2] ≈ ifelse(isub1 == isub2, 1.0, -0.5)
        end
    end

    #m_120degrees_q0    = compute_m_120degrees(q0)
    #m_120degrees_sqrt3 = compute_m_120degrees(sqrt3)
    #@test m_120degrees_q0 ≈ m_120degrees_sqrt3 ≈ 1.0

end
L = 3
test_compute_m(L)


function convert_spins_to_array(spins::Vector{HeisenbergSpin})
    N = length(spins)
    spins_array = Array{Float64,2}(undef, 3, N)
    for i in 1:N, j=1:3
        spins_array[j,i] = spins[i][j]
    end
    spins_array
end

function test_compute_vector_chirality(L)

    num_spins = 3L^2
    q0 = read_spin_config("q0.txt",num_spins)
    sqrt3 = read_spin_config("sqrt3.txt",num_spins)

    utriangles = read_triangles("utriangles.txt",num_spins)
    dtriangles = read_triangles("dtriangles.txt",num_spins)

    q0_new = convert_spins_to_array(q0)
    sqrt3_new = convert_spins_to_array(sqrt3)

    num_utriangles = length(utriangles)
    num_dtriangles = length(dtriangles)
    @assert num_utriangles == num_dtriangles

    ferro_vc = compute_ferro_vector_chirality(q0_new,utriangles,dtriangles)
    af_vc    = compute_af_vector_chirality(q0_new,utriangles,dtriangles)
 
    q0_ferro = (num_utriangles+num_dtriangles) * (3*sqrt(3)/2) / num_spins
    q0_ferro = q0_ferro^2 / 3
    q0_af    = 0
    @test isapprox(ferro_vc,q0_ferro) 
    @test isapprox(af_vc,q0_af)
 
    ferro_vc = compute_ferro_vector_chirality(sqrt3_new,utriangles,dtriangles)
    af_vc    = compute_af_vector_chirality(sqrt3_new,utriangles,dtriangles)
    sqrt3_ferro = 0
    sqrt3_af   = (num_utriangles+num_dtriangles) * (3*sqrt(3)/2) / num_spins
    sqrt3_af = sqrt3_af^2 / 3
    #@test isapprox(ferro_vc,sqrt3_ferro) 
    @test isapprox(af_vc,sqrt3_af)
    
    ferro_vc2, af_vc2, vc_corr = compute_vector_chiralities(sqrt3_new,utriangles,dtriangles)

    @test isapprox(ferro_vc, ferro_vc2, rtol=0.0, atol=1e-8)
    @test isapprox(af_vc, af_vc2, rtol=0.0, atol=1e-8)
end
L = 3
test_compute_vector_chirality(L)


# octopolar_v2 and order_parameter are former implementations measuring order parameters.
function octopolar_v2(spins,num_spins::Int64,num_temps::Int64)

    T = zeros(num_temps)

    for i in 1:num_temps
        temp = zeros(3^3)

        for j in 1:num_spins

            spin = spins[i][j]
            index = 1

            for (a,b,c) in Iterators.product(1:3,1:3,1:3)
                temp[index] +=  spin[a]*spin[b]*spin[c] - (spin[a]*delta(b,c) + spin[b]*delta(c,a) + spin[c]*delta(a,b))/5
                index += 1
            end
        end
        T[i] = sum(temp.^2)

    end
    return T/(num_spins^2)
end


# It looks like wrong,need to be fixed.
function order_parameter(spins::Vector{Vector{Tuple{Float64,Float64,Float64}}},
                         num_spins::Int64, num_temps::Int64, q::Tuple{Float64,Float64}, unit_cell::Vector{Tuple{Float64,Float64}})
    
    #unit_cell = make_kagome(num_spins) 

    M2_AF = zeros(num_temps) 

    for temp in 1:num_temps
        gp1 = spins[temp][1:3:num_spins] 
        gp2 = spins[temp][2:3:num_spins] 
        gp3 = spins[temp][3:3:num_spins] 
        
        temp_vec1 = (0.,0.,0.)
        temp_vec2 = (0.,0.,0.)
        temp_vec3 = (0.,0.,0.)

        for i in 1:Int(num_spins/3)

            phase_iqr = -im*dot(q,unit_cell[i])
   
            temp_vec1 = temp_vec1 .+ gp1[i] .* exp(phase_iqr) 
            temp_vec2 = temp_vec2 .+ gp2[i] .* exp(phase_iqr) 
            temp_vec3 = temp_vec3 .+ gp2[i] .* exp(phase_iqr) 
        end

        M2_AF[temp] = norm(temp_vec1)^2 + norm(temp_vec2)^2 + norm(temp_vec3)^2
    end
    
    return M2_AF
end


function test_compute_m2_af()
    
    num_spins = 6
    spins = [fill((0.,0.,0.),num_spins)]
    for i in 1:num_spins
        theta = 10*rand()
        spins[1][i] = (cos(theta),sin(theta),0.)
    end

   
    triangles = [(1,2,3),(4,5,6)]
    m2_af = compute_m2_af(spins[1],triangles)
    

end


function test_compute_T2_op(spins::Vector{Vector{HeisenbergSpin}})
    
    num_spins = length(spins[1])
    T2_op = compute_T2_op(spins[1],num_spins)
    @test isapprox(T2_op,octopolar_v2(spins,num_spins,1)[1])

end


function test_compute_1dcorr(spins::Vector{IsingSpin})
   
    num_spins = length(spins)

    corr_rough = zeros(Float64,num_spins) 
    for i in 1:num_spins
        for j in 1:num_spins
            if num_spins < i+j < 2*num_spins
                corr_rough[j] += spins[i]*spins[mod(i+j,num_spins)]
                continue
            elseif i+j == 2*num_spins
                corr_rough[num_spins] += spins[num_spins]^2
                continue
            end
        corr_rough[j] += spins[i]*spins[i+j]
        end
    end
    prepend!(corr_rough,corr_rough[end])
    pop!(corr_rough)
    corr_rough /= num_spins

    corr_wfft  = compute_1dcorr(spins) 
  
    #println("rough: ",corr_rough)
    #println("wfft: ",corr_wfft)
  
    @test all(isapprox.(corr_rough,corr_wfft))
 
    rdata = [(i-1) for i in 1:num_spins]
    param_init = [0.5,0.5]
    corr_length_wfft  = compute_corr_length(corr_wfft ,rdata,param_init)
    corr_length_rough = compute_corr_length(corr_rough,rdata,param_init)

    println(corr_length_wfft)
    println(corr_length_rough)
  
    @test isapprox(corr_length_wfft,corr_length_rough)

end


    

test_compute_m2_af()
test_compute_T2_op(x_model(6))

num_spins = 100
test_ising_spins = [rand(Int8.([-1,1])) for i in 1:num_spins]
#test_compute_1dcorr(test_ising_spins)
