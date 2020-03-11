using Test
using Random
include("loop_update.jl")

println("unit test results")

function test_estimate_plane()
    num_spins = 4
    spins = fill((0.,0.,0.), num_spins) 
    for i=1:num_spins 
        theta = 10 * rand() 
        spins[i] = (cos(theta), sin(theta), 0.0) 
    end
    
    @test isapprox(estimate_plane(spins), [0., 0., 1.]) || isapprox(estimate_plane(spins), [0., 0., -1.])
end

function test_estimate_axes()
    num_spins = 4
    spins = fill((0.,0.,0.), num_spins)
    for i=1:num_spins
        theta = 10 * rand()
        spins[i] = (cos(theta), sin(theta), 0.0)
    end
    
    normal_vec = estimate_plane(spins)
    spin       = spins[1]
    x_axis     = estimate_axes(spin,normal_vec)[1]
    @test isapprox(x_axis,collect(spins[1]))
end

function test_paint_black()
  
    num_spins = 100
    num_ref = 20
    indices = [i for i=1:num_ref]
    colors  = [red for i=1:num_spins]
    paint_black(colors,indices)

    @test colors[1:num_ref] == [black for i=1:20]
end

    
function test_paint_rbg_differently()

    num_spins = 6
    spins = fill((0.,0.,0.),num_spins)
    
    for i=1:num_spins
        theta = i*(2*pi)/3         
        spins[i] = (cos(theta),sin(theta),0.)
    end
    
    reference_system = spins[1:3]
    reference_index  = [1,2,3]
    colors = [red for i=1:num_spins]
    colors = paint_rbg_differently(spins,reference_system,reference_index,colors) 
    
    @test colors == [black,black,black,red,blue,green] 
end

function test_find_breaking_triangle()
    num_spins = 3
    Jij = []
    push!(Jij, (1,2,1.,1.,1.,1))
    push!(Jij, (1,3,1.,1.,1.,1))
    push!(Jij, (2,3,1.,1.,1.,1))
    
    model = JModel(num_spins,Jij)
    updater = SingleSpinFlipUpdater(model)
    colors = [red,red,blue]

    @test find_breaking_triangle(updater,colors) == [black,black,black]
end

function test_parallel_flip()
    
    reference_system = [(cos(i*2pi/3),sin(i*2pi/3),0.) for i=1:3]
    spin = reference_system[2] 
    color = red
   
    new_spin = parallel_flip(spin,color,reference_system)
    @test isapprox(collect(new_spin),collect(reference_system[3]))

end

function test_mk_new_spins_on_loop()

    spins_ref = [(cos(i*2pi/3),sin(i*2pi/3),0.) for i=1:3]
    colors_on_loop = (blue,green)
    
    new_spins_on_loop = mk_new_spins_on_loop(spins_ref,colors_on_loop,spins_ref)
    expected_spins = [spins_ref[1],spins_ref[3],spins_ref[2]]
    @test all(isapprox.(collect.(new_spins_on_loop),collect.(expected_spins)))
end

function ring_plus_one_model()
    # 1D system of 4 spins with a periodic boundary condition and alternating red and blue colors
    # One black site is connected to spin 1
    num_spins = 5
    loop_length = 4

    Jij = []
    for ispin=1:loop_length-1
        push!(Jij, (SpinIndex(ispin), SpinIndex(ispin+1), 0.1, 0.5, 1.,1))
    end
    push!(Jij, (SpinIndex(1), SpinIndex(loop_length), 0.1, 0.5, 1.,1))
    push!(Jij, (SpinIndex(1), SpinIndex(num_spins), 0.2, 0.5, 100.,1))
    model = JModel(num_spins, Jij)
    colors = [red, blue, red, blue, black]
    return model, colors
end

function all_to_all_model(num_spins)
    Jij = []
    for ispin=1:num_spins
        for jspin=ispin+1:num_spins
            push!(Jij, (SpinIndex(ispin), SpinIndex(jspin), rand(), rand(), rand(),1))
        end
    end
    model = JModel(num_spins, Jij)
    colors = rand([red, blue, green, black], num_spins)
    return model, colors
end

function test_find_loop(model::JModel, colors::Array{Color})
    u = SingleSpinFlipUpdater(model)
    num_spins = model.num_spins

    work = fill(0, num_spins)
    colors_on_loop = (red, blue)
    start_spin_idx = findfirst(x->x==colors_on_loop[1], colors)
    spin_idx_on_loop = find_loop(u, colors, colors_on_loop, start_spin_idx, num_spins, work, true)

    loop_length = length(spin_idx_on_loop)

    println("loop length: ", loop_length)
    @assert loop_length > 2
    @assert mod(loop_length, 2) == 0

    # Check two colors are alternating on the loop.
    for i in 1:loop_length
        println(i, " ", colors[spin_idx_on_loop[i]])
        if mod(i,2) == 1
            @assert colors[spin_idx_on_loop[i]] == colors[spin_idx_on_loop[1]]
        else
            @assert colors[spin_idx_on_loop[i]] == colors[spin_idx_on_loop[2]]
        end
    end
   
    spins = fill((0.,0.,1.), num_spins)
    new_spins_on_loop = fill((0.,0.,-1.), loop_length)
    dE = compute_dE_loop(u, spin_idx_on_loop, spins, new_spins_on_loop, work, true)

    new_spins = copy(spins)
    for i in 1:loop_length
        new_spins[spin_idx_on_loop[i]] = new_spins_on_loop[i]
    end
    dE_ref = compute_energy(model, new_spins) - compute_energy(model, spins)

    println("dE: ", dE)
    @test isapprox(dE, dE_ref)
end

function test_metropolis_method()
    
    beta = 1.
    dE   = 1.
    
    spins = [(0.,0.,0.) for i=1:10]
    spin_idx_on_loop  = [rand(1:10) for i=1:5]
    new_spins_on_loop = [(1.,1.,1.) for i=1:5]

    dE = metropolis_method(beta,dE,spins,spin_idx_on_loop,new_spins_on_loop) 
    @test 
    
test_estimate_plane()
test_estimate_axes()
test_paint_black()
test_paint_rbg_differently()
test_find_breaking_triangle()
test_parallel_flip()
test_mk_new_spins_on_loop()

Random.seed!(10)
model, colors = ring_plus_one_model()
println("ring_plus_one_model results")
test_find_loop(model, colors)

Random.seed!(10)
model, colors = all_to_all_model(20)
println("all_to_all_model results")
test_find_loop(model, colors)

test_metropolis_method()
