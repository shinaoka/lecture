include("mcmc.jl")

using LinearAlgebra

function estimate_plane(spins::AbstractArray{HeisenbergSpin})
    #=
    Estimate a coplanar plane and return its (normalized) normal vector
    =#

    num_spins = length(spins)

    # Normal vector perpendicular to the plane
    normal_vec = zeros(Float64, 3)

    for i=1:num_spins
        for j=i+1:num_spins
           cv = cross(collect(spins[i]), collect(spins[j]))
           if dot(cv, normal_vec) < 0
               cv *=  -1
           end
           normal_vec += cv
        end
    end
    normal_vec /= norm(normal_vec)

    return normal_vec
end
