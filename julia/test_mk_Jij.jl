using Test
using Benchmark
include("mk_Jij.jl")
include("mk_Jij_v2.jl")


function test_mk_Jij_kagome_nn(L)

    J1 = (1.,1.,1.)
    J2 = (1.,1.,1.)
    len1 = 2*cos(pi/3)
    len2 = 2*sin(pi/3)

    # test nearest-neighbor interaction v1 and v2.
    # Jij_v1 consists (site1,site2,Jx,Jy,Jz,nn_flag)
    Jij_v1_nn = mk_interaction(L,J1,J2,len1,0.0)
    Jij_v2_nn = mk_Jij_kagome_nn(L,1.,1.,1.,1)
 
    @benchmark mk_interaction(L,J1,J2,len1,0.)
    @benchmark mk_Jij_kagome_nn(L,1.,1.,1.,1)
end

test_mk_Jij_kagome_nn(100)

function test_mk_Jij_kagome_nnn(L)

    # test next-nearest-neighbor interaction v1 and v2.
    # Jij_v1 consists (site1,site2,Jx,Jy,Jz,nn_flag)
    Jij_v1_nnn = mk_interaction(L,J1,J2,len1,len2)
    Jij_v2_nnn = mk_Jij_kagome_nnn(L,J1,J2)

    @benchmark mk_interaction(L,J1,J2,len1,len2)
    @benchmark mk_(L,J1,J2,len1,len2)
end
