include("classical_mc.jl")
using .ClassicalMC

using MPI
using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "input"
        help = "input file"
        required = true
end
args = parse_args(ARGS, s)

# Initialize MPI environment
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

if rank == 0
    @time solve(args["input"], comm)
else
    solve(args["input"], comm)
end

MPI.Finalize()
