include("classical_mc.jl")
using .ClassicalMC

using MPI
using ArgParse

s = ArgParseSettings()
@add_arg_table! s begin
    "input"
        help = "input file"
        required = true
    "--nsplit"
        help = "Number of MPI subgroups"
        arg_type = Int
        default = 1
end
args = parse_args(ARGS, s)

# Initialize MPI environment
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nproc = MPI.Comm_size(comm)
nsplit = args["nsplit"]
input = args["input"]
basename = splitext(input)[1]
if splitext(input)[2] != ".ini"
    error("Extension of input file must be ini!")
end

if nsplit == 1
    mycomm = comm
    prefix = basename * "_"
else
    if mod(nproc, nsplit) != 0
        error("nproc cannot be divided by nsplit!")
    end
    group_size = nproc÷nsplit
    ig = rank÷group_size + 1
    prefix = basename * "_split$(ig)_"
    mycomm = MPI.Comm_split(MPI.COMM_WORLD, ig, rank)
end

big_prime_num = 3571
if rank == 0
    @time solve(args["input"], mycomm, prefix, rank*big_prime_num)
else
    solve(args["input"], mycomm, prefix, rank*big_prime_num)
end

MPI.Finalize()
