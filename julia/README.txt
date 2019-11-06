# Install OpenMPI
# (for Mac)
brew install openmpi

# (for Windows) TO DO
# ??

# Run install_package.jl once to install all required packages
julia install_package.jl

# Generate input files
julia mk_input.jl

# Run Monte Carlo simulations (with 2 CPU cores in this case)
mpirun -np 2 julia main.jl 1d.ini
