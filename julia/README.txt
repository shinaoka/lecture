# Run install_package.jl once to install all required packages
julia install_package.jl

# Generate input files
python mk_input.py

# Run Monte Carlo simulations
julia main.jl 1d.ini
