#=
Copyright 2019, Chris Coey and contributors

functions and caches for polyhedral approximation algorithms
=#

module Solvers

using Printf

import Aspasia.Cones

abstract type Solver end

# separation cuts algorithm
include("separation_cuts/solver.jl")

get_status(solver::Solver) = solver.status
get_solve_time(solver::Solver) = solver.solve_time
get_num_iters(solver::Solver) = solver.num_iters

end
