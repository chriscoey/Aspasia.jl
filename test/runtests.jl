
using Test
using Printf
import Random
import Aspasia
import MathOptInterface
const MOI = MathOptInterface

# OA model solver:

using GLPK
oam_solver = GLPK.Optimizer
oam_options = (
    # method = :Simplex,
    msg_lev = GLPK.MSG_ON,
    # msg_lev = GLPK.MSG_ERR,
    tol_int = 1e-9,
    tol_bnd = 1e-9,
    tol_dj = 1e-9,
    mip_gap = 1e-9,
    )

# using Clp
# oam_solver = Clp.Optimizer
# oam_options = (PrimalTolerance = 1e-8, DualTolerance = 1e-8, LogLevel = 0)

# using Gurobi
# oam_solver = Gurobi.Optimizer
# oam_options = (
#     OutputFlag = 1,
#     FeasibilityTol = 1e-9,
#     OptimalityTol = 1e-9,
#     IntFeasTol = 1e-9,
#     MIPGap = 1e-8,
# )

@testset "Aspasia tests" begin

@info("starting MathOptInterface tests")
verbose = false
@testset "MOI tests" begin
    include(joinpath(@__DIR__, "moi.jl"))
    test_moi(verbose, oam_solver, oam_options)
end

end
