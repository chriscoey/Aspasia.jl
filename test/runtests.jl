
using Test
import MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
import Aspasia

import GLPK
oa_solver = MOI.OptimizerWithAttributes(
    GLPK.Optimizer,
    "msg_lev" => 0,
    "tol_int" => 1e-9,
    "tol_bnd" => 1e-7,
    "mip_gap" => 0.0,
)
config = MOI.Test.Config(
    atol = 1e-4,
    rtol = 1e-4,
    exclude = Any[
        # MOI.ConstraintPrimal,
        MOI.ConstraintDual,
        MOI.ConstraintBasisStatus,
        MOI.DualObjectiveValue,
    ],
)

aspasia = Aspasia.Optimizer()
# MOI.set(aspasia, MOI.Silent(), true)
# MOI.set(aspasia, MOI.RawOptimizerAttribute("tol_feas"), 1e-6)
MOI.set(aspasia, MOI.RawOptimizerAttribute("oa_solver"), oa_solver)

opt = MOI.Utilities.CachingOptimizer(
    MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    MOI.Bridges.full_bridge_optimizer(aspasia, Float64),
)

excludes = String[
    # not implemented:
    "test_attribute_SolverVersion",
    # invalid model:
    "test_constraint_ZeroOne_bounds_3",
    "test_linear_VectorAffineFunction_empty_row",
    # CachingOptimizer does not throw if optimizer not attached:
    "test_model_copy_to_UnsupportedAttribute",
    "test_model_copy_to_UnsupportedConstraint",
]

includes = String["test_conic_SecondOrderCone",]

MOI.Test.runtests(
    opt,
    config,
    include = includes,
    # exclude = excludes,
)

# using Clp
# oa_solver = Clp.Optimizer
# oam_options = (PrimalTolerance = 1e-8, DualTolerance = 1e-8, LogLevel = 0)

# using Gurobi
# oa_solver = Gurobi.Optimizer
# oam_options = (
#     OutputFlag = 1,
#     FeasibilityTol = 1e-9,
#     OptimalityTol = 1e-9,
#     IntFeasTol = 1e-9,
#     MIPGap = 1e-8,
# )
