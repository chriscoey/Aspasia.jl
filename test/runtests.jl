
using Test
import MathOptInterface
const MOI = MathOptInterface
import Aspasia

import GLPK
oa_solver = MOI.OptimizerWithAttributes(
    GLPK.Optimizer,
    "msg_lev" => 0,
    "tol_int" => 1e-9,
    "tol_bnd" => 1e-9,
    "mip_gap" => 0.0,
)

function _run_moi_tests(use_iterative_method::Bool, oa_solver)
    config = MOI.Test.Config(
        atol = 1e-3,
        rtol = 1e-3,
        exclude = Any[
            MOI.ConstraintDual,
            MOI.ConstraintBasisStatus,
            MOI.DualObjectiveValue,
        ],
    )

    aspasia = Aspasia.Optimizer()
    MOI.set(aspasia, MOI.Silent(), true)
    MOI.set(aspasia, MOI.RawOptimizerAttribute("oa_solver"), oa_solver)
    MOI.set(
        aspasia,
        MOI.RawOptimizerAttribute("use_iterative_method"),
        use_iterative_method,
    )
    opt = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.Bridges.full_bridge_optimizer(aspasia, Float64),
    )

    excludes = String[
        # not implemented:
        "test_attribute_SolverVersion",
        # invalid model:
        "test_constraint_ZeroOne_bounds_3",
        # CachingOptimizer does not throw if optimizer not attached:
        "test_model_copy_to_UnsupportedAttribute",
        "test_model_copy_to_UnsupportedConstraint",
        # other:
        "test_objective_qp",
        "test_quadratic",
    ]
    if use_iterative_method
        push!(excludes, "test_unbounded")
    end
    # excludes = String[]

    # includes = String["test_conic_SecondOrderCone", "test_conic_PositiveSemidefiniteCone"]
    includes = String[]

    return MOI.Test.runtests(opt, config, include = includes, exclude = excludes)
end

println("starting MOI tests")
@testset "MOI tests" begin
    println("starting iterative method tests")
    @testset "iterative method" _run_moi_tests(true, oa_solver)
    println("starting OA solver driven method tests")
    @testset "OA solver driven method" _run_moi_tests(false, oa_solver)
    println("finished MOI tests")
end
