#=
MOI.Test linear and conic tests
=#

using Test
import MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
import Aspasia

function test_moi(verbose, oam_solver, oam_options)
    optimizer = MOIU.CachingOptimizer(MOIU.UniversalFallback(
        MOIU.Model{Float64}()), Aspasia.Optimizer(
        oam_solver, oam_options, verbose = verbose,
        time_limit = 2e1,
        # tol_rel_opt = 1e-5,
        # tol_abs_opt = 1e-5,
        # tol_feas = 1e-5,
        ))

    config = MOIT.Config(
        atol = 1e-3,
        rtol = 1e-3,
        exclude = Any[MOI.VariableBasisStatus, MOI.ConstraintBasisStatus],
        )

    @testset "conic tests" begin
        MOIT.runtests(
            optimizer, config,
            include = ["conic"],
            exclude = ["linear"],
            )
    end

    return
end
