#=
Copyright 2019, Chris Coey and contributors

algorithm for mixed-integer conic outer approximation with K* cuts
=#

mutable struct SepCutsSolver <: Solver
    oa_model::JuMP.Model
    cones::Vector{Cones.Cone}

    verbose::Bool
    tol_rel_opt::Float64
    tol_abs_opt::Float64
    tol_feas::Float64
    max_iters::Int
    time_limit::Float64

    status::Symbol
    num_iters::Int
    solve_time::Float64

    function SepCutsSolver(
        oa_model::JuMP.Model,
        cones::Vector{Cones.Cone},
        ;
        verbose::Bool = true,
        tol_rel_opt = 1e-6,
        tol_abs_opt = 1e-7,
        tol_feas = 1e-7,
        max_iters::Int = 500,
        time_limit::Float64 = 3e2,
        )
        solver = new()

        solver.oa_model = oa_model

        solver.verbose = verbose
        solver.tol_rel_opt = tol_rel_opt
        solver.tol_abs_opt = tol_abs_opt
        solver.tol_feas = tol_feas
        solver.max_iters = max_iters
        solver.time_limit = time_limit

        solver.status = :SolveNotCalled
        solver.num_iters = 0
        solver.solve_time = NaN

        return solver
    end
end

function solve(solver::SepCutsSolver)
    solver.status = :SolveCalled
    start_time = time()

    # TODO for now assuming unbounded, but if get rays then need to cut them off

    # @printf("\n%5s %12s %12s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s\n",
    #     "iter", "p_obj", "d_obj", "abs_gap", "rel_gap",
    #     "x_feas", "y_feas", "z_feas", "tau", "kap", "mu",
    #     "gamma", "alpha",
    #     )
    # flush(stdout)

    while true
        @show solver.num_iters

        JuMP.optimize!(solver.oa_model)
        oa_model_status = JuMP.termination_status(solver.oa_model)
        @show oa_model_status
        if oa_model_status == JuMP.OPTIMAL
            # TODO update obj bounds?
        elseif oa_model_status == JuMP.INFEASIBLE
            solver.verbose && println("infeasibility detected; terminating")
        elseif oa_model_status == JuMP.INFEASIBLE_OR_UNBOUNDED
            solver.verbose && println("infeasibility or unboundedness detected; terminating")
        else
            error("OA solver status not handled")
        end

        if solver.verbose
            # @printf("%5d %12.4e %12.4e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n",
            #     solver.num_iters, solver.primal_obj, solver.dual_obj, solver.gap, solver.rel_gap,
            #     solver.x_feas, solver.y_feas, solver.z_feas, solver.tau, solver.kap, solver.mu,
            #     stepper.prev_gamma, stepper.prev_alpha,
            #     )
            # flush(stdout)
        end

        # TODO check convergence

        if solver.num_iters == solver.max_iters
            solver.verbose && println("iteration limit reached; terminating")
            solver.status = :IterationLimit
            break
        end

        if time() - start_time >= solver.time_limit
            solver.verbose && println("time limit reached; terminating")
            solver.status = :TimeLimit
            break
        end

        is_cut_off = false
        for cone in solver.cones
            cuts = check_feas_get_cuts(cone)
            @show length(cuts)
            if !isempty(cuts)
                is_cut_off = true
                for cut in cuts
                    @show JuMP.value(cut)
                    JuMP.@constraint(solver.model, cut >= 0.0)
                end
            end
        end
        if !is_cut_off
            println("point was not cut off")
            solver.status = :Optimal
            break
        end

        solver.num_iters += 1
    end

    solver.solve_time = time() - start_time

    solver.verbose && println("\nstatus is $(solver.status) after $(solver.num_iters) iterations and $(trunc(solver.solve_time, digits=3)) seconds\n")

    return
end
