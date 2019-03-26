#=
Copyright 2019, Chris Coey and contributors

algorithm for mixed-integer conic outer approximation with separation K* cuts (no conic subproblem)
=#

mutable struct SepCutsSolver <: Solver
    verbose::Bool
    tol_rel_opt::Float64
    tol_abs_opt::Float64
    tol_feas::Float64
    max_iters::Int
    time_limit::Float64

    approx_model::MOI.ModelLike
    con_sets::Vector{MOI.AbstractVectorSet}
    con_funs::Vector{MOI.AbstractVectorFunction}

    status::Symbol
    num_iters::Int
    solve_time::Float64
    obj_value::Float64
    obj_bound::Float64

    function SepCutsSolver(
        approx_model::MOI.ModelLike,
        con_sets::Vector{MOI.AbstractVectorSet},
        con_funs::Vector{MOI.AbstractVectorFunction};
        verbose::Bool = true,
        tol_rel_opt = 1e-6,
        tol_abs_opt = 1e-7,
        tol_feas = 1e-7,
        max_iters::Int = 500,
        time_limit::Float64 = 3e2,
        )
        solver = new()

        solver.verbose = verbose
        solver.tol_rel_opt = tol_rel_opt
        solver.tol_abs_opt = tol_abs_opt
        solver.tol_feas = tol_feas
        solver.max_iters = max_iters
        solver.time_limit = time_limit

        solver.approx_model = approx_model
        solver.con_sets = con_sets
        solver.con_funs = con_funs

        solver.status = :SolveNotCalled
        solver.num_iters = 0
        solver.solve_time = NaN
        solver.obj_value = NaN
        solver.obj_bound = NaN

        return solver
    end
end

function solve(solver::SepCutsSolver)
    solver.status = :SolveCalled
    start_time = time()

    # add initial fixed OA cuts eg variable bounds
    for k in eachindex(solver.con_sets)
        con_set = solver.con_sets[k]
        con_fun = solver.con_funs[k]

        cuts = Aspasia.get_init_cuts(con_set)
        for cut in cuts
            if con_fun isa MOI.VectorOfVariables
                cut_terms = MOI.ScalarAffineTerm.(cut, con_fun.variables)
                cut_constant = 0.0
            else
                cut_terms = [MOI.ScalarAffineTerm(cut[t.output_index] * t.scalar_term.coefficient,
                    t.scalar_term.variable_index) for t in con_fun.terms]
                cut_constant = dot(cut, con_fun.constants)
            end

            cut_expr = MOI.ScalarAffineFunction(cut_terms, 0.0)
            cut_ref = MOI.add_constraint(solver.approx_model, cut_expr, MOI.GreaterThan(-cut_constant))
            # TODO store cut_ref with the constraint so can access values/duals
        end
    end

    # @printf("\n%5s %12s %12s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s\n",
    #     "iter", "p_obj", "d_obj", "abs_gap", "rel_gap",
    #     "x_feas", "y_feas", "z_feas", "tau", "kap", "mu",
    #     "gamma", "alpha",
    #     )
    # flush(stdout)

    while true
        MOI.optimize!(solver.approx_model)

        approx_model_status = MOI.get(solver.approx_model, MOI.TerminationStatus())
        if approx_model_status == MOI.OPTIMAL
            solver.obj_bound = MOI.get(solver.approx_model, MOI.ObjectiveBound())
        elseif approx_model_status == MOI.DUAL_INFEASIBLE
            # TODO need to cut off the ray
        elseif approx_model_status == MOI.INFEASIBLE
            solver.verbose && println("infeasibility detected; terminating")
            solver.status = :PrimalInfeasible
            break
        elseif approx_model_status == MOI.INFEASIBLE_OR_UNBOUNDED
            @show approx_model_status
            # solver.verbose && println("infeasibility or unboundedness detected; terminating")
            # solver.status = :ApproxSolverStatusNotHandled
            solver.verbose && println("LP solver returned infeasible or unbounded - we are assuming infeasibility, but we could be wrong")
            solver.status = :PrimalInfeasible
            break
        else
            error("OA solver status not handled")
        end

        # if solver.verbose
            # @printf("%5d %12.4e %12.4e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n",
            #     solver.num_iters, solver.primal_obj, solver.dual_obj, solver.gap, solver.rel_gap,
            #     solver.x_feas, solver.y_feas, solver.z_feas, solver.tau, solver.kap, solver.mu,
            #     stepper.prev_gamma, stepper.prev_alpha,
            #     )
            # flush(stdout)
        # end

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
        for k in eachindex(solver.con_sets)
            con_set = solver.con_sets[k]
            con_fun = solver.con_funs[k]
            # TODO do not include constant if primal solution is a ray?
            fun_val = MOIU.evalvariables(vi -> MOI.get(solver.approx_model, MOI.VariablePrimal(), vi), con_fun)
            @assert !any(isnan, fun_val)

            cuts = Aspasia.get_sep_cuts(fun_val, con_set)
            for cut in cuts
                if con_fun isa MOI.VectorOfVariables
                    cut_terms = MOI.ScalarAffineTerm.(cut, con_fun.variables)
                    cut_constant = 0.0
                else
                    cut_terms = [MOI.ScalarAffineTerm(cut[t.output_index] * t.scalar_term.coefficient,
                        t.scalar_term.variable_index) for t in con_fun.terms]
                    cut_constant = dot(cut, con_fun.constants)
                end

                cut_expr = MOI.ScalarAffineFunction(cut_terms, 0.0)
                cut_ref = MOI.add_constraint(solver.approx_model, cut_expr, MOI.GreaterThan(-cut_constant))
                # TODO store cut_ref with the constraint so can access values/duals
                is_cut_off = true
            end
        end
        if !is_cut_off
            solver.verbose && println("no cuts were added; terminating")
            solver.status = :Optimal
            solver.obj_value = MOI.get(solver.approx_model, MOI.ObjectiveValue())
            break
        end

        solver.num_iters += 1
    end

    solver.solve_time = time() - start_time

    solver.verbose && println("\nstatus is $(solver.status) after $(solver.num_iters) iterations and $(trunc(solver.solve_time, digits=3)) seconds\n")

    return
end
