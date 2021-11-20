"""
outer approximation algorithm for mixed-integer conic problems
"""

# TODO use SAF,VAF, and abbrevs for ScalarAffineTerm etc
function MOI.optimize!(opt::Optimizer)
    solver.status = :SolveCalled
    start_time = time()

    # add initial fixed OA cuts eg variable bounds
    for k in eachindex(solver.con_sets)
        con_set = solver.con_sets[k]
        con_fun = solver.con_funs[k]

        cuts = Cuts.get_init_cuts(con_set)
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
            cut_ref = MOI.add_constraint(solver.oa_model, cut_expr, MOI.GreaterThan(-cut_constant))
            # TODO store cut_ref with the constraint so can access values/duals
        end

        solver.num_cuts += length(cuts)
    end

    if solver.verbose
        println("initial cut count: $(solver.num_cuts)")
        @printf("\n%5s %8s %12s\n",
            "iter", "cuts", "obj bound",
            )
        flush(stdout)
    end

    while true
        MOI.optimize!(solver.oa_model)

        oa_model_status = MOI.get(solver.oa_model, MOI.TerminationStatus())
        if oa_model_status == MOI.OPTIMAL
            solver.obj_bound = MOI.get(solver.oa_model, MOI.ObjectiveBound())
        elseif oa_model_status == MOI.DUAL_INFEASIBLE
            # TODO need to cut off the ray
        elseif oa_model_status == MOI.INFEASIBLE
            solver.verbose && println("infeasibility detected; terminating")
            solver.status = :PrimalInfeasible
            break
        elseif oa_model_status == MOI.INFEASIBLE_OR_UNBOUNDED
            @show oa_model_status
            # solver.verbose && println("infeasibility or unboundedness detected; terminating")
            # solver.status = :ApproxAlgorithmStatusNotHandled
            solver.verbose && println("LP solver returned infeasible or unbounded - we are assuming infeasibility, but we could be wrong")
            solver.status = :PrimalInfeasible
            break
        else
            error("OA solver status not handled")
        end

        if solver.verbose
            @printf("%5d %8d %12.4e\n",
                solver.num_iters, solver.num_cuts, solver.obj_bound,
                )
            flush(stdout)
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
        for k in eachindex(solver.con_sets)
            con_set = solver.con_sets[k]
            con_fun = solver.con_funs[k]
            # TODO do not include constant if primal solution is a ray?
            fun_val = MOIU.evalvariables(vi -> MOI.get(solver.oa_model, MOI.VariablePrimal(), vi), con_fun)
            @assert !any(isnan, fun_val)

            cuts = Cuts.get_sep_cuts(fun_val, con_set)
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
                cut_ref = MOI.add_constraint(solver.oa_model, cut_expr, MOI.GreaterThan(-cut_constant))
                # TODO store cut_ref with the constraint so can access values/duals
                is_cut_off = true
            end

            solver.num_cuts += length(cuts)
        end
        if !is_cut_off
            solver.verbose && println("no cuts were added; terminating")
            solver.status = :Optimal
            solver.obj_value = MOI.get(solver.oa_model, MOI.ObjectiveValue())
            break
        end

        solver.num_iters += 1
    end

    solver.solve_time = time() - start_time

    solver.verbose && println("\nstatus is $(solver.status) after $(solver.num_iters) iterations and $(trunc(solver.solve_time, digits=3)) seconds\n")

    return
end
