# Outer approximation algorithm

mutable struct Optimizer <: MOI.AbstractOptimizer
    verbose::Bool
    tol_feas::Float64
    oa_solver::Union{Nothing, MOI.OptimizerWithAttributes}
    use_iterative_method::Bool

    oa_opt::Union{Nothing, MOI.ModelLike}
    con_vars::Vector{Vector{VI}}
    con_sets::Vector{MOI.AbstractVectorSet}
    num_cuts::Int

    # for iterative method
    num_iters::Int
    status::Symbol

    # for OA-solver-driven method
    has_discrete_constraint::Bool
    dummy_discrete_constraint::Union{Nothing, Tuple{VI, CI}}
    num_callbacks::Int

    function Optimizer()
        opt = new()
        opt.verbose = true
        opt.tol_feas = 1e-7
        opt.oa_solver = nothing
        opt.use_iterative_method = false
        return _empty(opt)
    end
end

function _empty(opt::Optimizer)
    opt.oa_opt = nothing
    opt.con_vars = Vector{VI}[]
    opt.con_sets = MOI.AbstractVectorSet[]
    opt.num_cuts = 0
    opt.num_iters = 0
    opt.status = :Empty
    opt.has_discrete_constraint = false
    opt.dummy_discrete_constraint = nothing
    opt.num_callbacks = 0
    return opt
end

function MOI.optimize!(opt::Optimizer)
    # add initial cuts
    for k in eachindex(opt.con_sets)
        add_init_cuts(k, opt)
    end
    if opt.verbose
        println("added $(opt.num_cuts) initial cuts; starting OA solver")
    end

    if opt.use_iterative_method
        iterative_method(opt)
    else
        oa_solver_driven_method(opt)
    end

    if opt.verbose
        oa_status = MOI.get(opt.oa_opt, MOI.TerminationStatus())
        println("OA solver finished with status $oa_status, after $(opt.num_cuts) cuts")
        if opt.use_iterative_method
            println("iterative method used $(opt.num_iters) iterations\n")
        else
            println("OA solver driven method used $(opt.num_callbacks) callbacks\n")
        end
    end
    return
end

# TODO handle statuses properly
function iterative_method(opt::Optimizer)
    while true
        MOI.optimize!(opt.oa_opt)

        oa_opt_status = MOI.get(opt.oa_opt, MOI.TerminationStatus())
        if oa_opt_status == MOI.OPTIMAL
            obj_bound = MOI.get(opt.oa_opt, MOI.ObjectiveBound())
            if opt.verbose
                @printf("%5d %8d %12.4e\n", opt.num_iters, opt.num_cuts, obj_bound)
            end
        elseif oa_opt_status == MOI.DUAL_INFEASIBLE
            if opt.verbose
                println("OA solver status is $oa_opt_status; relaxation likely unbounded")
            end
            opt.status = :DualInfeasibleRelaxation
            break
        elseif oa_opt_status == MOI.INFEASIBLE
            if opt.verbose
                println("infeasibility detected; terminating")
            end
            opt.status = :PrimalInfeasible
            break
        elseif oa_opt_status == MOI.INFEASIBLE_OR_UNBOUNDED
            if opt.verbose
                println("OA solver status is $oa_opt_status; assuming infeasibility")
            end
            opt.status = :PrimalInfeasible
            break
        else
            error("OA solver status $oa_opt_status not handled")
        end

        # TODO
        # if opt.num_iters == opt.max_iters
        #     opt.verbose && println("iteration limit reached; terminating")
        #     opt.status = :IterationLimit
        #     break
        # end

        # if time() - start_time >= opt.time_limit
        #     opt.verbose && println("time limit reached; terminating")
        #     opt.status = :TimeLimit
        #     break
        # end

        is_cut_off = any(add_sep_cuts(k, opt) for k in eachindex(opt.con_sets))
        if !is_cut_off
            if opt.verbose
                println("no cuts were added; terminating")
            end
            opt.status = :Optimal
            break
        end

        opt.num_iters += 1
        flush(stdout)
    end
    return
end

function oa_solver_driven_method(opt::Optimizer)
    if !opt.has_discrete_constraint
        # no discrete constraints, so add one to force OA solver to use lazy callbacks
        if opt.verbose
            println("model has no discrete constraints; adding a dummy discrete constraint")
        end
        opt.dummy_discrete_constraint = MOI.add_constrained_variable(opt, MOI.Integer())
    end

    function lazy_callback(cb)
        opt.num_callbacks += 1
        is_cut_off = any(add_sep_cuts(k, opt, cb) for k in eachindex(opt.con_sets))
        if !is_cut_off && opt.verbose
            println("no cuts were added during callback")
        end
    end
    MOI.set(opt.oa_opt, MOI.LazyConstraintCallback(), lazy_callback)

    MOI.optimize!(opt.oa_opt)

    if !isnothing(opt.dummy_discrete_constraint)
        # delete dummy discrete constraint
        (vi, ci) = opt.dummy_discrete_constraint
        MOI.delete(opt.oa_opt, ci)
        MOI.delete(opt.oa_opt, vi)
        opt.dummy_discrete_constraint = nothing
    end

    return
end

function add_init_cuts(k::Int, opt::Optimizer)
    con_vars = opt.con_vars[k]
    con_set = opt.con_sets[k]

    cuts = Cuts.get_init_cuts(con_set)
    for cut in cuts
        cut_expr = SAF(MOI.ScalarAffineTerm.(cut, con_vars), 0.0)
        MOI.add_constraint(opt.oa_opt, cut_expr, MOI.GreaterThan(0.0))
    end

    opt.num_cuts += length(cuts)
    return nothing
end

function add_sep_cuts(k::Int, opt::Optimizer, cb = nothing)
    vars = opt.con_vars[k]
    point = _get_var_primal(vars, opt.oa_opt, cb)
    @assert !any(isnan, point)

    cuts = Cuts.get_sep_cuts(point, opt.con_sets[k], opt.tol_feas)
    for cut in cuts
        cut_expr = SAF(MOI.ScalarAffineTerm.(cut, vars), 0.0)
        _add_constraint(cut_expr, opt.oa_opt, cb)
    end

    opt.num_cuts += length(cuts)
    return !isempty(cuts)
end

function _add_constraint(cut_expr::SAF, oa_opt::MOI.ModelLike, ::Nothing)
    return MOI.add_constraint(oa_opt, cut_expr, MOI.GreaterThan(0.0))
end

function _add_constraint(cut_expr::SAF, oa_opt::MOI.ModelLike, cb)
    return MOI.submit(oa_opt, MOI.LazyConstraint(cb), cut_expr, MOI.GreaterThan(0.0))
end

function _get_var_primal(vars::Vector{VI}, oa_opt::MOI.ModelLike, ::Nothing)
    return MOI.get(oa_opt, MOI.VariablePrimal(), vars)
end

function _get_var_primal(vars::Vector{VI}, oa_opt::MOI.ModelLike, cb)
    return MOI.get(oa_opt, MOI.CallbackVariablePrimal(cb), vars)
end
