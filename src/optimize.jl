# Outer approximation algorithm

mutable struct Optimizer <: MOI.AbstractOptimizer
    verbose::Bool
    tol_feas::Float64
    oa_solver::Union{Nothing, MOI.OptimizerWithAttributes}

    oa_opt::Union{Nothing, MOI.ModelLike}
    con_vars::Vector{Vector{MOI.VariableIndex}}
    con_sets::Vector{MOI.AbstractVectorSet}
    has_discrete_constraint::Bool
    dummy_discrete_constraint::Union{Nothing, Tuple{VI, CI}}

    num_callbacks::Int
    num_cuts::Int

    function Optimizer()
        opt = new()
        opt.verbose = true
        opt.tol_feas = 1e-7
        opt.oa_solver = nothing
        return _empty(opt)
    end
end

function _empty(opt::Optimizer)
    opt.oa_opt = nothing
    opt.con_vars = Vector{MOI.VariableIndex}[]
    opt.con_sets = MOI.AbstractVectorSet[]
    opt.has_discrete_constraint = false
    opt.dummy_discrete_constraint = nothing
    opt.num_callbacks = 0
    opt.num_cuts = 0
    return opt
end

function MOI.optimize!(opt::Optimizer)
    if !opt.has_discrete_constraint
        # no discrete constraints, so add one to force MIP solver to use lazy callbacks
        if opt.verbose
            println("model has no discrete constraints; adding a dummy discrete constraint")
        end
        opt.dummy_discrete_constraint = MOI.add_constrained_variable(opt, MOI.Integer())
    end

    for k in eachindex(opt.con_sets)
        add_init_cuts(k, opt)
    end
    if opt.verbose
        println("added $(opt.num_cuts) initial cuts; starting MIP solver\n")
    end

    function lazy_callback(cb)
        opt.num_callbacks += 1
        is_cut_off = any(add_cb_cuts(k, cb, opt) for k in eachindex(opt.con_sets))
        if !is_cut_off && opt.verbose
            println("no cuts were added during callback")
        end
    end
    MOI.set(opt.oa_opt, MOI.LazyConstraintCallback(), lazy_callback)

    MOI.optimize!(opt.oa_opt)

    if opt.verbose
        mip_status = MOI.get(opt.oa_opt, MOI.TerminationStatus())
        println(
            "MIP solver finished with status $(mip_status), after " *
            "$(opt.num_callbacks) callbacks and $(opt.num_cuts) cuts\n",
        )
    end

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
        cut_expr = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(cut, con_vars), 0.0)
        MOI.add_constraint(opt.oa_opt, cut_expr, MOI.GreaterThan(0.0))
    end

    opt.num_cuts += length(cuts)
    return nothing
end

function add_cb_cuts(k::Int, cb, opt::Optimizer)
    con_vars = opt.con_vars[k]
    con_set = opt.con_sets[k]
    var_vals = MOI.get(opt.oa_opt, MOI.CallbackVariablePrimal(cb), con_vars)
    @assert !any(isnan, var_vals)

    cuts = Cuts.get_sep_cuts(var_vals, con_set, opt.tol_feas)
    for cut in cuts
        cut_expr = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(cut, con_vars), 0.0)
        MOI.add_constraint(opt.oa_opt, cut_expr, MOI.GreaterThan(0.0))
    end

    opt.num_cuts += length(cuts)
    return !isempty(cuts)
end
