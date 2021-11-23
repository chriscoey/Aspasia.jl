# Outer approximation algorithm

mutable struct Optimizer <: MOI.AbstractOptimizer
    verbose::Bool
    tol_feas::Float64
    oa_solver::Union{Nothing, MOI.OptimizerWithAttributes}

    oa_model::Union{Nothing, MOI.ModelLike}
    con_vars::Vector{Vector{MOI.VariableIndex}}
    con_sets::Vector{MOI.AbstractVectorSet}

    num_callbacks::Int
    num_cuts::Int

    function Optimizer()
        opt = new()
        opt.verbose = true
        opt.tol_feas = 1e-6 # TODO??
        opt.oa_solver = nothing
        opt.oa_model = nothing
        opt.con_vars = Vector{MOI.VariableIndex}[]
        opt.con_sets = MOI.AbstractVectorSet[]
        opt.num_callbacks = 0
        opt.num_cuts = 0
        return opt
    end
end

function MOI.optimize!(opt::Optimizer)
    for k in eachindex(opt.con_sets)
        add_init_cuts(k, solver)
    end
    if opt.verbose
        println("added $(opt.num_cuts) initial cuts; starting MIP solver\n")
    end

    function lazy_callback(cb)
        opt.num_callbacks += 1
        sol = MOI.get(opt.oa_model, MOI.CallbackVariablePrimal(cb))

        is_cut_off = any(add_cb_cuts(k, cb, solver) for k in eachindex(opt.con_sets))
        if !is_cut_off && opt.verbose
            println("no cuts were added during callback")
        end
    end
    MOI.set(opt.oa_model, MOI.LazyConstraintCallback(), lazy_callback)

    MOI.optimize!(opt.oa_model)

    if opt.verbose
        mip_status = MOI.get(opt.oa_model, MOI.TerminationStatus())
        println(
            "MIP solver finished with status $(mip_status), after " *
            "$(opt.num_callbacks) callbacks and $(opt.num_cuts) cuts\n",
        )
    end

    return
end

function add_init_cuts(k::Int, opt::Optimizer)
    con_vars = opt.con_vars[k]
    con_set = opt.con_sets[k]

    cuts = Cuts.get_init_cuts(con_set)
    for cut in cuts
        cut_expr = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(cut, con_vars), 0.0)
        MOI.add_constraint(opt.oa_model, cut_expr, MOI.GreaterThan(0.0))
    end

    opt.num_cuts += length(cuts)
    return nothing
end

function add_cb_cuts(k::Int, cb, opt::Optimizer)
    con_vars = opt.con_vars[k]
    con_set = opt.con_sets[k]
    var_vals = MOI.get(opt.oa_model, MOI.VariablePrimal(), con_vars)
    @assert !any(isnan, var_vals)

    cuts = Cuts.get_sep_cuts(var_vals, con_set)
    for cut in cuts
        cut_expr = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(cut, con_vars), 0.0)
        MOI.add_constraint(opt.oa_model, cut_expr, MOI.GreaterThan(0.0))
    end

    opt.num_cuts += length(cuts)
    return !isempty(cuts)
end
