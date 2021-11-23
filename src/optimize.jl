#=
outer approximation algorithm
=#

"""
$(TYPEDEF)

A MathOptInterface optimizer type for Aspasia.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    verbose::Bool
    tol_feas::Float64
    oa_solver::Union{Nothing, MOI.OptimizerWithAttributes}

    oa_model::Union{Nothing, MOI.ModelLike}
    con_funs::Vector{MOI.AbstractVectorFunction}
    con_sets::Vector{MOI.AbstractVectorSet}

    num_callbacks::Int
    num_cuts::Int

    function Optimizer()
        opt = new()
        opt.verbose = true
        opt.tol_feas = 1e-6 # TODO??
        opt.oa_solver = nothing
        opt.oa_model = nothing
        opt.con_funs = MOI.AbstractVectorFunction[]
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
    con_set = opt.con_sets[k]
    con_fun = opt.con_funs[k]

    cuts = Cuts.get_init_cuts(con_set)
    for cut in cuts
        if con_fun isa MOI.VectorOfVariables
            cut_terms = MOI.ScalarAffineTerm.(cut, con_fun.variables)
            cut_constant = 0.0
        else
            cut_terms = [
                MOI.ScalarAffineTerm(
                    cut[t.output_index] * t.scalar_term.coefficient,
                    t.scalar_term.variable_index,
                ) for t in con_fun.terms
            ]
            cut_constant = dot(cut, con_fun.constants)
        end

        cut_expr = MOI.ScalarAffineFunction(cut_terms, 0.0)
        cut_ref = MOI.add_constraint(opt.oa_model, cut_expr, MOI.GreaterThan(-cut_constant))
    end

    opt.num_cuts += length(cuts)
    return nothing
end

function add_cb_cuts(k::Int, cb, opt::Optimizer)
    con_set = opt.con_sets[k]
    con_fun = opt.con_funs[k]
    # TODO do not include constant if primal solution is a ray?
    fun_val =
        MOIU.evalvariables(vi -> MOI.get(opt.oa_model, MOI.VariablePrimal(), vi), con_fun)
    @assert !any(isnan, fun_val)

    cuts = Cuts.get_sep_cuts(fun_val, con_set)
    is_cut_off = false
    for cut in cuts
        if con_fun isa MOI.VectorOfVariables
            cut_terms = MOI.ScalarAffineTerm.(cut, con_fun.variables)
            cut_constant = 0.0
        else
            cut_terms = [
                MOI.ScalarAffineTerm(
                    cut[t.output_index] * t.scalar_term.coefficient,
                    t.scalar_term.variable_index,
                ) for t in con_fun.terms
            ]
            cut_constant = dot(cut, con_fun.constants)
        end

        cut_expr = MOI.ScalarAffineFunction(cut_terms, 0.0)
        cut_ref = MOI.add_constraint(opt.oa_model, cut_expr, MOI.GreaterThan(-cut_constant))
        # TODO store cut_ref with the constraint so can access values/duals
        is_cut_off = true
    end

    opt.num_cuts += length(cuts)
    return is_cut_off
end
