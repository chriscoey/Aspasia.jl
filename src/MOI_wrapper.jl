# MathOptInterface wrapper of Aspasia solver

# supported sets needing outer approximation
const ApproxSet = Union{MOI.SecondOrderCone, MOI.PositiveSemidefiniteConeTriangle}

function _oa_model(model::Optimizer)
    if isnothing(model.oa_model)
        if isnothing(model.oa_solver)
            error("No MIP solver specified (set `oa_solver` attribute)")
        end
        model.oa_model = MOI.instantiate(model.oa_solver, with_bridge_type = Float64)
        if !MOI.supports(model.oa_model, MOI.LazyConstraintCallback())
            error("MIP solver does not support lazy constraint callbacks")
        end
    end
    return model.oa_model
end

MOI.is_empty(model::Optimizer) = (isnothing(model.oa_model) || MOI.is_empty(model.oa_model))

function MOI.empty!(model::Optimizer)
    model.oa_model = nothing
    # TODO etc
    # model.status = MOI.OPTIMIZE_NOT_CALLED
    return
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Aspasia"

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(model, src)
end

# TODO is this OK?
function MOI.get(model::Optimizer, attr::MOI.AbstractModelAttribute)
    return MOI.get(_oa_model(model), attr)
end

function MOI.supports(model::Optimizer, attr::MOI.AnyAttribute, args...)
    return MOI.supports(_oa_model(model), attr, args...)
end

function MOI.set(model::Optimizer, attr::MOI.AnyAttribute, args...)
    return MOI.set(_oa_model(model), attr, args...)
end

function MOI.supports_constraint(
    model::Optimizer,
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet},
)
    return MOI.supports_constraint(_oa_model(model), F, S)
end

MOI.supports_constraint(::Optimizer, ::Type{<:Union{VV, VAF}}, ::Type{<:ApproxSet}) = true

MOI.add_variable(model::Optimizer) = MOI.add_variable(_oa_model(model))

function MOI.add_constraint(
    model::Optimizer,
    func::F,
    set::S,
) where {F <: MathOptInterface.AbstractFunction, S <: MathOptInterface.AbstractSet}
    if MOI.supports_constraint(_oa_model(model), F, S)
        return MOI.add_constraint(_oa_model(model), func, set)
    end
    # store for outer approximation
    push!(model.con_funs, func)
    push!(model.con_sets, set)
    println("\napprox set\n")
    return #TODO
end

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(model, Symbol(param.name))
end

function MOI.supports(::Optimizer, param::MOI.RawOptimizerAttribute)
    return (Symbol(param.name) in fieldnames(Optimizer))
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(model, Symbol(param.name), value)
end

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    # MOI.set(_oa_model(model), MOI.Silent, value)
    model.verbose = !value
    return
end

MOI.get(model::Optimizer, ::MOI.Silent) = !model.verbose

MOI.get(::Optimizer, ::MOI.DualStatus) = MOI.NO_SOLUTION
