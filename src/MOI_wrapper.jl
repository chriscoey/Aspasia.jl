# MathOptInterface wrapper of Aspasia solver

const DiscreteSet = Union{MOI.Integer, MOI.ZeroOne, MOI.Semiinteger{Float64}}

function _oa_opt(opt::Optimizer)
    if isnothing(opt.oa_opt)
        if isnothing(opt.oa_solver)
            error("No MIP solver specified (set `oa_solver` attribute)")
        end
        opt.oa_opt = MOI.instantiate(opt.oa_solver, with_bridge_type = Float64)
        if !MOI.supports(opt.oa_opt, MOI.LazyConstraintCallback())
            error("MIP solver does not support lazy constraint callbacks")
        end
    end
    return opt.oa_opt
end

MOI.is_empty(opt::Optimizer) = (isnothing(opt.oa_opt) || MOI.is_empty(opt.oa_opt))

MOI.empty!(opt::Optimizer) = _empty(opt)

MOI.get(::Optimizer, ::MOI.SolverName) = "Aspasia"

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(opt::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(opt, src)
end

# TODO is this OK?
# function MOI.get(opt::Optimizer, attr::MOI.AbstractModelAttribute)
#     return MOI.get(_oa_opt(opt), attr)
# end
function MOI.get(opt::Optimizer, attr::MOI.AnyAttribute, args...)
    return MOI.get(_oa_opt(opt), attr, args...)
end
function MOI.get(model::Optimizer, attr::MOI.AnyAttribute, idxs::Vector)
    return MOI.get.(model, attr, idxs)
end

function MOI.supports(opt::Optimizer, attr::MOI.AnyAttribute, args...)
    return MOI.supports(_oa_opt(opt), attr, args...)
end

function MOI.set(opt::Optimizer, attr::MOI.AnyAttribute, args...)
    return MOI.set(_oa_opt(opt), attr, args...)
end

function MOI.supports_constraint(
    opt::Optimizer,
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet},
)
    return MOI.supports_constraint(_oa_opt(opt), F, S)
end

# NOTE slack bridge should handle VAF in cone
# MOI.supports_constraint(::Optimizer, ::Type{<:Union{VV, VAF}}, ::Type{<:ApproxSet}) = true
MOI.supports_constraint(::Optimizer, ::Type{VV}, ::Type{<:ApproxSet}) = true

MOI.add_variable(opt::Optimizer) = MOI.add_variable(_oa_opt(opt))

function MOI.add_constraint(
    opt::Optimizer,
    func::F,
    set::S,
) where {F <: MathOptInterface.AbstractFunction, S <: MathOptInterface.AbstractSet}
    if S <: DiscreteSet
        opt.has_discrete_constraint = true
    end
    if MOI.supports_constraint(_oa_opt(opt), F, S)
        return MOI.add_constraint(_oa_opt(opt), func, set)
    end
    # store for outer approximation
    push!(opt.con_vars, func.variables)
    push!(opt.con_sets, set)
    return MOI.ConstraintIndex{F, S}(length(opt.con_vars))
end

function MOI.get(opt::Optimizer, param::MOI.RawOptimizerAttribute)
    return getproperty(opt, Symbol(param.name))
end

function MOI.supports(::Optimizer, param::MOI.RawOptimizerAttribute)
    return (Symbol(param.name) in fieldnames(Optimizer))
end

function MOI.set(opt::Optimizer, param::MOI.RawOptimizerAttribute, value)
    return setproperty!(opt, Symbol(param.name), value)
end

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(opt::Optimizer, ::MOI.Silent, value::Bool)
    # MOI.set(_oa_opt(opt), MOI.Silent, value)
    opt.verbose = !value
    return
end

MOI.get(opt::Optimizer, ::MOI.Silent) = !opt.verbose

MOI.get(::Optimizer, ::MOI.DualStatus) = MOI.NO_SOLUTION
