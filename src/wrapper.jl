#=
MathOptInterface wrapper of Aspasia solver
=#

const VI = MOI.VariableIndex
const SAF = MOI.ScalarAffineFunction{Float64}
const VV = MOI.VectorOfVariables
const VAF = MOI.VectorAffineFunction{Float64}

# sets not needing polyhedral approximation
const LinearScalarSet = Union{
    MOI.EqualTo{Float64},
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.Interval{Float64},
    }
const LinearVectorSet = Union{
    MOI.Zeros,
    MOI.Nonnegatives,
    MOI.Nonpositives,
    }
const LinearSet = Union{LinearScalarSet, LinearVectorSet}

# cones needing polyhedral approximation
const NonlinearSet = Union{
    MOI.SecondOrderCone,
    MOI.PositiveSemidefiniteConeTriangle,
    # MOI.RotatedSecondOrderCone,
    # MOI.ExponentialCone,
    # MOI.PowerCone{Float64},
    # MOI.GeometricMeanCone,
    }

export Optimizer

"""
$(TYPEDEF)

A MathOptInterface optimizer type for Aspasia.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    verbose::Bool
    tol_rel_opt::Float64
    tol_abs_opt::Float64
    tol_feas::Float64
    max_iters::Int
    time_limit::Float64

    oa_solver
    oa_solver_options

    oa_model::MOI.ModelLike
    con_sets::Vector{MOI.AbstractVectorSet}
    con_funs::Vector{MOI.AbstractVectorFunction}

    status::Symbol
    num_iters::Int
    num_cuts::Int
    solve_time::Float64
    obj_value::Float64
    obj_bound::Float64

    function Optimizer(
        oa_solver,
        oa_solver_options,
        verbose::Bool,
        max_iters::Int,
        time_limit::Float64,
        tol_rel_opt::Float64,
        tol_abs_opt::Float64,
        tol_feas::Float64,
        )
        opt = new()

        opt.oa_solver = oa_solver
        opt.oa_solver_options = oa_solver_options

        opt.verbose = verbose
        opt.max_iters = max_iters
        opt.time_limit = time_limit
        opt.tol_rel_opt = tol_rel_opt
        opt.tol_abs_opt = tol_abs_opt
        opt.tol_feas = tol_feas

        opt.status = :NotLoaded
        opt.num_iters = -1
        opt.num_cuts = -1
        opt.solve_time = NaN
        opt.obj_value = NaN
        opt.obj_bound = NaN

        return opt
    end
end

Optimizer(
    oa_solver,
    oa_solver_options;
    verbose::Bool = true,
    max_iters::Int = 500,
    time_limit::Float64 = 3.6e3, # TODO should be Inf
    tol_rel_opt::Float64 = 1e-4,
    tol_abs_opt::Float64 = 1e-4,
    tol_feas::Float64 = 1e-4,
    ) = Optimizer(oa_solver, oa_solver_options, verbose, max_iters, time_limit,
        tol_rel_opt, tol_abs_opt, tol_feas)

MOI.get(::Optimizer, ::MOI.SolverName) = "Aspasia"

MOI.is_empty(opt::Optimizer) = (opt.status == :NotLoaded)

MOI.empty!(opt::Optimizer) = (opt.status = :NotLoaded)

MOI.supports(
    ::Optimizer,
    ::Union{MOI.ObjectiveSense, MOI.ObjectiveFunction{<:Union{VI, SAF}}},
    ) = true

MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:Union{VI, SAF}},
    ::Type{<:LinearScalarSet},
    ) = true

MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:Union{VV, VAF}},
    ::Type{<:Union{LinearScalarSet, NonlinearSet}},
    ) = true

# build polyhedral approximation MOI model and list of cone objects
function MOI.copy_to(
    opt::Optimizer,
    src::MOI.ModelLike;
    copy_names::Bool = false,
    )
    @assert !copy_names
    idx_map = MOIU.IndexMap()

    oa_model = opt.oa_solver(; opt.oa_solver_options...)

    # variables
    n = MOI.get(src, MOI.NumberOfVariables()) # columns of A
    x = MOI.add_variables(oa_model, n)
    for (vj, xj) in zip(MOI.get(src, MOI.ListOfVariableIndices()), x)
        idx_map[vj] = xj
    end

    # objective function
    obj_type = MOI.get(src, MOI.ObjectiveFunctionType())
    obj = MOI.get(src, MOI.ObjectiveFunction{obj_type}())
    MOI.set(oa_model, MOI.ObjectiveFunction{obj_type}(), obj)

    obj_sense = MOI.get(src, MOI.ObjectiveSense())
    MOI.set(oa_model, MOI.ObjectiveSense(), obj_sense)

    # constraints
    con_sets = MOI.AbstractVectorSet[]
    con_funs = MOI.AbstractVectorFunction[]


# TODO use filter
# https://github.com/jump-dev/MathOptInterface.jl/blob/0431632728feb2bde3bbd49bcc2a6004abdab632/src/Utilities/copy.jl#L438
#
    i = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        if S <: LinearSet
            # equality and orthant cone constraints
            MOIU.copyconstraints!(oa_model, src, false, idx_map, F, S)
        else
            # constraints requiring polyhedral approximation
            for ci in MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
                i += 1
                idx_map[ci] = MOI.ConstraintIndex{F, S}(i)
                si = MOI.get(src, MOI.ConstraintSet(), ci)
                fi = MOI.get(src, MOI.ConstraintFunction(), ci)
                push!(con_sets, si)
                push!(con_funs, fi)
            end
        end
    end

    opt.oa_model = oa_model
    opt.con_sets = con_sets
    opt.con_funs = con_funs
    opt.status = :SolveNotCalled

    return idx_map
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
MOI.set(opt::Optimizer, ::MOI.Silent, value::Bool) = (opt.verbose = !value)
MOI.get(opt::Optimizer, ::MOI.Silent) = !opt.verbose

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.set(opt::Optimizer, ::MOI.TimeLimitSec, value::Real) =
    (opt.time_limit = value)
MOI.set(opt::Optimizer, ::MOI.TimeLimitSec, ::Nothing) =
    (opt.time_limit = Inf)

function MOI.get(opt::Optimizer, ::MOI.TimeLimitSec)
    isfinite(opt.time_limit) && return opt.time_limit
    return
end

function MOI.get(opt::Optimizer, ::MOI.SolveTimeSec)
    if opt.status in (:NotLoaded, :Loaded)
        error("solve has not been called")
    end
    return opt.solve_time
end

MOI.get(opt::Optimizer, ::MOI.RawStatusString) = string(opt.status)

function MOI.get(opt::Optimizer, ::MOI.TerminationStatus)
    if opt.status in (:NotLoaded, :Loaded)
        return MOI.OPTIMIZE_NOT_CALLED
    elseif opt.status == :Optimal
        return MOI.OPTIMAL
    elseif opt.status == :PrimalInfeasible
        return MOI.INFEASIBLE
    elseif opt.status == :DualInfeasible
        return MOI.DUAL_INFEASIBLE
    elseif opt.status == :IterationLimit
        return MOI.ITERATION_LIMIT
    elseif opt.status == :TimeLimit
        return MOI.TIME_LIMIT
    else
        @warn("Aspasia status $(opt.status) not handled")
        return MOI.OTHER_ERROR
    end
end

MOI.get(opt::Optimizer, ::MOI.ObjectiveValue) = opt.obj_value
MOI.get(opt::Optimizer, ::MOI.ObjectiveBound) = opt.obj_bound

function MOI.get(opt::Optimizer, ::MOI.ResultCount)
    if opt.status in (:Optimal, :PrimalInfeasible, :DualInfeasible)
        return 1
    end
    return 0
end

function MOI.get(opt::Optimizer, ::MOI.PrimalStatus)
    if opt.status == :Optimal
        return MOI.FEASIBLE_POINT
    elseif opt.status == :PrimalInfeasible
        return MOI.INFEASIBLE_POINT
    elseif opt.status == :DualInfeasible
        return MOI.INFEASIBILITY_CERTIFICATE
    else
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

function MOI.get(opt::Optimizer, ::MOI.DualStatus)
    if opt.status == :Optimal
        return MOI.FEASIBLE_POINT
    elseif opt.status == :PrimalInfeasible
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif opt.status == :DualInfeasible
        return MOI.INFEASIBLE_POINT
    else
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

MOI.get(
    opt::Optimizer,
    ::MOI.VariablePrimal,
    vi::VI,
    ) = MOI.get(opt.oa_model, MOI.VariablePrimal(), vi)

MOI.get(
    opt::Optimizer,
    a::MOI.VariablePrimal,
    vi::Vector{VI},
    ) = MOI.get.(opt, a, vi)

MOI.get(
    opt::Optimizer,
    ::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{F, S}
    ) where {F <: MOI.AbstractFunction, S <: LinearSet} =
    MOI.get(opt.oa_model, MOI.ConstraintDual(), ci)

# function MOI.get(
#     opt::Optimizer,
#     ::MOI.ConstraintDual,
#     ci::MOI.ConstraintIndex{F, S},
#     ) where {F <: MOI.AbstractFunction, S <: NonlinearSet}
#     # dual vector corresponding to the constraint's polyhedral approximation is combination of the cuts weighted by their duals
#     # TODO use the refs of the cuts for the constraint
#     return []
# end

MOI.get(
    opt::Optimizer,
    a::MOI.ConstraintDual,
    ci::Vector{MOI.ConstraintIndex},
    ) = MOI.get.(opt, a, ci)

MOI.get(
    opt::Optimizer,
    ::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{F, S},
    ) where {F <: MOI.AbstractFunction, S <: LinearSet} =
    MOI.get(opt.oa_model, MOI.ConstraintPrimal(), ci)

MOI.get(
    opt::Optimizer,
    ::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{F, S},
    ) where {F <: MOI.AbstractFunction, S <: NonlinearSet} =
    MOIU.evalvariables(vi -> MOI.get(opt.oa_model, MOI.VariablePrimal(), vi),
    opt.con_funs[ci.value])

MOI.get(
    opt::Optimizer,
    a::MOI.ConstraintPrimal,
    ci::Vector{MOI.ConstraintIndex},
    ) = MOI.get.(opt, a, ci)
