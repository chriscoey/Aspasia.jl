#=
Copyright 2019, Chris Coey and contributors

MathOptInterface wrapper of Aspasia solver
=#

mutable struct Optimizer <: MOI.AbstractOptimizer
    verbose::Bool
    tol_rel_opt::Float64
    tol_abs_opt::Float64
    tol_feas::Float64
    max_iters::Int
    time_limit::Float64

    approx_solver
    approx_options

    approx_model::MOI.ModelLike
    con_sets::Vector{MOI.AbstractVectorSet}
    con_funs::Vector{MOI.AbstractVectorFunction}

    idx_map::MOIU.IndexMap

    status::Symbol
    solve_time::Float64
    obj_value::Float64
    obj_bound::Float64

    function Optimizer(approx_solver, approx_options, verbose::Bool, max_iters::Int, time_limit::Float64, tol_rel_opt::Float64, tol_abs_opt::Float64, tol_feas::Float64)
        opt = new()

        opt.approx_solver = approx_solver
        opt.approx_options = approx_options

        opt.verbose = verbose
        opt.max_iters = max_iters
        opt.time_limit = time_limit
        opt.tol_rel_opt = tol_rel_opt
        opt.tol_abs_opt = tol_abs_opt
        opt.tol_feas = tol_feas

        opt.status = :NotLoaded
        opt.solve_time = NaN
        opt.obj_value = NaN
        opt.obj_bound = NaN

        return opt
    end
end

Optimizer(
    approx_solver,
    approx_options;
    verbose::Bool = false,
    max_iters::Int = 500,
    time_limit::Float64 = 3.6e3, # TODO should be Inf
    tol_rel_opt::Float64 = 1e-6,
    tol_abs_opt::Float64 = 1e-7,
    tol_feas::Float64 = 1e-7,
    ) = Optimizer(approx_solver, approx_options, verbose, max_iters, time_limit, tol_rel_opt, tol_abs_opt, tol_feas)

MOI.get(::Optimizer, ::MOI.SolverName) = "Aspasia"

MOI.is_empty(opt::Optimizer) = (opt.status == :NotLoaded)
MOI.empty!(opt::Optimizer) = (opt.status = :NotLoaded) # TODO empty the data and results? keep options?

# TODO support convex quad obj if OA solver does
MOI.supports(::Optimizer, ::Union{
    MOI.ObjectiveSense,
    MOI.ObjectiveFunction{MOI.SingleVariable},
    MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}},
    }) = true

# TODO don't restrict to Float64 type
# TODO allow SOC/quadratic constraints if approx solver allows specifying?
SupportedFun = Union{
    MOI.SingleVariable, MOI.ScalarAffineFunction{Float64},
    MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64},
    }

LinearSet = Union{
    MOI.EqualTo{Float64},
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.Interval{Float64},
    MOI.Zeros,
    MOI.Nonnegatives,
    MOI.Nonpositives
    }

NonlinearSet = Union{ # cones to use polyhedral approximation on
    # MOI.SecondOrderCone,
    # MOI.RotatedSecondOrderCone,
    # MOI.ExponentialCone,
    # MOI.PowerCone{Float64},
    # MOI.GeometricMeanCone,
    MOI.PositiveSemidefiniteConeTriangle,
    # MOI.LogDetConeTriangle,
    WSOSPolyInterpCone,
    WSOSPolyInterpMatCone,
    }

MOI.supports_constraint(::Optimizer, ::Type{<:SupportedFun}, ::Type{<:LinearSet}) = true
MOI.supports_constraint(::Optimizer, ::Type{<:SupportedFun}, ::Type{<:NonlinearSet}) = true

# build polyhedral approximation MOI model and list of cone objects
function MOI.copy_to(
    opt::Optimizer,
    src::MOI.ModelLike;
    copy_names::Bool = false,
    warn_attributes::Bool = true,
    )
    @assert !copy_names
    idx_map = MOIU.IndexMap()

    approx_model = opt.approx_solver(; opt.approx_options...)

    # variables
    n = MOI.get(src, MOI.NumberOfVariables()) # columns of A
    x = MOI.add_variables(approx_model, n)
    for (vj, xj) in zip(MOI.get(src, MOI.ListOfVariableIndices()), x)
        idx_map[vj] = xj
    end

    # objective function
    obj_type = MOI.get(src, MOI.ObjectiveFunctionType())
    obj = MOI.get(src, MOI.ObjectiveFunction{obj_type}())
    MOI.set(approx_model, MOI.ObjectiveFunction{obj_type}(), obj)

    obj_sense = MOI.get(src, MOI.ObjectiveSense())
    MOI.set(approx_model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # constraints
    con_sets = MOI.AbstractVectorSet[]
    con_funs = MOI.AbstractVectorFunction[]

    i = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        if S <: LinearSet
            # equality and orthant cone constraints
            MOIU.copyconstraints!(approx_model, src, false, idx_map, F, S)
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

    opt.approx_model = approx_model
    opt.con_sets = con_sets
    opt.con_funs = con_funs
    opt.status = :Loaded
    opt.idx_map = idx_map

    return idx_map
end

function MOI.optimize!(opt::Optimizer)
    solver = Solvers.SepCutsSolver(opt.approx_model, opt.con_sets, opt.con_funs) # TODO tols
    Solvers.solve(solver)

    opt.status = Solvers.get_status(solver)
    opt.solve_time = Solvers.get_solve_time(solver)
    opt.obj_value = Solvers.get_obj_value(solver)
    opt.obj_bound = Solvers.get_obj_bound(solver)

    return
end

# function MOI.free!(opt::Optimizer) # TODO ?

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

MOI.get(opt::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex) = MOI.get(opt.approx_model, MOI.VariablePrimal(), opt.idx_map[vi])

MOI.get(opt::Optimizer, a::MOI.VariablePrimal, vi::Vector{MOI.VariableIndex}) = MOI.get.(opt, a, vi)

MOI.get(opt::Optimizer, ::MOI.ConstraintDual, ci::MOI.ConstraintIndex{F, S}) where {F <: MOI.AbstractFunction, S <: LinearSet} =
    MOI.get(opt.approx_model, MOI.ConstraintDual(), opt.idx_map[ci]) # scalar constraint so in approx_model

function MOI.get(opt::Optimizer, ::MOI.ConstraintDual, ci::MOI.ConstraintIndex{F, S}) where {F <: MOI.AbstractFunction, S <: NonlinearSet}
    # dual vector corresponding to the constraint's polyhedral approximation is combination of the cuts weighted by their duals
    # TODO use the refs of the cuts for the constraint
    return []
end

MOI.get(opt::Optimizer, a::MOI.ConstraintDual, ci::Vector{MOI.ConstraintIndex}) = MOI.get.(opt, a, ci)

MOI.get(opt::Optimizer, ::MOI.ConstraintPrimal, ci::MOI.ConstraintIndex{F, S}) where {F <: MOI.AbstractFunction, S <: LinearSet} =
    MOI.get(opt.approx_model, MOI.ConstraintPrimal(), opt.idx_map[ci]) # scalar constraint so in approx_model

MOI.get(opt::Optimizer, ::MOI.ConstraintPrimal, ci::MOI.ConstraintIndex{F, S}) where {F <: MOI.AbstractFunction, S <: NonlinearSet} =
    MOIU.evalvariables(vi -> MOI.get(opt.approx_model, MOI.VariablePrimal(), vi), opt.con_funs[opt.idx_map[ci].value])

MOI.get(opt::Optimizer, a::MOI.ConstraintPrimal, ci::Vector{MOI.ConstraintIndex}) = MOI.get.(opt, a, ci)
