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

    approx_solver::MOI.AbstractOptimizer

    approx_model::MOI.ModelLike
    constraint_idxs::Vector{MOI.ConstraintIndex}

    status::Symbol
    solve_time::Float64
    primal_obj::Float64
    dual_obj::Float64

    function Optimizer(approx_solver::MOI.AbstractOptimizer, verbose::Bool, max_iters::Int, time_limit::Float64, tol_rel_opt::Float64, tol_abs_opt::Float64, tol_feas::Float64)
        opt = new()
        opt.approx_solver
        opt.verbose = verbose
        opt.max_iters = max_iters
        opt.time_limit = time_limit
        opt.tol_rel_opt = tol_rel_opt
        opt.tol_abs_opt = tol_abs_opt
        opt.tol_feas = tol_feas
        opt.status = :NotLoaded
        return opt
    end
end

Optimizer(
    approx_solver::MOI.AbstractOptimizer;
    verbose::Bool = false,
    max_iters::Int = 500,
    time_limit::Float64 = 3.6e3, # TODO should be Inf
    tol_rel_opt::Float64 = 1e-6,
    tol_abs_opt::Float64 = 1e-7,
    tol_feas::Float64 = 1e-7,
    ) = Optimizer(approx_solver, verbose, max_iters, time_limit, tol_rel_opt, tol_abs_opt, tol_feas)

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

ApproxCones = ( # cones to use polyhedral approximation on
    # MOI.SecondOrderCone,
    # MOI.RotatedSecondOrderCone,
    # MOI.ExponentialCone,
    # MOI.PowerCone{Float64},
    # MOI.GeometricMeanCone,
    MOI.PositiveSemidefiniteConeTriangle,
    # MOI.LogDetConeTriangle,
    WSOSPolyInterpCone,
    WSOSPolyInterpMatCone,
    )

# TODO don't restrict to Float64 type
# TODO allow SOC/quadratic constraints if approx solver allows specifying?
SupportedFuns = Union{
    MOI.SingleVariable, MOI.ScalarAffineFunction{Float64},
    MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64},
    }

SupportedSets = Union{
    MOI.EqualTo{Float64}, MOI.Zeros,
    MOI.GreaterThan{Float64}, MOI.Nonnegatives,
    MOI.LessThan{Float64}, MOI.Nonpositives,
    MOI.Interval{Float64},
    ApproxCones...
    }

MOI.supports_constraint(::Optimizer, ::Type{<:SupportedFuns}, ::Type{<:SupportedSets}) = true

linear_scalar_funs = (MOI.SingleVariable, MOI.ScalarAffineFunction{Float64})
linear_scalar_sets = (MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Interval{Float64})
linear_vector_funs = (MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64})
linear_vector_sets = (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives)

# build polyhedral approximation MOI model and list of cone objects
function MOI.copy_to(
    opt::Optimizer,
    src::MOI.ModelLike;
    copy_names::Bool = false,
    warn_attributes::Bool = true,
    )
    @assert !copy_names
    idx_map = Dict{MOI.Index, MOI.Index}()

    approx_model = opt.approx_solver()

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
    @assert obj_sense == MOI.MIN_SENSE # TODO generalize
    MOI.set(approx_model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # constraints
    get_src_cons(F, S) = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    get_con_fun(con_idx) = MOI.get(src, MOI.ConstraintFunction(), con_idx)
    get_con_set(con_idx) = MOI.get(src, MOI.ConstraintSet(), con_idx)

    # equality and orthant cone constraints
    function add_linear_constraints(F, S)
        for ci in get_src_cons(F, S)
            old_fun = get_con_fun(ci)
            new_fun = # TODO .... # TODO need the function to use new variable indices
            idx_map[ci] = MOI.add_constraint(approx_model, new_fun, get_con_set(ci))
        end
    end
    for F in linear_scalar_funs, S in linear_scalar_sets # scalar constraints
        add_linear_constraints(F, S)
    end
    for F in linear_vector_funs, S in linear_vector_sets # vector constraints
        add_linear_constraints(F, S)
    end

    # nonpolyhedral cone constraints
    constraint_idxs = MOI.ConstraintIndex[]
    # TODO decide whether to add new variable vector equal to VAF, so cuts are only on variables
    for F in linear_vector_funs, S in ApproxCones, ci in get_src_cons(F, S)
        push!(constraint_idxs)
        idx_map[ci] = MOI.ConstraintIndex{F, S}(length(constraint_idxs))
    end

    opt.approx_model = approx_model
    opt.constraint_idxs = constraint_idxs
    opt.status = :Loaded

    return idx_map
end

function MOI.optimize!(opt::Optimizer)
    solver = Solvers.SepCutsSolver(opt.approx_model, opt.cones) # TODO tols
    Solvers.solve(solver)

    opt.status = Solvers.get_status(solver)
    opt.solve_time = Solvers.get_solve_time(solver)
    opt.primal_obj = Solvers.get_primal_obj(solver)
    opt.dual_obj = Solvers.get_dual_obj(solver)

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

function MOI.get(opt::Optimizer, ::MOI.ObjectiveValue)
    if opt.obj_sense == MOI.MIN_SENSE
        return opt.primal_obj
    # elseif opt.obj_sense == MOI.MAX_SENSE
    #     return -opt.primal_obj
    else
        error("no objective sense is set")
    end
end

function MOI.get(opt::Optimizer, ::MOI.ObjectiveBound)
    if opt.obj_sense == MOI.MIN_SENSE
        return opt.dual_obj
    # elseif opt.obj_sense == MOI.MAX_SENSE
    #     return -opt.dual_obj
    else
        error("no objective sense is set")
    end
end

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





# TODO same variables and same linear constraints, but for nonpolyhedral constraint need to combine cut information

MOI.get(opt::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex) = opt.x[vi.value]
MOI.get(opt::Optimizer, a::MOI.VariablePrimal, vi::Vector{MOI.VariableIndex}) = MOI.get.(opt, a, vi)

function MOI.get(opt::Optimizer, ::MOI.ConstraintDual, ci::MOI.ConstraintIndex{F, S}) where {F <: MOI.AbstractFunction, S <: MOI.AbstractScalarSet}
    # scalar set
    i = ci.value
    if i <= opt.num_eq_constrs
        # constraint is an equality
        return opt.y[opt.constr_offset_eq[i] + 1]
    else
        # constraint is conic
        i -= opt.num_eq_constrs
        return opt.z[opt.constr_offset_cone[i] + 1]
    end
end
function MOI.get(opt::Optimizer, ::MOI.ConstraintDual, ci::MOI.ConstraintIndex{F, S}) where {F <: MOI.AbstractFunction, S <: MOI.AbstractVectorSet}
    # vector set
    i = ci.value
    if i <= opt.num_eq_constrs
        # constraint is an equality
        os = opt.constr_offset_eq
        return opt.y[(os[i] + 1):os[i + 1]]
    else
        # constraint is conic
        i -= opt.num_eq_constrs
        os = opt.constr_offset_cone
        return opt.z[(os[i] + 1):os[i + 1]]
    end
end
MOI.get(opt::Optimizer, a::MOI.ConstraintDual, ci::Vector{MOI.ConstraintIndex}) = MOI.get.(opt, a, ci)

function MOI.get(opt::Optimizer, ::MOI.ConstraintPrimal, ci::MOI.ConstraintIndex{F, S}) where {F <: MOI.AbstractFunction, S <: MOI.AbstractScalarSet}
    # scalar set
    i = ci.value
    if i <= opt.num_eq_constrs
        # constraint is an equality
        return opt.constr_prim_eq[opt.constr_offset_eq[i] + 1]
    else
        # constraint is conic
        i -= opt.num_eq_constrs
        return opt.constr_prim_cone[opt.constr_offset_cone[i] + 1]
    end
end
function MOI.get(opt::Optimizer, ::MOI.ConstraintPrimal, ci::MOI.ConstraintIndex{F, S}) where {F <: MOI.AbstractFunction, S <: MOI.AbstractVectorSet}
    # vector set
    i = ci.value
    if i <= opt.num_eq_constrs
        # constraint is an equality
        os = opt.constr_offset_eq
        return opt.constr_prim_eq[(os[i] + 1):os[i + 1]]
    else
        # constraint is conic
        i -= opt.num_eq_constrs
        os = opt.constr_offset_cone
        return opt.constr_prim_cone[(os[i] + 1):os[i + 1]]
    end
end
MOI.get(opt::Optimizer, a::MOI.ConstraintPrimal, ci::Vector{MOI.ConstraintIndex}) = MOI.get.(opt, a, ci)
