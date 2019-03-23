#=
Copyright 2019, Chris Coey and contributors

definitions of conic sets not already defined by MathOptInterface
and functions for converting between Aspasia and MOI cone definitions
=#

export WSOSPolyInterpCone

struct WSOSPolyInterpCone <: MOI.AbstractVectorSet
    dimension::Int
    ipwt::Vector{Matrix{Float64}}
    is_dual::Bool
end
WSOSPolyInterpCone(dimension::Int, ipwt::Vector{Matrix{Float64}}) = WSOSPolyInterpCone(dimension, ipwt, false)

export WSOSPolyInterpMatCone

struct WSOSPolyInterpMatCone <: MOI.AbstractVectorSet
    R::Int
    U::Int
    ipwt::Vector{Matrix{Float64}}
    is_dual::Bool
end
WSOSPolyInterpMatCone(R::Int, U::Int, ipwt::Vector{Matrix{Float64}}) = WSOSPolyInterpMatCone(R, U, ipwt, false)

MOIOtherCones = (
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

# MOI cones for which no transformation is needed
# cone_from_moi(s::MOI.SecondOrderCone) = Cones.EpiNormEucl(MOI.dimension(s))
# cone_from_moi(s::MOI.RotatedSecondOrderCone) = Cones.EpiPerSquare(MOI.dimension(s))
# cone_from_moi(s::MOI.ExponentialCone) = Cones.HypoPerLog()
# cone_from_moi(s::MOI.GeometricMeanCone) = (l = MOI.dimension(s) - 1; Cones.HypoGeomean(fill(inv(l), l)))
# cone_from_moi(s::MOI.PowerCone{Float64}) = Cones.EpiPerPower(inv(s.exponent))
cone_from_moi(s::WSOSPolyInterpCone) = Cones.WSOSPolyInterp(s.dimension, s.ipwt, s.is_dual)
cone_from_moi(s::WSOSPolyInterpMatCone) = Cones.WSOSPolyInterpMat(s.R, s.U, s.ipwt, s.is_dual)
cone_from_moi(s::MOI.AbstractVectorSet) = error("MOI set $s is not recognized")

function build_var_cone(fi::MOI.VectorOfVariables, si::MOI.AbstractVectorSet, dim::Int, q::Int)
    IGi = (q + 1):(q + dim)
    VGi = -ones(dim)
    conei = cone_from_moi(si)
    return (IGi, VGi, conei)
end

function build_constr_cone(fi::MOI.VectorAffineFunction{Float64}, si::MOI.AbstractVectorSet, dim::Int, q::Int)
    IGi = [q + vt.output_index for vt in fi.terms]
    VGi = [-vt.scalar_term.coefficient for vt in fi.terms]
    Ihi = (q + 1):(q + dim)
    Vhi = fi.constants
    conei = cone_from_moi(si)
    return (IGi, VGi, Ihi, Vhi, conei)
end

# MOI cones requiring transformations (eg rescaling, changing order)
# TODO later remove if MOI gets scaled triangle sets
const rt2 = sqrt(2)
const rt2i = inv(rt2)

svec_scale(dim) = [(i == j ? 1.0 : rt2) for i in 1:round(Int, sqrt(0.25 + 2.0 * dim) - 0.5) for j in 1:i]
svec_unscale(dim) = [(i == j ? 1.0 : rt2i) for i in 1:round(Int, sqrt(0.25 + 2.0 * dim) - 0.5) for j in 1:i]

# PSD cone: convert from smat to svec form (scale off-diagonals)
function build_var_cone(fi::MOI.VectorOfVariables, si::MOI.PositiveSemidefiniteConeTriangle, dim::Int, q::Int)
    IGi = (q + 1):(q + dim)
    VGi = -svec_scale(dim)
    conei = Cones.PosSemidef(dim)
    return (IGi, VGi, conei)
end

function build_constr_cone(fi::MOI.VectorAffineFunction{Float64}, si::MOI.PositiveSemidefiniteConeTriangle, dim::Int, q::Int)
    scalevec = svec_scale(dim)
    IGi = [q + vt.output_index for vt in fi.terms]
    VGi = [-vt.scalar_term.coefficient * scalevec[vt.output_index] for vt in fi.terms]
    Ihi = (q + 1):(q + dim)
    Vhi = scalevec .* fi.constants
    conei = Cones.PosSemidef(dim)
    return (IGi, VGi, Ihi, Vhi, conei)
end

# # logdet cone: convert from smat to svec form (scale off-diagonals)
# function build_var_cone(fi::MOI.VectorOfVariables, si::MOI.LogDetConeTriangle, dim::Int, q::Int)
#     IGi = (q + 1):(q + dim)
#     VGi = vcat(-1.0, -1.0, -svec_scale(dim - 2))
#     conei = Cones.HypoPerLogdet(dim)
#     return (IGi, VGi, conei)
# end
#
# function build_constr_cone(fi::MOI.VectorAffineFunction{Float64}, si::MOI.LogDetConeTriangle, dim::Int, q::Int)
#     scalevec = vcat(1.0, 1.0, svec_scale(dim - 2))
#     IGi = [q + vt.output_index for vt in fi.terms]
#     VGi = [-vt.scalar_term.coefficient * scalevec[vt.output_index] for vt in fi.terms]
#     Ihi = (q + 1):(q + dim)
#     Vhi = scalevec .* fi.constants
#     conei = Cones.HypoPerLogdet(dim)
#     return (IGi, VGi, Ihi, Vhi, conei)
# end
