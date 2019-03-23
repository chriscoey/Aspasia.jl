#=
Copyright 2019, Chris Coey and contributors

functions and caches for cones
=#

module Cones

using LinearAlgebra
import JuMP

abstract type Cone end

# include("epinorminf.jl")
# include("epinormeucl.jl")
# include("epipersquare.jl")
# include("hypoperlog.jl")
# include("epiperpower.jl")
# include("epipersumexp.jl")
# include("hypogeomean.jl")
# include("epinormspectral.jl")
include("semidefinite.jl")
# include("hypoperlogdet.jl")
include("wsospolyinterp.jl")
include("wsospolyinterpmat.jl")

use_dual(cone::Cone) = cone.use_dual
dimension(cone::Cone) = cone.dim

# for PSD-like cones, use Λ operator with eigendecomposition to get cuts
PSDLikeCone = Union{PosSemidef, WSOSPolyInterp, WSOSPolyInterpMat}
function check_feas_get_cuts(cone::PSDLikeCone)
    tmp = JuMP.value(cone.Λ)
    F = eigen!(tmp, 1:5) # TODO play with only getting negative eigenvalues (what lower bound?) vs smallest k eigenvalues
    println("infeas is $(F.values[1])")
    return [F.vectors[:, k]' * cone.Λ * F.vectors[:, k] for k in eachindex(F.values) if F.values[k] < 1e-5]
end


end
