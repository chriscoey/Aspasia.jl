#=
Copyright 2019, Chris Coey, Lea Kapelevich and contributors

interpolation-based weighted-sum-of-squares (multivariate) polynomial matrix cone parametrized by interpolation points ipwt
=#

mutable struct WSOSPolyInterpMat <: Cone
    use_dual::Bool
    R::Int
    U::Int
    Ls::Vector{Int}
    ipwt::Vector{Matrix{Float64}}
    Λs::Vector{Symmetric{Matrix{JuMP.AffExpr}}}
    # tmp::Vector{Symmetric{Matrix{Float64}}}

    function WSOSPolyInterpMat(x::Vector{JuMP.AffExpr}, R::Int, U::Int, ipwt::Vector{Matrix{Float64}}, is_dual::Bool = false)
        cone = new()
        cone.use_dual = !is_dual # using dual cone oracle
        cone.R = R
        cone.U = U
        cone.Ls = [size(ipwtj, 2) for ipwtj in ipwt]
        cone.ipwt = ipwt
        cone.Λs = [construct_Λ(cone, x, ipwt[j], cone.Ls[j]) for j in eachindex(ipwt)]
        # cone.tmp = [Matrix{Float64}(undef, R * L, R * L) for L in cone.Ls]
        return cone
    end
end

# construct Λ, a symmetric matrix of JuMP affine expressions (only once when building cone object)
function construct_Λ(cone::WSOSPolyInterpMat, x::Vector{JuMP.AffExpr}, ipwtj::Matrix{Float64}, L::Int)
    side = cone.R * L
    Λ = Symmetric(Matrix{JuMP.AffExpr}(undef, side, side), :L)
    u0 = 0
    i0 = 0
    for p in 1:cone.R
        j0 = 0
        for q in 1:(p - 1)
            # off diagonal block, so rescale by rt2i TODO maybe don't define cone using rescaling
            for i in 1:side, j in 1:side
                Λ[i0 + i, j0 + j] = Λ[i0 + j, j0 + i] = sum(rt2i * ipwtj[u, i] * ipwtj[u, j] * x[u0 + u] for u in eachindex(x))
            end
            u0 += cone.U
            j0 += L
        end
        # on diagonal block, so only do lower triangle
        for i in 1:side, j in 1:i
            Λ[i0 + i, j0 + j] = sum(ipwtj[u, i] * ipwtj[u, j] * x[u0 + u] for u in eachindex(x))
        end
        i0 += L
    end
    return Λ
end
