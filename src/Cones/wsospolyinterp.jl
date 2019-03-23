#=
Copyright 2019, Chris Coey, Lea Kapelevich and contributors

interpolation-based weighted-sum-of-squares (multivariate) polynomial cone parametrized by interpolation points ipwt

definition from "Sum-of-squares optimization without semidefinite programming" by D. Papp and S. Yildiz, available at https://arxiv.org/abs/1712.01792
=#

mutable struct WSOSPolyInterp <: Cone
    use_dual::Bool
    U::Int
    Ls::Vector{Int}
    ipwt::Vector{Matrix{Float64}}
    Λs::Vector{Symmetric{Matrix{JuMP.AffExpr}}}
    # tmp::Vector{Symmetric{Matrix{Float64}}}

    function WSOSPolyInterp(x::Vector{JuMP.AffExpr}, U::Int, ipwt::Vector{Matrix{Float64}}, is_dual::Bool = false)
        cone = new()
        cone.use_dual = !is_dual # using dual cone oracle
        cone.U = U
        cone.Ls = [size(ipwtj, 2) for ipwtj in ipwt]
        cone.ipwt = ipwt
        cone.Λs = [construct_Λ(cone, x, ipwt[j], cone.Ls[j]) for j in eachindex(ipwt)]
        # cone.tmp = [Matrix{Float64}(undef, L, L) for L in cone.Ls]
        return cone
    end
end

# construct Λ, a symmetric matrix of JuMP affine expressions (only once when building cone object)
function construct_Λ(cone::WSOSPolyInterp, x::Vector{JuMP.AffExpr}, ipwtj::Matrix{Float64}, L::Int)
    Λ = Symmetric(Matrix{JuMP.AffExpr}(undef, L, L), :L)
    for i in 1:L, j in 1:i
        Λ[i, j] = sum(ipwtj[u, i] * ipwtj[u, j] * x[u] for u in eachindex(x))
    end
    return Λ
end
