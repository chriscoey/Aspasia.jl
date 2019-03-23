#=
Copyright 2019, Chris Coey and contributors

row-wise lower triangle (smat space) of positive semidefinite matrix cone
(smat space) W \in S^n : 0 >= eigmin(W)
(see equivalent MathOptInterface PositiveSemidefiniteConeTriangle definition)
=#

mutable struct PosSemidef <: Cone
    use_dual::Bool
    side::Int
    Λ::Symmetric{Matrix{JuMP.AffExpr}}
    # tmp::Symmetric{Matrix{Float64}}

    function PosSemidef(x::Vector{JuMP.AffExpr}, is_dual::Bool = false)
        cone = new()
        cone.use_dual = is_dual
        cone.side = round(Int, sqrt(0.25 + 2.0 * length(x)) - 0.5)
        cone.Λ = construct_Λ(cone, x)
        # cone.tmp = Matrix{Float64}(undef, cone.side, cone.side)
        return cone
    end
end

# construct Λ, a symmetric matrix of JuMP affine expressions (only once when building cone object)
function construct_Λ(cone::PosSemidef, x::Vector{JuMP.AffExpr})
    Λ = Symmetric(Matrix{JuMP.AffExpr}(undef, cone.side, cone.side), :L)
    k = 1
    for i in 1:cone.side, j in 1:i
        Λ[i, j] = x[k]
        k += 1
    end
    return Λ
end
