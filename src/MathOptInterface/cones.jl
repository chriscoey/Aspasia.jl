#=
Copyright 2019, Chris Coey and contributors

definitions of conic sets not already defined by MathOptInterface
and oracles for polyhedral approximation
=#

const rt2 = sqrt(2)
const rt2i = inv(rt2)

# lower triangle mat to vec
function mat_to_vec!(vec::Vector{Float64}, mat::AbstractMatrix{Float64})
    k = 1
    m = size(mat, 1)
    for i in 1:m, j in 1:i # fill from lower triangle
        vec[k] = mat[i, j]
        k += 1
    end
    return vec
end

# vec to lower triangle mat
function vec_to_mat!(mat::AbstractMatrix{Float64}, vec::Vector{Float64})
    k = 1
    m = size(mat, 1)
    for i in 1:m, j in 1:i # fill lower triangle
        mat[i, j] = vec[k]
        k += 1
    end
    return mat
end

#=
PSD cone
cuts are <V*V^T, X> for eigenvectors V corresponding to negative eigenvalues of matrix point X
=#
function get_cuts(x::Vector{Float64}, cone::MOI.PositiveSemidefiniteConeTriangle)
    L = cone.side_dimension
    X = Matrix{Float64}(undef, L, L)
    vec_to_mat!(X, x)

    F = eigen!(Symmetric(X, :L), 1:min(L, 5)) # TODO try only getting negative eigenvalues (what lower bound?) vs smallest k eigenvalues

    cuts = Vector{Float64}[]
    for k in eachindex(F.values)
        if F.values[k] >= 1e-5
            continue
        end
        V = F.vectors[:, k]
        W = V * V'
        cut = similar(x)
        mat_to_vec!(cut, W) # TODO syrk
        push!(cuts, cut)

        @show dot(x, cut)
    end

    return cuts
end

#=
interpolation-based weighted-sum-of-squares (multivariate) polynomial cone parametrized by interpolation points ipwt
cuts are <V*V^T, Λ> for eigenvectors V corresponding to negative eigenvalues of Λ
Λ is given by # TODO
=#
export WSOSPolyInterpCone

struct WSOSPolyInterpCone <: MOI.AbstractVectorSet
    dimension::Int
    ipwt::Vector{Matrix{Float64}}
    is_dual::Bool
end
WSOSPolyInterpCone(dimension::Int, ipwt::Vector{Matrix{Float64}}) = WSOSPolyInterpCone(dimension, ipwt, false)
dimension(cone::WSOSPolyInterpCone) = cone.dimension

function get_cuts(x::Vector{Float64}, cone::WSOSPolyInterpCone)
    cuts = Vector{Float64}[]
    U = cone.dimension
    for w in eachindex(cone.ipwt)
        P = cone.ipwt[w]
        L = size(P, 2)
        X = Symmetric(P' * Diagonal(x) * P)
        # X = Matrix{Float64}(undef, L, L)
        # for i in 1:L, j in 1:i
        #     X[i, j] = sum(Wt[u, i] * Wt[u, j] * x[u] for u in eachindex(x))
        # end

        F = eigen!(X, 1:min(L, 5)) # TODO try only getting negative eigenvalues (what lower bound?) vs smallest k eigenvalues

        for k in eachindex(F.values)
            if F.values[k] >= 1e-5
                continue
            end
            V = F.vectors[:, k]
            cut = [] # TODO V * P' * Diagonal(vars) * P * V' # TODO syrk
            push!(cuts, cut)
        end
    end
    return cuts
end

#=
interpolation-based weighted-sum-of-squares (multivariate) polynomial matrix cone parametrized by interpolation points ipwt
cuts are <V*V^T, Λ> for eigenvectors V corresponding to negative eigenvalues of Λ
Λ is given by # TODO
=#
export WSOSPolyInterpMatCone

struct WSOSPolyInterpMatCone <: MOI.AbstractVectorSet
    R::Int
    U::Int
    ipwt::Vector{Matrix{Float64}}
    is_dual::Bool
end
WSOSPolyInterpMatCone(R::Int, U::Int, ipwt::Vector{Matrix{Float64}}) = WSOSPolyInterpMatCone(R, U, ipwt, false)
dimension(cone::WSOSPolyInterpMatCone) = cone.dimension

function get_cuts(x::Vector{Float64}, cone::WSOSPolyInterpMatCone)
    cuts = Vector{Float64}[]
    R = cone.R
    U = cone.U
    for w in eachindex(cone.ipwt)
        P = cone.ipwt[w]
        L = size(P, 2)
        RL = R * L
        X = Matrix{Float64}(undef, RL, RL)
        # TODO use matrix multiplications below instead of scalar, maybe
        u0 = 0
        i0 = 0
        for p in 1:R
            j0 = 0
            for q in 1:(p - 1)
                # off diagonal block, so rescale by rt2i TODO maybe don't define cone using rescaling
                for i in 1:RL, j in 1:RL
                    X[i0 + i, j0 + j] = X[i0 + j, j0 + i] = sum(rt2i * P[u, i] * P[u, j] * x[u0 + u] for u in eachindex(x))
                end
                u0 += U
                j0 += L
            end
            # on diagonal block, so only do lower triangle
            for i in 1:RL, j in 1:i
                X[i0 + i, j0 + j] = sum(P[u, i] * P[u, j] * x[u0 + u] for u in eachindex(x))
            end
            i0 += L
        end

        F = eigen!(Symmetric(X, :L), 1:min(RL, 5)) # TODO try only getting negative eigenvalues (what lower bound?) vs smallest k eigenvalues

        for k in eachindex(F.values)
            if F.values[k] >= 1e-5
                continue
            end
            V = F.vectors[:, k]
            cut = [] # TODO V * ... * V' # TODO syrk
            push!(cuts, cut)
        end
    end
    return cuts
end
