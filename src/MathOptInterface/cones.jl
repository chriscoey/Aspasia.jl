#=
Copyright 2019, Chris Coey and contributors

definitions of conic sets not already defined by MathOptInterface
and oracles for polyhedral approximation

TODO refactor eigendecomps and cuts on lambda
=#

const rt2 = sqrt(2)
const rt2i = inv(rt2)

# lower triangle mat to vec
function mat_to_vec_scale!(vec::Vector{Float64}, mat::AbstractMatrix{Float64})
    k = 1
    m = size(mat, 1)
    for i in 1:m # fill vector from lower triangle, adding off-diagonals twice
        for j in 1:(i - 1)
            vec[k] = 2.0 * mat[i, j]
            k += 1
        end
        vec[k] = mat[i, i]
        k += 1
    end
    return vec
end

# vec to lower triangle mat
function vec_to_mat_noscale!(mat::AbstractMatrix{Float64}, vec::Vector{Float64})
    k = 1
    m = size(mat, 1)
    for i in 1:m, j in 1:i # fill lower triangle from vector
        mat[i, j] = vec[k]
        k += 1
    end
    return mat
end

# get (i,j) index in symmetric matrix from lower triangle index
vec_to_mat_idx(k::Int) = div(1 + isqrt(8 * k - 7), 2)

# get lower triangle vector index from (i,j) index in symmetric matrix
mat_to_vec_idx(i::Int, j::Int) = (i < j) ? (div((j - 1) * j, 2) + i) : (div((i - 1) * i, 2) + j)


#=
second-order cone
cuts are (1, -y/||y||_2) for (z, y) outside the cone
TODO extended formulation
=#
function get_init_cuts(cone::MOI.SecondOrderCone)
    cuts = Vector{Float64}[]
    dim = MOI.dimension(cone)

    # z is nonnegative
    cut = zeros(dim)
    cut[1] = 1.0
    push!(cuts, cut)

    # z is at least absolute value of each y_i
    for i in 2:dim
        cut = zeros(dim)
        cut[1] = 1.0
        cut[i] = 1.0
        push!(cuts, cut)
        cut = copy(cut)
        cut[i] = -1.0
        push!(cuts, cut)
    end

    return cuts
end

function get_sep_cuts(x::Vector{Float64}, cone::MOI.SecondOrderCone)
    dim = MOI.dimension(cone)
    y_norm = norm(x[i] for i in 2:dim)

    if y_norm - x[1] > 1e-7
        cut = similar(x)
        cut[1] = 1
        for i in 2:dim
            cut[i] = -x[i] / y_norm
        end

        if dot(x, cut) < -1e-6 # TODO tolerance option
            return [cut]
        end
    end

    return Vector{Float64}[]
end

#=
PSD cone
cuts are V*V^T for eigenvectors V corresponding to negative eigenvalues of matrix point X
=#
function get_init_cuts(cone::MOI.PositiveSemidefiniteConeTriangle)
    cuts = Vector{Float64}[]
    L = cone.side_dimension
    dim = MOI.dimension(cone)

    # diagonal variables are nonnegative
    for i in 1:L
        cut = zeros(dim)
        cut[mat_to_vec_idx(i, i)] = 1.0
        push!(cuts, cut)
    end

    # 2x2 principal minors satisfy linearizations of PSD/RSOC condition
    for i in 1:L, j in 1:(i - 1)
        cut = zeros(dim)
        cut[mat_to_vec_idx(i, i)] = 1.0
        cut[mat_to_vec_idx(j, j)] = 1.0
        ij_idx = mat_to_vec_idx(i, j)
        cut[ij_idx] = 2.0
        push!(cuts, cut)
        cut = copy(cut)
        cut[ij_idx] = -2.0
        push!(cuts, cut)
    end

    return cuts
end

function get_sep_cuts(x::Vector{Float64}, cone::MOI.PositiveSemidefiniteConeTriangle)
    cuts = Vector{Float64}[]
    L = cone.side_dimension
    X = Matrix{Float64}(undef, L, L)
    vec_to_mat_noscale!(X, x)

    F = eigen!(Symmetric(X, :L), 1:min(L, 5)) # TODO try only getting negative eigenvalues (what lower bound?) vs smallest k eigenvalues

    for k in eachindex(F.values)
        if F.values[k] >= -1e-7 # TODO tolerance option
            continue
        end

        V = F.vectors[:, k]
        W = V * V' # TODO syrk
        cut = similar(x)
        mat_to_vec_scale!(cut, W)

        if dot(x, cut) < -1e-6 # TODO tolerance option
            push!(cuts, cut)
        end
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

function get_init_cuts(cone::WSOSPolyInterpCone)
    cuts = Vector{Float64}[]
    U = cone.dimension

    for w in eachindex(cone.ipwt)
        P = cone.ipwt[w]
        L = size(P, 2)

        PiPj = Matrix{Vector{Float64}}(undef, L, L)
        for i in 1:L, j in 1:i
            PiPj[i, j] = P[:, i] .* P[:, j]
        end

        # diagonal of Λ is nonnegative
        for i in 1:L
            push!(cuts, PiPj[i, i])
        end

        # 2x2 principal minors of Λ satisfy linearizations of PSD/RSOC condition
        for i in 1:L, j in 1:(i - 1)
            cut = PiPj[i, i] + PiPj[j, j] + 2.0 * PiPj[i, j]
            push!(cuts, cut)
            cut = PiPj[i, i] + PiPj[j, j] - 2.0 * PiPj[i, j]
            push!(cuts, cut)
        end
    end

    return cuts
end

# TODO maybe faster if store the PiPj matrix calculated for initial cuts and re-use here
function get_sep_cuts(x::Vector{Float64}, cone::WSOSPolyInterpCone)
    cuts = Vector{Float64}[]
    U = cone.dimension

    for w in eachindex(cone.ipwt)
        P = cone.ipwt[w]
        L = size(P, 2)
        X = Symmetric(P' * Diagonal(x) * P)
        # X = Matrix{Float64}(undef, L, L)
        # for i in 1:L, j in 1:i
        #     X[i, j] = sum(P[u, i] * P[u, j] * x[u] for u in eachindex(x))
        # end

        F = eigen!(X, 1:min(L, 5)) # TODO try only getting negative eigenvalues (what lower bound?) vs smallest k eigenvalues
        # F = eigen!(X, -Inf, -1e-7)
        # F = eigen!(X)

        @show F.values

        for k in eachindex(F.values)
            if F.values[k] >= -1e-5
                continue
            end

            V = F.vectors[:, k]
            PV = P * V
            cut = abs2.(PV) # V' * P' * Diagonal(vars) * P * V

            @show dot(x, cut)

            if dot(x, cut) < -1e-5 # TODO tolerance option
                push!(cuts, cut)
            end
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

function get_sep_cuts(x::Vector{Float64}, cone::WSOSPolyInterpMatCone)
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
            if F.values[k] >= 1e-6
                continue
            end

            V = F.vectors[:, k]
            cut = [] # TODO V * ... * V' # TODO syrk
            push!(cuts, cut)
        end
    end

    return cuts
end
