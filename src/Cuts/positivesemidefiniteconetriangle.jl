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

function get_sep_cuts(
    x::Vector{Float64},
    cone::MOI.PositiveSemidefiniteConeTriangle,
    tol::Float64,
)
    cuts = Vector{Float64}[]
    L = cone.side_dimension
    X = Matrix{Float64}(undef, L, L)
    vec_to_mat_noscale!(X, x)

    # TODO try only getting negative eigenvalues (what lower bound?) vs smallest k eigenvalues
    F = eigen!(Symmetric(X, :L), 1:min(L, 5))

    for k in eachindex(F.values)
        if F.values[k] >= -0.1 * tol
            continue
        end

        V = F.vectors[:, k]
        W = V * V'
        cut = similar(x)
        mat_to_vec_scale!(cut, W)

        if dot(x, cut) < -tol
            push!(cuts, cut)
        end
    end

    return cuts
end
