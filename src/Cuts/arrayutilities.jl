
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

# # get (i,j) index in symmetric matrix from lower triangle index
# vec_to_mat_idx(k::Int) = div(1 + isqrt(8 * k - 7), 2)

# get lower triangle vector index from (i,j) index in symmetric matrix
mat_to_vec_idx(i::Int, j::Int) = (i < j) ? (div((j - 1) * j, 2) + i) : (div((i - 1) * i, 2) + j)