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
