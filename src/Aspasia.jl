#=
Copyright 2019, Chris Coey, Lea Kapelevich and contributors
=#

module Aspasia

# submodules
include("Cones/Cones.jl")
include("Solvers/Solvers.jl")

using LinearAlgebra
using SparseArrays
import MathOptInterface
const MOI = MathOptInterface
import JuMP

include("MathOptInterface/cones.jl")
include("MathOptInterface/wrapper.jl")

end
