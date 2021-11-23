"""
A Julia package for generic conic outer approximation methods.
"""
module Aspasia

using DocStringExtensions
using Printf

include("Cuts/Cuts.jl")

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const VI = MOI.VariableIndex
const SAF = MOI.ScalarAffineFunction{Float64}
const VV = MOI.VectorOfVariables
const VAF = MOI.VectorAffineFunction{Float64}

include("optimize.jl")
include("MOI_wrapper.jl")

end
