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

include("wrapper.jl")
include("optimize.jl")

end
