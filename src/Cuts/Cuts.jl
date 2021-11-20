"""
Cut oracles for MathOptInterface cones.
"""
module Cuts

using DocStringExtensions
using LinearAlgebra
import LinearAlgebra.copytri!

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

include("secondordercone.jl")
include("positivesemidefiniteconetriangle.jl")

end
