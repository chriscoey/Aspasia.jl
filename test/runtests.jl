#=
Copyright 2019, Chris Coey, Lea Kapelevich and contributors
=#

import Aspasia
const ASP = Aspasia

# import Hypatia
# const HYP = Hypatia
# const CO = HYP.Cones
# const MU = HYP.ModelUtilities

import Random
using LinearAlgebra
using SparseArrays
using Test

import MathOptInterface
const MOI = MathOptInterface

using Clp # approx solver


include(joinpath(@__DIR__, "MathOptInterface.jl"))

@testset "Aspasia tests" begin

@info("starting MathOptInterface tests")
verbose = false
approx_solver = Clp.Optimizer
@testset "MOI tests" begin
    test_moi(verbose, approx_solver)
end


end


# @testset "namedpoly dual" begin
#     include("../../Hypatia/examples/namedpoly/jump.jl")
#
#     (polyname, deg) = (:caprasse, 4)
#     (x, f, dom, truemin) = getpolydata(polyname)
#
#     # if use_wsos
#         model = build_JuMP_namedpoly_WSOS(x, f, dom, d = deg, primal_wsos = false)
#     # else
#     #     model = build_JuMP_namedpoly_PSD(x, f, dom, d = deg)
#     # end
#
#     JuMP.optimize!(model)
#
#     term_status = JuMP.termination_status(model)
#     primal_obj = JuMP.objective_value(model)
#     dual_obj = JuMP.objective_bound(model)
#     pr_status = JuMP.primal_status(model)
#     du_status = JuMP.dual_status(model)
#
#     @test term_status == MOI.OPTIMAL
#     @test pr_status == MOI.FEASIBLE_POINT
#     @test du_status == MOI.FEASIBLE_POINT
#     @test primal_obj ≈ dual_obj atol = 1e-4 rtol = 1e-4
#     @test primal_obj ≈ truemin atol = 1e-4 rtol = 1e-4
# end
