using JuMP, Ipopt
using Polyhedra, CDDLib, Plots
import GLPK;
solver = GLPK.Optimizer;
using SparseArrays, LinearAlgebra, Random
rng = MersenneTwister(1425)

include("algorithmes.jl")
include("directions.jl")
include("vectices.jl")
include("examples.jl")
include("plotting_checking.jl")

function test_example_1(dim::Int; explicit=false, ϵ=0.05)
    Xbar, SCCounter, V, Vunused, P, hr, Vused = algorith1(
        example1(; dim, explicit)...,
        chooseDirection=fixed_direction,
        chooseVertex=adj_direction,
        ϵ=ϵ,
        newDirection=adj_direction
    )
    check_example_1(Xbar)
    p = plot_xbar(Xbar)
    vscodedisplay(Xbar)
    p
end

function test_example_3(a::Int; ϵ=0.05)
    Xbar, SCCounter, V, Vunused, P, hr, Vused = algorith1(
        example3(a=a)...,
        chooseDirection=fixed_direction,
        chooseVertex=random_vertex,
        ϵ=ϵ,
        # newDirection=adj_direction
    )
    check_example_3(Xbar, a)
    p = plot_xbar(Xbar)
    vscodedisplay(Xbar)
    p
end



