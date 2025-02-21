module jIIR2AG

using JuMP
using Polynomials
using AdderGraphs

mutable struct Iir2AdderGraphs
    coefficients_a::Tuple{Int, Int}
    coefficients_b::Tuple{Int, Int, Int}
    shifts::Tuple{Int, Int}
    addergraph_a::AdderGraph
    addergraph_b::AdderGraph
end

function Iir2AdderGraphs()
    return Iir2AdderGraphs((0,0), (0,0,0), (0,0), AdderGraph(Vector{AdderNode}()), AdderGraph(Vector{AdderNode}()))
end

function Iir2AdderGraphs(
        coefficients_a::Tuple{Int, Int},
        coefficients_b::Tuple{Int, Int, Int},
        shifts::Tuple{Int, Int}
    )
    return Iir2AdderGraphs(coefficients_a, coefficients_b, shifts, AdderGraph(Vector{AdderNode}()), AdderGraph(Vector{AdderNode}()))
end

# Includes
include("iir2addergraphs.jl")
include("utils.jl")
include("addergraph.jl")
include("iir_ilp1.jl")
include("iir_ilp2.jl")
include("iir.jl")
include("mcm_ilp1.jl")
include("mcm_ilp2.jl")
#include("plots.jl")
include("presolve.jl")
include("design.jl")

export Iir2AdderGraphs
export design_second_order_iir
export all_coefficients
export get_specifications
export isempty
export get_addergraph_a
export get_addergraph_b
export get_coefficients_a
export get_coefficients_b
export get_shift_a
export get_shift_b

end # module
