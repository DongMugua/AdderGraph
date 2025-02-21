module AdderGraphs

abstract type AbstractInputEdge end

# OriginAdder: AdderNode(1, Vector{InputEdge}(), 0, 0, -1.0, 0)
mutable struct AdderNode
    value::Int
    inputs::Vector{AbstractInputEdge}
    msb_out::Int
    lsb_out::Int
    error_out::Float64
    nb_full_adders::Int
end

mutable struct InputEdge <: AbstractInputEdge
    input_adder::AdderNode
    shift::Int
    is_negative::Bool
end

function origin_addernode()
    return AdderNode(1, Vector{InputEdge}(), 0, 0, -1.0, 0)
end

function AdderNode(value::Int, inputs::Vector{InputEdge})
    return AdderNode(value, inputs, 0, 0, -1.0, 0)
end


mutable struct AdderGraph
    origin::AdderNode
    constants::Vector{AdderNode}
    outputs::Vector{Int}
    full_adders_computed::Bool
    error_in::Float64
    errors_computed::Bool
    data_paths::Bool
end

function AdderGraph(c::Vector{AdderNode}, v::Vector{Int})
    return AdderGraph(origin_addernode(), c, v, false, 0.0, false, false)
end
function AdderGraph()
    return AdderGraph(Vector{AdderNode}(), Vector{Int}())
end
function AdderGraph(v::Vector{Int})
    return AdderGraph(Vector{AdderNode}(), v)
end
function AdderGraph(c::Vector{AdderNode})
    return AdderGraph(c, Vector{Int}())
end


export AdderGraph
export AdderNode
export InputEdge

include("utils.jl")
include("inputedge.jl")
include("addernode.jl")
include("addergraph.jl")

# Utils
export odd


# InputEdge
export get_input_addernode
export get_input_shift
export is_negative_input
export get_input_addernode_value


# AdderGraph
export length
export isempty
export get_nodes
export push_node!
export push_output!
export get_outputs
export get_error_in
export set_error_in!
export are_data_paths_set
export set_data_paths!
export are_errors_computed
export set_errors_computed!
export are_full_adders_computed
export set_full_adders_computed!
export get_origin
export get_lsb_in
export get_msb_in
export set_lsb_in!
export set_msb_in!
export get_total_nb_full_adders
export get_addernode_by_value
export read_addergraph
export write_addergraph


# AdderNode
export get_value
export get_input_edges
export get_adder_msb
export get_adder_lsb
export get_adder_error
export set_adder_msb!
export set_adder_lsb!
export set_adder_error!
export get_input_addernodes
export get_input_addernode_values
export get_input_shifts
export are_negative_inputs
export get_adder_msb_in
export get_adder_lsb_in
export get_nb_full_adders
export set_nb_full_adders!
export produce_addernode
export get_depth



end # module
