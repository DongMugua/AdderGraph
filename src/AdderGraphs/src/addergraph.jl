import Base.length
import Base.isempty

function length(addergraph::AdderGraph)
    return length(addergraph.constants)
end

function isempty(addergraph::AdderGraph)
    return (length(addergraph.constants)+length(addergraph.outputs)==0)
end

function get_nodes(addergraph::AdderGraph)
    return addergraph.constants
end

function push_node!(addergraph::AdderGraph, addernode::AdderNode)
    push!(addergraph.constants, addernode)
    return addergraph
end

function push_node!(addergraph::AdderGraph, inputs::Vector{InputEdge})
    push_node!(addergraph, AdderNode(_compute_value_from_inputs(inputs), inputs))
end

function push_output!(addergraph::AdderGraph, output_value::Int)
    push!(addergraph.outputs, output_value)
    return addergraph
end

function get_outputs(addergraph::AdderGraph)
    return addergraph.outputs
end


function get_error_in(addergraph::AdderGraph)
    return get_adder_error(addergraph.origin)
end

function set_error_in!(addergraph::AdderGraph, error_in::Float64)
    set_adder_error!(addergraph.origin, error_in)
    return addergraph
end

function are_data_paths_set(addergraph::AdderGraph)
    return addergraph.data_paths
end

function set_data_paths!(addergraph::AdderGraph, data_paths_value::Bool=true)
    addergraph.data_paths=data_paths_value
    return addergraph
end

function are_errors_computed(addergraph::AdderGraph)
    return addergraph.errors_computed
end

function set_errors_computed!(addergraph::AdderGraph, errors_computed_value::Bool=true)
    addergraph.errors_computed = errors_computed_value
    return addergraph
end

function are_full_adders_computed(addergraph::AdderGraph)
    return addergraph.full_adders_computed
end

function set_full_adders_computed!(addergraph::AdderGraph, full_adders_computed_value::Bool=true)
    addergraph.full_adders_computed=full_adders_computed_value
    return addergraph
end

function get_origin(addergraph::AdderGraph)
    return addergraph.origin
end

function get_lsb_in(addergraph::AdderGraph)
    return get_adder_lsb(get_origin(addergraph))
end
function get_msb_in(addergraph::AdderGraph)
    return get_adder_msb(get_origin(addergraph))
end

function set_lsb_in!(addergraph::AdderGraph, lsb::Int)
    set_adder_lsb!(get_origin(addergraph), lsb)
    return addergraph
end
function set_msb_in!(addergraph::AdderGraph, msb::Int)
    set_adder_msb!(get_origin(addergraph), msb)
    return addergraph
end

function get_total_nb_full_adders(addergraph::AdderGraph)
    return sum(get_nb_full_adders(addernode) for addernode in get_nodes(addergraph))
end

function get_addernode_by_value(addergraph::AdderGraph, value::Int)
    if get_value(get_origin(addergraph)) == value
        return get_origin(addergraph)
    end
    for addernode in get_nodes(addergraph)
        if get_value(addernode) == value
            return addernode
        end
    end
    return nothing
end


function read_addergraph(s::String)
    addergraph = AdderGraph()
    adders_registers_outputs = split(s[3:(end-2)], "},{")
    for val in adders_registers_outputs
        if startswith(val, "'A'")
            node_details = split(val, ",")
            node_value = parse(Int, node_details[2][2:(end-1)])
            node_inputs_signed = Vector{Int}()
            node_input_shifts = Vector{Int}()
            for current_input in 1:div(length(node_details)-3, 3)
                push!(node_inputs_signed, parse(Int, node_details[3*current_input+1][2:(end-1)]))
                push!(node_input_shifts, parse(Int, node_details[3*current_input+3]))
            end
            node_inputs = Vector{Int}(abs.(node_inputs_signed))
            node_subtraction = Vector{Bool}(abs.(div.(sign.(node_inputs_signed).-1, 2)))
            push_node!(addergraph, AdderNode(node_value, [InputEdge(get_addernode_by_value(addergraph, node_inputs[i]), node_input_shifts[i], node_subtraction[i]) for i in 1:length(node_subtraction)]))
        elseif startswith(val, "'O'")
            push_output!(addergraph, parse(Int, split(val, ",")[2][2:(end-1)]))
        end
    end
    return addergraph
end


function write_addergraph(addergraph::AdderGraph)
    depth_by_value = Dict([1 => 0])
    maximum_depth = 0
    for addernode in get_nodes(addergraph)
        depth_by_value[get_value(addernode)] = maximum(depth_by_value[get_input_addernode_values(addernode)[i]] for i in 1:length(get_input_addernode_values(addernode)))+1
        if depth_by_value[get_value(addernode)] > maximum_depth
            maximum_depth = depth_by_value[get_value(addernode)]
        end
    end
    output_values = Vector{Int}()
    adderstring = "{"
    firstcoefnocomma = true
    coefficients = get_outputs(addergraph)
    for coefind in 1:length(coefficients)
        coef = coefficients[coefind]
        if coef != 0
            value = odd(abs(coef))
            if !(abs(coef) in output_values)
                push!(output_values, abs(coef))
                shift = round(Int, log2(abs(coef)/value))
                if !firstcoefnocomma
                    adderstring *= ","
                else
                    firstcoefnocomma = false
                end
                adderstring *= "{'O',[$(abs(coef))],$(maximum_depth+1),[$(value)],$(depth_by_value[value]),$(shift)}"
            end
        end
    end
    for addernode in get_nodes(addergraph)
        adderstring *= ",{'A',[$(get_value(addernode))],$(depth_by_value[get_value(addernode)])"
        for input_edge in get_input_edges(addernode)
            adderstring *= ",[$((-1)^(is_negative_input(input_edge))*get_input_addernode_value(input_edge))],"
            adderstring *= "$(depth_by_value[get_input_addernode_value(input_edge)]),$(get_input_shift(input_edge))"
        end
        adderstring *= "}"
    end
    adderstring *= "}"
    return adderstring
end
