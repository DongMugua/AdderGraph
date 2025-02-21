HypernodeID = Int
HyperarcID = Int
import Base.isempty

mutable struct Hyperarc
    id::HyperarcID
    input_nodes::Set{HypernodeID}
    output_nodes::Set{HypernodeID}
    weight::Int
end

mutable struct Hypernode
    id::HypernodeID
    value::Int
end

mutable struct Hypergraph
    hypernodes::Dict{HypernodeID, Hypernode}
    hypernodes_arc_relation::Dict{HypernodeID, Set{HyperarcID}}
    hyperarcs::Dict{HyperarcID, Hyperarc}
    nodes_number::Int
    arcs_number::Int
    node_id::Int
    arc_id::Int
end

function add_node!(h::Hypergraph, value::Int)
    h.nodes_number += 1
    h.node_id += 1
    node = Hypernode(h.node_id, value)
    h.hypernodes[h.node_id] = node
    return h.node_id
end


function add_hyperarc!(h::Hypergraph,
                       input_nodes::Set{HypernodeID},
                       output_nodes::Set{HypernodeID},
                       weight::Int)
    h.arcs_number += 1
    h.arc_id += 1
    hyperarc = Hyperarc(h.arcs_number, input_nodes, output_nodes, weight)
    for node_id in input_nodes
        push!(get!(h.hypernodes_arc_relation, node_id, Set{HyperarcID}()), h.arc_id)
    end
    for node_id in output_nodes
        push!(get!(h.hypernodes_arc_relation, node_id, Set{HyperarcID}()), h.arc_id)
    end
    h.hyperarcs[h.arc_id] = hyperarc
    return h.arc_id
end


function add_input_node_hyperarc!(h::Hypergraph, node_id::HypernodeID, arc_id::HyperarcID)
    push!(get!(h.hypernodes_arc_relation, node_id, Set{HyperarcID}()), arc_id)
    push!(h.hyperarcs[arc_id].input_nodes, node_id)
    return h
end

function add_output_node_hyperarc!(h::Hypergraph, node_id::HypernodeID, arc_id::HyperarcID)
    push!(get!(h.hypernodes_arc_relation, node_id, Set{HyperarcID}()), arc_id)
    push!(h.hyperarcs[arc_id].output_nodes, node_id)
    return h
end


function remove_hyperarc!(h::Hypergraph, arc_id::HyperarcID)
    for (hypernode_id, relation_set) in h.hypernodes_arc_relation
        delete!(relation_set, arc_id)
    end
    delete!(h.hyperarcs, arc_id)
    h.arcs_number -= 1
    return h
end

function remove_hypernode!(h::Hypergraph, node_id::HypernodeID)
    for arc_id in h.hypernodes_arc_relation[node_id]
        delete!(h.hyperarcs[arc_id].input_nodes, node_id)
        delete!(h.hyperarcs[arc_id].output_nodes, node_id)
        if isempty(h.hyperarcs[arc_id].input_nodes) || isempty(h.hyperarcs[arc_id].output_nodes)
            delete!(h.hyperarcs, arc_id)
            h.arcs_number -= 1
        end
    end
    delete!(h.hypernodes_arc_relation, node_id)
    delete!(h.hypernodes, node_id)
    h.nodes_number -= 1
    return h
end

function get_node_value(h::Hypergraph, node_id::HypernodeID)
    return h.hypernodes[node_id].value
end


#TODO
function remove_isolated_nodes!(h::Hypergraph)
    return h
end


function get_arcs(h::Hypergraph, node_id::HypernodeID)
    return h.hypernodes_arc_relation[node_id]
end

function get_input_nodes(h::Hypergraph, arc_id::HyperarcID)
    return h.hyperarcs[arc_id].input_nodes
end

function get_output_nodes(h::Hypergraph, arc_id::HyperarcID)
    return h.hyperarcs[arc_id].output_nodes
end


function get_previous_nodes(h::Hypergraph, node_id::HypernodeID)
    V = Set{HypernodeID}()
    for arc_id in get_arcs(h, node_id)
        input_nodes = get_input_nodes(h, arc_id)
        if !(node_id in input_nodes)
            append!(V, input_nodes)
        end
    end

    return V
end

function get_next_nodes(h::Hypergraph, node_id::HypernodeID)
    V = Set{HypernodeID}()
    for arc_id in get_arcs(h, node_id)
        output_nodes = get_output_nodes(h, arc_id)
        if !(node_id in output_nodes)
            append!(V, output_nodes)
        end
    end

    return V
end


function isempty(h::Hypergraph)
    if h.nodes_number == 0
        return true
    end
    return false
end


function get_max_value(h::Hypergraph)
    if isempty(h)
        return nothing
    end
    maximum_value = maximum(x -> get_node_value(h, x), keys(h.hypernodes))
end


function aset(u::Int, v::Int, wordlength::Int)
    if u > v
        u, v = v, u
    end
    X = Set{Int}()
    val_max = 1 << wordlength
    push!(X, odd(u+v))
    if u != v
        push!(X, odd(abs(u - v)))
        s = 1
        tmp_plus = v << s + u
        while tmp_plus < val_max
            push!(X, tmp_plus)
            push!(X, abs(v << s - u))
            push!(X, u << s + v)
            push!(X, abs(u << s - v))
            s += 1
            tmp_plus = v << s + u
        end
        s_old = s
        tmp_minus = abs(v << s - u)
        while tmp_minus < val_max
            push!(X, tmp_minus)
            s += 1
            tmp_minus = abs(v << s - u)
        end
        s = s_old
        tmp_plus = u << s + v
        while tmp_plus < val_max
            push!(X, tmp_plus)
            push!(X, abs(u << s - v))
            s += 1
            tmp_plus = u << s + v
        end
        tmp_minus = abs(u << s - v)
        while tmp_minus < val_max
            push!(X, tmp_minus)
            s += 1
            tmp_minus = abs(u << s - v)
        end
    else
        s = 1
        tmp_plus = v << s + u
        while tmp_plus < val_max
            push!(X, tmp_plus)
            push!(X, abs(v << s - u))
            s += 1
            tmp_plus = v << s + u
        end
        tmp_minus = abs(v << s - u)
        while tmp_minus < val_max
            push!(X, tmp_minus)
            s += 1
            tmp_minus = abs(v << s - u)
        end
    end

    return sort!(collect(X))
end


"""
    construct_hypertree(wordlength::Int, layers::Int = 3)

Return the hypergraph that will be used by ILP2.
"""
function construct_hypertree(wordlength::Int, layers::Int = 3)
    h = Hypergraph(Dict{HypernodeID, Hypernode}(),
                   Dict{HypernodeID, Set{HyperarcID}}(),
                   Dict{HyperarcID, Hyperarc}(),
                   0, 0, 0, 0)
    #
    node_id_root = add_node!(h, 1)
    nodes_id_layer = Vector{Vector{HypernodeID}}()
    push!(nodes_id_layer, Vector{HypernodeID}())
    # First layer
    current_layer = 1
    push!(nodes_id_layer[current_layer], add_node!(h, 1))
    for value in [2^w+minusplusone for w in 2:(wordlength-1) for minusplusone in [-1,1]]
        push!(nodes_id_layer[current_layer], add_node!(h, value))
    end
    push!(nodes_id_layer[current_layer], add_node!(h, 2^wordlength-1))
    add_hyperarc!(h, Set([node_id_root]), Set([nodes_id_layer[current_layer][1]]), 0)
    for node_id in nodes_id_layer[current_layer][2:end]
        add_hyperarc!(h, Set([node_id_root]), Set([node_id]), 1)
    end

    # Values in first layer are always done in first layer:
    # No point of doing those values another way
    values_of_first_layer = [get_node_value(h, node_id) for node_id in nodes_id_layer[1]]

    # Second layer
    if layers >= 2
        push!(nodes_id_layer, Vector{HypernodeID}())
        previous_layer = current_layer
        current_layer += 1
        values_of_current_layer = Dict{Int, HypernodeID}()
        # Start by recopying the previous layer
        for node_id in nodes_id_layer[previous_layer]
            value = get_node_value(h, node_id)
            new_node_id = add_node!(h, value)
            values_of_current_layer[value] = new_node_id
            push!(nodes_id_layer[current_layer], new_node_id)
            add_hyperarc!(h, Set([node_id]), Set([new_node_id]), 0)
        end
        # Add new nodes and arcs
        # skip (1,1), in general skip what was done in previous-previous layer
        for i in 1:1
            for j in 2:length(nodes_id_layer[previous_layer])
                id_i = nodes_id_layer[previous_layer][i]
                id_j = nodes_id_layer[previous_layer][j]
                X = aset(get_node_value(h, id_i), get_node_value(h, id_j), wordlength)
                for value in X
                    if !(value in values_of_first_layer)
                        node_id = get(values_of_current_layer, value, nothing)
                        if node_id === nothing
                            node_id = add_node!(h, value)
                            values_of_current_layer[value] = node_id
                            push!(nodes_id_layer[current_layer], node_id)
                        end
                        add_hyperarc!(h, Set([id_i, id_j]), Set([node_id]), 1)
                    end
                end
            end
        end
        for i in 2:length(nodes_id_layer[previous_layer])
            for j in i:length(nodes_id_layer[previous_layer])
                id_i = nodes_id_layer[previous_layer][i]
                id_j = nodes_id_layer[previous_layer][j]
                X = aset(get_node_value(h, id_i), get_node_value(h, id_j), wordlength)
                for value in X
                    if !(value in values_of_first_layer)
                        node_id = get(values_of_current_layer, value, nothing)
                        if node_id === nothing
                            node_id = add_node!(h, value)
                            values_of_current_layer[value] = node_id
                            push!(nodes_id_layer[current_layer], node_id)
                        end
                        add_hyperarc!(h, Set([id_i, id_j]), Set([node_id]), 1)
                    end
                end
            end
        end
        # Third layer
        if layers >= 3
            push!(nodes_id_layer, Vector{HypernodeID}())
            previous_layer = current_layer
            current_layer += 1
            values_of_current_layer = Dict{Int, HypernodeID}()
            # Start by recopying the previous layer
            for node_id in nodes_id_layer[previous_layer]
                value = get_node_value(h, node_id)
                new_node_id = add_node!(h, value)
                values_of_current_layer[value] = new_node_id
                push!(nodes_id_layer[current_layer], new_node_id)
                add_hyperarc!(h, Set([node_id]), Set([new_node_id]), 0)
            end
            # Add new nodes and arcs
            # skip what was done in previous-previous layer
            for i in 1:length(nodes_id_layer[previous_layer])
                for j in i:length(nodes_id_layer[previous_layer])
                    id_i = nodes_id_layer[previous_layer][i]
                    id_j = nodes_id_layer[previous_layer][j]
                    if !(id_i in nodes_id_layer[current_layer-2]) || !(id_j in nodes_id_layer[current_layer-2])
                        X = aset(get_node_value(h, id_i), get_node_value(h, id_j), wordlength)
                        for value in X
                            if !(value in values_of_first_layer)
                                node_id = get(values_of_current_layer, value, nothing)
                                if node_id === nothing
                                    node_id = add_node!(h, value)
                                    values_of_current_layer[value] = node_id
                                    push!(nodes_id_layer[current_layer], node_id)
                                end
                                add_hyperarc!(h, Set([id_i, id_j]), Set([node_id]), 1)
                            end
                        end
                    end
                end
            end
        end
    end

    # More layers?

    return h
end
