import Base.isempty

AdderID = Int
EdgeID = Int

mutable struct Edge
    id::EdgeID
    inputs::Tuple{AdderID, AdderID}
    output::AdderID
    weight::Int
    input_shift::Int
    output_shift::Int
    subtraction::Bool
end


mutable struct Adder
    id::AdderID
    value::Int
    depth::Int
end


mutable struct CompleteAdderGraph
    adders::Dict{AdderID, Adder}
    adders_edges_relation::Dict{AdderID, Set{EdgeID}}
    depth_adders_relation::Dict{Int, Set{AdderID}}
    edges::Dict{EdgeID, Edge}
    adders_number::Int
    edges_number::Int
    next_adder_id::Int
    next_edge_id::Int
    wordlength::Int
    maximum_value::Int
end


################################################################################


function add_adder!(AG::CompleteAdderGraph, value::Int, depth::Int)
    AG.adders_number += 1
    adder = Adder(AG.next_adder_id, value, depth)
    AG.adders[AG.next_adder_id] = adder
    push!(get!(AG.depth_adders_relation, depth, Set{AdderID}()), AG.next_adder_id)
    AG.next_adder_id += 1
    return AG.next_adder_id-1
end


function add_edge!(AG::CompleteAdderGraph,
                   inputs::Tuple{AdderID, AdderID},
                   output::AdderID,
                   weight::Int,
                   input_shift::Int,
                   output_shift::Int,
                   subtraction::Bool)
    AG.edges_number += 1
    #println(get_adder_value(AG, output))
    edge = Edge(AG.edges_number, inputs, output, weight, input_shift, output_shift, subtraction)
    for input_adder in inputs
        push!(get!(AG.adders_edges_relation, input_adder, Set{EdgeID}()), AG.next_edge_id)
    end
    push!(get!(AG.adders_edges_relation, output, Set{EdgeID}()), AG.next_edge_id)
    AG.edges[AG.next_edge_id] = edge
    AG.next_edge_id += 1
    return AG.next_edge_id-1
end


function remove_edge!(AG::CompleteAdderGraph, edge_id::EdgeID)
    for (adder_id, relation_set) in AG.adders_edges_relation
        delete!(relation_set, edge_id)
    end
    delete!(AG.edges, edge_id)
    AG.edges_number -= 1
    return AG
end


function remove_adder!(AG::CompleteAdderGraph, adder_id::AdderID)
    for edge_id in AG.adders_edges_relation[adder_id]
        delete!(AG.edges, edge_id)
        AG.edges_number -= 1
    end
    delete!(AG.adders_edges_relation, adder_id)
    delete!(AG.depth_adders_relation, adder_id)
    delete!(AG.adders, adder_id)
    AG.adders_number -= 1
    return AG
end


function get_adder_value(AG::CompleteAdderGraph, adder_id::AdderID)
    return AG.adders[adder_id].value
end


function get_adder_depth(AG::CompleteAdderGraph, adder_id::AdderID)
    return AG.adders[adder_id].depth
end


#TODO
function remove_isolated_adders!(AG::CompleteAdderGraph)
    return AG
end


function get_edges(AG::CompleteAdderGraph, adder_id::AdderID)
    return AG.adders_edges_relation[adder_id]
end


function get_weight(AG::CompleteAdderGraph, edge_id::EdgeID)
    return AG.edges[edge_id].weight
end


function get_input_adders(AG::CompleteAdderGraph, edge_id::EdgeID)
    return AG.edges[edge_id].inputs
end


function get_output_adder(AG::CompleteAdderGraph, edge_id::EdgeID)
    return AG.edges[edge_id].output
end


function get_previous_adders(AG::CompleteAdderGraph, adder_id::AdderID)
    adders = Set{AdderID}()
    for edge_id in get_edges(AG, adder_id)
        inputs = get_input_adders(AG, edge_id)
        append!(adders, Set(inputs))
    end
    return adders
end


function get_next_adders(AG::CompleteAdderGraph, adder_id::AdderID)
    adders = Set{AdderID}()
    for edge_id in get_edges(AG, adder_id)
        output = get_output_adder(AG, edge_id)
        append!(adders, Set(output))
    end
    return adders
end


function isempty(AG::CompleteAdderGraph)
    if AG.adders_number == 0
        return true
    end
    return false
end


function get_max_value(AG::CompleteAdderGraph)
    if isempty(AG)
        return nothing
    end
    return maximum(x -> get_adder_value(AG, x), keys(AG.adders))
end


function get_adder_id_at_depth(AG::CompleteAdderGraph, depth::Int)
    return AG.depth_adders_relation[depth]
end


function get_values_at_depth(AG::CompleteAdderGraph, depth::Int)
    return Set{Int}([get_adder_value(AG, adder_id) for adder_id in AG.depth_adders_relation[depth]])
end


function get_values_at_depth_dict(AG::CompleteAdderGraph, depth::Int)
    return Dict{Int, AdderID}([get_adder_value(AG, adder_id) => adder_id for adder_id in AG.depth_adders_relation[depth]])
end


function get_number_of_adders(AG::CompleteAdderGraph)
    return AG.adders_number
end


################################################################################


# A(u,u) = A(1,1)*u
function compute_successors!(AG::CompleteAdderGraph, adder_id::AdderID)
    #depth_adders_relation
    depth = get_adder_depth(AG, adder_id)
    next_depth = depth + 1
    adder_value = get_adder_value(AG, adder_id)
    dict_value_id = get_values_at_depth_dict(AG, next_depth)
    new_adder_id = get(dict_value_id, adder_value, nothing)
    if new_adder_id === nothing
        new_adder_id = add_adder!(AG, adder_value, next_depth)
    end
    add_edge!(AG, (adder_id, adder_id), new_adder_id, 0, 0, 1, false)
    c_max = AG.maximum_value
    tmp_max_j = floor(log2(c_max/adder_value - 1))
    if tmp_max_j >= 0
        max_j = Int(tmp_max_j)
        for j in 1:max_j
            value = adder_value * ((1 << j) - 1)
            new_adder_id = get(dict_value_id, value, nothing)
            if new_adder_id === nothing
                new_adder_id = add_adder!(AG, value, next_depth)
            end
            add_edge!(AG, (adder_id, adder_id), new_adder_id, 1, j, 0, true)
            value = adder_value * ((1 << j) + 1)
            new_adder_id = get(dict_value_id, value, nothing)
            if new_adder_id === nothing
                new_adder_id = add_adder!(AG, value, next_depth)
            end
            add_edge!(AG, (adder_id, adder_id), new_adder_id, 1, j, 0, false)
        end
        value = (1 << (max_j+1) - 1) * adder_value
        if value <= c_max
            new_adder_id = get(dict_value_id, value, nothing)
            if new_adder_id === nothing
                new_adder_id = add_adder!(AG, value, next_depth)
            end
            add_edge!(AG, (adder_id, adder_id), new_adder_id, 1, max_j+1, 0, true)
        end
    end

    return AG
end


function compute_successors!(AG::CompleteAdderGraph, adder_id_1::AdderID, adder_id_2::AdderID)
    #depth_adders_relation
    depth = get_adder_depth(AG, adder_id_1)
    @assert depth == get_adder_depth(AG, adder_id_2)
    next_depth = depth + 1

    dict_value_id = get_values_at_depth_dict(AG, next_depth)
    adder_value_1 = get_adder_value(AG, adder_id_1)
    adder_value_2 = get_adder_value(AG, adder_id_2)
    c_max = AG.maximum_value

    value = div(adder_value_1 + adder_value_2, 2)
    output_shift = 1
    while mod(value, 2) == 0
        value = div(value, 2)
        output_shift += 1
    end
    new_adder_id = get(dict_value_id, value, nothing)
    if new_adder_id === nothing
        new_adder_id = add_adder!(AG, value, next_depth)
    end
    add_edge!(AG, (adder_id_1, adder_id_2), new_adder_id, 1, 0, output_shift, false)

    value = div(abs(adder_value_1 - adder_value_2), 2)
    output_shift = 1
    while mod(value, 2) == 0
        value = div(value, 2)
        output_shift += 1
    end
    new_adder_id = get(dict_value_id, value, nothing)
    if new_adder_id === nothing
        new_adder_id = add_adder!(AG, value, next_depth)
    end
    add_edge!(AG, (adder_id_1, adder_id_2), new_adder_id, 1, 0, output_shift, true)

    max_j = floor(log2((c_max-adder_value_1)/adder_value_2))
    for j in 1:max_j
        value = adder_value_1 + (adder_value_2 << Int(j))
        new_adder_id = get(dict_value_id, value, nothing)
        if new_adder_id === nothing
            new_adder_id = add_adder!(AG, value, next_depth)
        end
        add_edge!(AG, (adder_id_1, adder_id_2), new_adder_id, 1, Int(j), 0, false)
    end
    max_j = floor(log2((c_max-adder_value_2)/adder_value_1))
    for j in 1:max_j
        value = adder_value_2 + (adder_value_1 << Int(j))
        new_adder_id = get(dict_value_id, value, nothing)
        if new_adder_id === nothing
            new_adder_id = add_adder!(AG, value, next_depth)
        end
        add_edge!(AG, (adder_id_2, adder_id_1), new_adder_id, 1, Int(j), 0, false)
    end

    max_j = floor(log2((c_max+adder_value_1)/adder_value_2))
    for j in 1:max_j
        value = abs(adder_value_1 - (adder_value_2 << Int(j)))
        new_adder_id = get(dict_value_id, value, nothing)
        if new_adder_id === nothing
            new_adder_id = add_adder!(AG, value, next_depth)
        end
        add_edge!(AG, (adder_id_1, adder_id_2), new_adder_id, 1, Int(j), 0, true)
    end
    max_j = floor(log2((c_max+adder_value_2)/adder_value_1))
    for j in 1:max_j
        value = abs(adder_value_2 - (adder_value_1 << Int(j)))
        new_adder_id = get(dict_value_id, value, nothing)
        if new_adder_id === nothing
            new_adder_id = add_adder!(AG, value, next_depth)
        end
        add_edge!(AG, (adder_id_2, adder_id_1), new_adder_id, 1, Int(j), 0, true)
    end

    return AG
end


function construct_addergraph(wordlength::Int, max_depth::Int)
    AG = CompleteAdderGraph(Dict{AdderID, Adder}(),
                            Dict{AdderID, Set{EdgeID}}(),
                            Dict{Int, Set{AdderID}}(),
                            Dict{EdgeID, Edge}(),
                            0, 0, 1, 1, wordlength, (1 << wordlength)-1)
    #
    adder_id_root = add_adder!(AG, 1, 0)
    # First layer
    current_layer = 1
    current_depth = 1
    adder_id = add_adder!(AG, 1, current_depth)
    add_edge!(AG, (adder_id_root, adder_id_root), adder_id, 0, 0, -1, false)
    for w in 2:(wordlength-1)
        adder_id = add_adder!(AG, 2^w+1, current_depth)
        add_edge!(AG, (adder_id_root, adder_id_root), adder_id, 1, w, 0, true)
    end
    for w in 2:wordlength
        adder_id = add_adder!(AG, 2^w-1, current_depth)
        add_edge!(AG, (adder_id_root, adder_id_root), adder_id, 1, w, 0, false)
    end

    # Values in first layer are always done in first layer:
    # No point of doing those values another way
    values_at_first_depth = sort(collect(get_values_at_depth(AG, 1)))


    # Second layer
    if max_depth >= 2
        previous_depth = current_depth
        current_depth += 1
        values_at_previous_depth = get_values_at_depth_dict(AG, previous_depth)
        # Start by copying the previous layer
        for adder_value in keys(values_at_previous_depth)
            adder_id = values_at_previous_depth[adder_value]
            new_adder_id = add_adder!(AG, adder_value, current_depth)
            add_edge!(AG, (adder_id, adder_id), new_adder_id, 0, 0, -1, false)
        end
        # Add new adders and edges
        # skip (1,1), in general skip what was done in previous-previous layer
        for (adder_value_1, adder_id_1) in values_at_previous_depth
            for (adder_value_2, adder_id_2) in values_at_previous_depth
                if adder_value_1 == 1 && adder_value_2 == 1
                    continue
                end
                if adder_value_1 == adder_value_2
                    compute_successors!(AG, adder_id_1)
                else
                    compute_successors!(AG, adder_id_1, adder_id_2)
                end
            end
        end
        # Remove adders which value can be done from (1,1)
        values_at_current_depth = get_values_at_depth_dict(AG, current_depth)
        for (adder_value, adder_id) in values_at_current_depth
            if adder_value in values_at_first_depth
                edges = copy(get_edges(AG, adder_id))
                for edge_id in edges
                    if get_weight(AG, edge_id) != 0
                        remove_edge!(AG, edge_id)
                    end
                end
            end
        end

        # Third layer
        if max_depth >= 3
            previous_depth = current_depth
            current_depth += 1
            values_at_previous_previous_depth = collect(keys(values_at_previous_depth))
            values_at_previous_depth = get_values_at_depth_dict(AG, previous_depth)
            # Start by copying the previous layer
            for adder_value in keys(values_at_previous_depth)
                adder_id = values_at_previous_depth[adder_value]
                new_adder_id = add_adder!(AG, adder_value, current_depth)
                add_edge!(AG, (adder_id, adder_id), new_adder_id, 0, 0, -1, false)
            end
            # Add new adders and edges
            # skip what was done in previous-previous layer
            for (adder_value_1, adder_id_1) in values_at_previous_depth
                for (adder_value_2, adder_id_2) in values_at_previous_depth
                    if adder_value_1 in values_at_previous_previous_depth && adder_value_2 in values_at_previous_previous_depth
                        continue
                    end
                    if adder_value_1 == adder_value_2
                        compute_successors!(AG, adder_id_1)
                    else
                        compute_successors!(AG, adder_id_1, adder_id_2)
                    end
                end
            end
            # Remove adders which value can be done from (1,1)
            values_at_current_depth = get_values_at_depth_dict(AG, current_depth)
            for (adder_value, adder_id) in values_at_current_depth
                if adder_value in values_at_first_depth
                    edges = copy(get_edges(AG, adder_id))
                    for edge_id in edges
                        if get_weight(AG, edge_id) != 0
                            remove_edge!(AG, edge_id)
                        end
                    end
                end
            end

            # More layers?

        end
    end

    # Return the CompleteAdderGraph
    #@info AG
    return AG
end
