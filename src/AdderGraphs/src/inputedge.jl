"""



"""
function get_input_addernode(input_edge::InputEdge)
    return input_edge.input_adder
end

function get_input_shift(input_edge::InputEdge)
    return input_edge.shift
end

function is_negative_input(input_edge::InputEdge)
    return input_edge.is_negative
end


function get_input_addernode_value(input_edge::InputEdge)
    return get_value(get_input_addernode(input_edge))
end
