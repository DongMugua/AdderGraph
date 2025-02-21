# Getters
function get_value(addernode::AdderNode)
    return addernode.value
end


"""



"""
function get_input_edges(addernode::AdderNode)
    return addernode.inputs
end


"""
    get_adder_msb(addernode::AdderNode)

Return the msb out of the adder `addernode`.
"""
function get_adder_msb(addernode::AdderNode)
    return addernode.msb_out
end


"""
    get_adder_lsb(addernode::AdderNode)

Return the lsb out of the adder `addernode`.
"""
function get_adder_lsb(addernode::AdderNode)
    return addernode.lsb_out
end


"""
    get_adder_error(addernode::AdderNode)

Return the error out of the adder `addernode` or `-1.0` if not specified yet.
"""
function get_adder_error(addernode::AdderNode)
    return addernode.error_out
end





# Setters
"""
    set_adder_msb(addernode::AdderNode, msb::Int)

Set the msb out of the adder `addernode` to `msb`.
"""
function set_adder_msb!(addernode::AdderNode, msb::Int)
    addernode.msb_out = msb
    return addernode
end


"""
    set_adder_lsb(addernode::AdderNode, lsb::Int)

Set the lsb out of the adder `addernode` to `lsb`.
"""
function set_adder_lsb!(addernode::AdderNode, lsb::Int)
    addernode.lsb_out = lsb
    return addernode
end


function set_adder_error!(addernode::AdderNode, error_out::Float64)
    addernode.error_out = error_out
    return addernode
end





# Extract more information
"""



"""
function get_input_addernodes(addernode::AdderNode)
    return get_input_addernode.(get_input_edges(addernode))
end


function get_input_addernode_values(addernode::AdderNode)
    return get_value.(get_input_addernode.(get_input_edges(addernode)))
end


function get_input_shifts(addernode::AdderNode)
    return get_input_shift.(get_input_edges(addernode))
end


function are_negative_inputs(addernode::AdderNode)
    return is_negative_input.(get_input_edges(addernode))
end


"""
    get_adder_msb_in(addernode::AdderNode)

Return the msb in of the adder `addernode`.
"""
function get_adder_msb_in(addernode::AdderNode)
    return get_input_shifts(addernode) .+ get_adder_msb.(get_input_addernodes(addernode))
end


"""
    get_adder_lsb_in(addernode::AdderNode)

Return the lsb in of the adder `addernode`.
"""
function get_adder_lsb_in(addernode::AdderNode)
    return get_input_shifts(addernode) .+ get_adder_lsb.(get_input_addernodes(addernode))
end



function get_nb_full_adders(addernode::AdderNode)
    return addernode.nb_full_adders
end


function set_nb_full_adders!(addernode::AdderNode, nb_full_adders::Int)
    addernode.nb_full_adders = nb_full_adders
    return addernode
end


# Generation
function produce_addernode(input_adders::Vector{AdderNode}, shifts::Vector{Int}, sign_switch::Vector{Bool})
    @assert length(input_values) == length(shifts)
    @assert length(shifts) == length(sign_switch)
    inputs = [InputEdge(input_values[i], shifts[i], sign_switch[i]) for i in 1:length(sign_switch)]
    addernode = AdderNode(_compute_value_from_inputs(inputs), inputs)
    return addernode
end



function get_depth(addernode::AdderNode)
    if get_value(addernode) == 1
        return 0
    end
    return maximum(get_depth.(get_input_addernodes(addernode)).+1)
end
