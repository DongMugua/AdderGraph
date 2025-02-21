function get_addergraph_a(iir2addergraphs::Iir2AdderGraphs)
    return Iir2AdderGraphs.addergraph_a
end

function get_addergraph_b(iir2addergraphs::Iir2AdderGraphs)
    return Iir2AdderGraphs.addergraph_b
end

function get_coefficients_a(iir2addergraphs::Iir2AdderGraphs)
    return Iir2AdderGraphs.coefficients_a
end

function get_coefficients_b(iir2addergraphs::Iir2AdderGraphs)
    return Iir2AdderGraphs.coefficients_b
end

function get_shift_a(iir2addergraphs::Iir2AdderGraphs)
    return Iir2AdderGraphs.shifts[1]
end

function get_shift_b(iir2addergraphs::Iir2AdderGraphs)
    return Iir2AdderGraphs.shifts[2]
end
