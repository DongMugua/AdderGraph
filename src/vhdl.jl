"""
    run_flopoco(flopoco_cmd::Cmd)

Run FloPoCo via command line and silence output
"""
function run_flopoco(flopoco_cmd::Cmd)
    run(pipeline(flopoco_cmd, stdout=devnull, stderr=devnull))
    return nothing
end


"""
    get_flopoco_cmd(...)

Generate FloPoCo command
"""
function get_flopoco_cmd(filename::String, coefficients_a::Vector{Float64}, coefficients_b::Vector{Float64};
                         lsbin::Int, lsbout::Int, msbin::Int=-1)
    return `flopoco FixIIR coeffa="$(join(coefficients_a, ":"))" coeffb="$(join(coefficients_b, ":"))" msbIn=$msbin lsbIn=$(lsbin) lsbOut=$(lsbin) outputFile="$(filename)"`
end

function get_flopoco_cmd(filename::String, coefficients_a::Vector{String}, coefficients_b::Vector{String}; lsbin::Int, lsbout::Int, msbin::Int=-1)
    coefficients_a .*= "0x"
    coefficients_b .*= "0x"
    return `flopoco FixIIR coeffa="$(join(coefficients_a, ":"))" coeffb="$(join(coefficients_b, ":"))" msbIn=$msbin lsbIn=$(lsbin) lsbOut=$(lsbin) outputFile="$(filename)"`
end


function get_flopoco_cmd(filename::String, coefficients_a::Vector{Int}, coefficients_b::Vector{Int}, shifta::Int, shiftb::Int; lsbin::Int, lsbout::Int, msbin::Int=-1, guardbits::Int=-1, wcpg_val::Float64=-1, wcpg_eps_val::Float64=-1)
    return `flopoco FixIIRShiftAdd coeffa="$(join(coefficients_a, ":"))" coeffb="$(join(coefficients_b, ":"))" shifta=$(shifta) shiftb=$(shiftb) msbIn=$msbin lsbIn=$(lsbin) lsbOut=$(lsbin) guardBits=$(guardbits) method="plain" $(wcpg_val!=-1 ? "H=$(wcpg_val)" : "") $(wcpg_eps_val!=-1 ? "Heps=$(wcpg_eps_val)" : "") outputFile="$(filename)"`
end

function get_flopoco_cmd(filename::String, addergraph_a::AdderGraph, addergraph_b::AdderGraph, shifta::Int, shiftb::Int; lsbin::Int, lsbout::Int, msbin::Int=-1, guardbits::Int=-1, wcpg_val::Float64=-1, wcpg_eps_val::Float64=-1)
    coefficients_a = get_outputs(addergraph_a)
    coefficients_b = get_outputs(addergraph_b)
    addergraph_a_str = write_addergraph(addergraph_a)
    addergraph_b_str = write_addergraph(addergraph_b)
    return `flopoco FixIIRShiftAdd coeffa="$(join(coefficients_a, ":"))" coeffb="$(join(coefficients_b, ":"))" shifta=$(shifta) shiftb=$(shiftb) grapha="$(addergraph_a_str)" graphb="$(addergraph_b_str)" msbIn=$(msbin) lsbIn=$(lsbin) lsbOut=$(lsbin) guardBits=$(guardbits) method="multiplierless" $(wcpg_val!=-1 ? "H=$(wcpg_val)" : "") $(wcpg_eps_val!=-1 ? "Heps=$(wcpg_eps_val)" : "") outputFile="$(filename)"`
end

"""
    write_vhdl(...)

Generate and run a FloPoCo command that produces a VHDL file
"""
function write_vhdl(filename::String, solution::Iir2AdderGraphs, args...; kwargs...)
    run_flopoco(get_flopoco_cmd(filename, solution.addergraph_a, solution.addergraph_b, solution.shifts[1], solution.shifts[2], args...; kwargs...))
    return nothing
end

function write_vhdl(args...; kwargs...)
    run_flopoco(get_flopoco_cmd(args...; kwargs...))
    return nothing
end
