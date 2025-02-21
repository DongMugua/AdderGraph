module jIIR2HW

using JuMP
using AdderGraphs
using jIIR2AG
using jMCM

# Includes
include("vhdl.jl")
include("design.jl")
include("wcpg.jl")

export design2vhdl
export get_flopoco_cmd
export wcpg
export wcpg_eps


end # module
