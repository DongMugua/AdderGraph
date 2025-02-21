using JuMP
using CPLEX
#using MATLAB
# include("../src/AdderGraphs/src/AdderGraphs.jl")
# using .AdderGraphs
# include("../src/jMCM/src/jMCM.jl")
# using .jMCM
# include("../src/jIIR2AG/src/jIIR2AG.jl")
# using .jIIR2AG
include("../src/includepkg.jl")


function transfer_function_hp0(val::Float64)
    # if mod(val, 2) > 2-0.00001 || mod(val, 2) < 0.00001
    #     return NaN
    # end
    #b0, b1, b2, a0, a1, a2 = 101.8, -203.3964, 101.596419342, 1.0, -1.968, 0.96800951
    #return abs((b0*exp(2*im*val*pi)+b1*exp(im*val*pi)+b2)/(a0*exp(2*im*val*pi)+a1*exp(im*val*pi)+a2))
    return abs(((exp(im*val*pi)-0.9981)*(exp(im*val*pi)-0.9999))/((exp(im*val*pi)-0.9683)*(exp(im*val*pi)-0.9997)))
end

function our_method()
    hp0 = ("hp0", 6, (0.0, 0.0), (0.041, 1.0), 0.0, 0.016)
    instancename, wordlength, stopband, passband, deltastop, deltapass = hp0
    fbands = [passband, stopband]
    abands = [(1-deltapass, 1+deltapass), (0.0, deltastop)]
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    solution = design_second_order_iir(model, fbands, abands, wordlength, avoid_internal_shifts=true, presolve_time_sec=600, verbose=true)

    return solution
end

function transfer_function_coefs(val, a1, a2, b0, b1, b2)
    return abs((exp(2*im*val*pi)*b0+exp(im*val*pi)*b1+b2)/(exp(2*im*val*pi)+exp(im*val*pi)*a1+a2))
end


function plot_solutions()
    x = 0:0.01:1
    y_origin = transfer_function_hp0.(x)
    str_hp0_init = ""
    for i in 1:length(x)
        str_hp0_init *= "($(x[i]), $(y_origin[i])) "
    end
    solution = our_method()
    a1, a2 = solution.coefficients_a
    b0, b1, b2 = solution.coefficients_b
    shifta, shiftb = solution.shifts
    y_ours = [transfer_function_coefs(val, a1/(2^shifta), a2/(2^shifta), b0/(2^shiftb), b1/(2^shiftb), b2/(2^shiftb)) for val in x]
    str_hp0_ours = ""
    for i in 1:length(x)
        str_hp0_ours *= "($(x[i]), $(y_ours[i])) "
    end

    str_plot = """
\\begin{figure}
\\centering
\\begin{tikzpicture}
	\\begin{axis}[
			ytick pos=left,
			xtick pos=bottom,
			legend style={at={(0.9,0.1)}, anchor=south east},
			ymin=0,
			xmin=0,
			xmax=1,
			width=\\linewidth,
			ylabel={Frequency response \$\\left|H\\left(e^{i\\omega}\\right)\\right|\$},
			xlabel={Normalized frequency \$\\omega\$},
		]
		\\addplot [
		legend entry=Ours,
		color=green,
		thick,
		] coordinates
		{$(str_hp0_ours)};
		\\addplot [
		legend entry=Original,
		color=blue,
		thick,
		%dashed,
		] coordinates
		{$(str_hp0_init)};
	\\end{axis}
\\end{tikzpicture}
\\caption{Frequency response of hp0}
\\label{fig:hp0}
\\end{figure}
"""
    println(str_plot)
    return nothing
end




# function our_method()
#     hp0 = ("hp0", 10, transfer_function_hp0, 1, 0.01)
#     instancename, wordlength, transfer_function, error_margin_percent, error_minimum = hp0
#     model = Model(CPLEX.Optimizer)
#     set_silent(model)
#     #set_optimizer_attributes(model, "CPX_PARAM_INTSOLLIM" => 1)
#     if wordlength >= 10
#         if wordlength >= 12
#             if wordlength > 16
#                 error("ILP1 not usable due to numerical instability")
#             else
#                 set_optimizer_attributes(model, "CPX_PARAM_EPINT" => 10^-8)
#             end
#         else
#             set_optimizer_attributes(model, "CPX_PARAM_EPINT" => 10^-7)
#         end
#     end
#     solutions = all_coefficients(model, transfer_function, error_margin_percent, error_minimum, wordlength, presolve_time_sec=3000, with_symmetry_breaking=true, verbose=true)
#
#     filename = (@__DIR__)*"/hp0.txt"
#     open(filename, "w") do f
#     write(f, "a1;a2;b0;b1;b2;shifta;shiftb\n")
#         for solution in solutions
#             write(f, "$(solution.coefficients_a[1]);$(solution.coefficients_a[2]);$(solution.coefficients_b[1]);$(solution.coefficients_b[2]);$(solution.coefficients_b[3]);$(solution.shifts[1]);$(solution.shifts[2])\n")
#         end
#     end
#
#     nbAddersMin = typemax(Int)
#     nbAddersDict = Dict{Set{Int}, Int}()
#     wordlength = 10
#     open(filename) do fread
#         lines = readlines(fread)
#         for line in lines[2:end]
#             a1_init, a2_init, b0_init, b1_init, b2_init, shifta_init, shiftb_init = parse.(Int, split(line, ";"))
#             a1, a2, b0, b1, b2, shifta, shiftb = odd.(abs.(parse.(Int, split(line, ";"))))
#             modela = Model(CPLEX.Optimizer)
#             modelb = Model(CPLEX.Optimizer)
#             set_silent(modela)
#             set_silent(modelb)
#             nbAdders = get!(nbAddersDict, Set{Int}([a1, a2]),
#                 length(mcm(modela, [a1, a2], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true))) +
#                        get!(nbAddersDict, Set{Int}([b0, b1, b2]),
#                 length(mcm(modelb, [b0, b1, b2], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true))) +
#                 sum(coef != 0 ? 1 : 0 for coef in [a1, a2, b0, b1, b2]) - 1
#             if nbAdders < nbAddersMin
#                 nbAddersMin = nbAdders
#                 println("$nbAdders --> $a1_init/(1<<$shifta_init), $a2_init/(1<<$shifta_init), $b0_init/(1<<$shiftb_init), $b1_init/(1<<$shiftb_init), $b2_init/(1<<$shiftb_init)")
#             end
#         end
#     end
#     #solution = design_second_order_iir(model, transfer_function, error_margin_percent, error_minimum, wordlength, presolve_time_sec=3000, avoid_internal_shifts=true, verbose=true, nb_adders_lb=5)
#     return solutions
# end
