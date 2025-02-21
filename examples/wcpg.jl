using JuMP
using CPLEX
using Polynomials
# include("../src/AdderGraphs/src/AdderGraphs.jl")
# using .AdderGraphs
# include("../src/jMCM/src/jMCM.jl")
# using .jMCM
# include("../src/jIIR2AG/src/jIIR2AG.jl")
# using .jIIR2AG
include("../src/includepkg.jl")

using Plots

function WCPG(a::Vector{Float64}, b::Vector{Float64})::Float64
    W = zeros(1)
    if ccall((:WCPG_tf, "libwcpg"), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cuintmax_t, Cuintmax_t, Cint), W, b, a, length(b), length(a), 0) == 0
        @warn "WCPG computation errored"
    end
    return W[1]
end


function WCPG_eps(a::Vector{Float64}, b::Vector{Float64})::Float64
    W = zeros(1)
    if ccall((:WCPG_tf, "libwcpg"), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cuintmax_t, Cuintmax_t, Cint), W, [1.0], a, 1, length(a), 0) == 0
        @warn "WCPG_eps computation errored"
    end
    return W[1]
end


function get_poles_coordinates(a::Vector{Float64})
    poles_complex = roots(Polynomial(push!(reverse(a),1)))

    return [(real(pole_complex), imag(pole_complex)) for pole_complex in poles_complex]
end


function plot_all_wcpg(;wcpgeps::Bool=false)
    wcpg_fun = WCPG
    if wcpgeps
        println("WCPG of the error filter")
        wcpg_fun = WCPG_eps
    end
    named_instances = [
        #("compilation", [4], (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        #("hp1", 4:10, (0.7, 1.0), (0.0, 0.3), 0.1, 0.1),
        ("hp0", 4:7, (0.041, 1.0), (0.0, 0.0), 0.016, 0.0),
        ("lp1x4", 4:6, (0.0, 0.3), (0.7, 1.0), 0.06, 0.06),
        ("lp1x0", 4:6, (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        ("lp1x1", 4:6, (0.0, 0.3), (0.7, 1.0), 0.09, 0.09),
        ("lp1x2", 4:6, (0.0, 0.3), (0.7, 1.0), 0.08, 0.08),
        ("lp1x3", 4:6, (0.0, 0.3), (0.7, 1.0), 0.07, 0.07),
        ("lp1x5", 4:6, (0.0, 0.3), (0.7, 1.0), 0.05, 0.05),
        ("lp1x6", 4:6, (0.0, 0.3), (0.7, 1.0), 0.04, 0.04),
        ("lp1x7", 4:6, (0.0, 0.3), (0.7, 1.0), 0.03, 0.03),
        ("lp1x8", 4:6, (0.0, 0.3), (0.7, 1.0), 0.02, 0.02),
        ("lp1x9", 4:6, (0.0, 0.3), (0.7, 1.0), 0.01, 0.01),
        ("lp2x0", 4:6, (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        ("lp2x1", 4:6, (0.0, 0.3), (0.65, 1.0), 0.1, 0.1),
        ("lp2x2", 4:6, (0.0, 0.3), (0.6, 1.0), 0.1, 0.1),
        ("lp2x3", 4:6, (0.0, 0.3), (0.55, 1.0), 0.1, 0.1),
        ("lp2x4", 4:6, (0.0, 0.3), (0.5, 1.0), 0.1, 0.1),
        ("lp2x5", 4:6, (0.0, 0.3), (0.45, 1.0), 0.1, 0.1),
        ("lp2x6", 4:6, (0.0, 0.3), (0.4, 1.0), 0.1, 0.1),
        ("lp2x7", 4:6, (0.0, 0.3), (0.35, 1.0), 0.1, 0.1),
        ("lp3x0", 4:6, (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        ("lp3x1", 4:6, (0.0, 0.35), (0.7, 1.0), 0.1, 0.1),
        ("lp3x2", 4:6, (0.0, 0.4), (0.7, 1.0), 0.1, 0.1),
        ("lp3x3", 4:6, (0.0, 0.45), (0.7, 1.0), 0.1, 0.1),
        ("lp3x4", 4:6, (0.0, 0.5), (0.7, 1.0), 0.1, 0.1),
        ("lp3x5", 4:6, (0.0, 0.55), (0.7, 1.0), 0.1, 0.1),
        ("lp3x6", 4:6, (0.0, 0.6), (0.7, 1.0), 0.1, 0.1),
        ("lp3x7", 4:6, (0.0, 0.65), (0.7, 1.0), 0.1, 0.1),
        ("lp4", 4:6, (0.0, 0.5), (0.9, 1.0), 0.1, 0.1),
    ]
    for instance in named_instances
        instancename, wordlengths, passband, stopband, deltapass, deltastop = instance
        for wordlength in wordlengths
            fbands = [passband, stopband]
            abands = [(1-deltapass, 1+deltapass), (0.0, deltastop)]

            model = Model(CPLEX.Optimizer)
            set_silent(model)
            if wordlength >= 10
                if wordlength >= 12
                    if wordlength > 16
                        error("ILP1 not usable due to numerical instability")
                    else
                        set_optimizer_attributes(model, "CPX_PARAM_EPINT" => 10^-8)
                    end
                else
                    set_optimizer_attributes(model, "CPX_PARAM_EPINT" => 10^-7)
                end
            end

            println("$instancename -- $wordlength")

            @time solutions = all_coefficients(model, fbands, abands, wordlength, with_symmetry_breaking=false, verbose=true)
            if !isempty(solutions)
                all_wcpgs = Vector{Float64}()
                all_poles_pos = Vector{Tuple{Float64, Float64}}()
                new_poles = 0
                for solution in solutions
                    tmp = length(all_poles_pos)
                    append!(all_poles_pos, get_poles_coordinates(collect(solution.coefficients_a)./(2^solution.shifts[1])))
                    new_poles = length(all_poles_pos) - tmp
                    push!(all_wcpgs, wcpg_fun(collect(solution.coefficients_a)./(2^solution.shifts[1]), collect(solution.coefficients_b)./(2^solution.shifts[2])))
                    for _ in 2:new_poles
                        push!(all_wcpgs, all_wcpgs[end])
                    end
                end

                p = scatter(first.(all_poles_pos), last.(all_poles_pos), zcolor = all_wcpgs, colormap = :colorwheel, xlims = (-1, 1), ylims = (-1, 1))#, legend = false)
                display(p)
                if wcpgeps
                    savefig((@__DIR__)*"/$(instancename)_$(wordlength)_eps.png")
                else
                    savefig((@__DIR__)*"/$(instancename)_$(wordlength).png")
                end
            end
        end
    end

    return nothing
end
