using CPLEX
using MATLAB
using Polynomials
using JuMP
using Plots

# include("../src/AdderGraphs/src/AdderGraphs.jl")
# using .AdderGraphs
# include("../src/jMCM/src/jMCM.jl")
# using .jMCM
# include("../src/jIIR2AG/src/jIIR2AG.jl")
# using .jIIR2AG
# include("../../src/includepkg.jl")
using AdderGraphs
using jIIR2AG
using jMCM
using jIIR2HW

include("utils.jl")
include("read_write.jl")



"""
    one_step_method(instance)

For a given instance this function generates filter coefficients and produces
FloPoCo calls for the synthesis.
"""
function one_step_method(instance; verbose::Bool=true)
    instancename, wordlength, passband, stopband, deltapass, deltastop = instance
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

    println(" --- $(instancename) --- $wordlength --- ")
    solution = design_second_order_iir(model, fbands, abands, wordlength, avoid_internal_shifts=true, use_big_m=true, presolve_time_sec=600, verbose=verbose)
    # has_values(model) && println(value.(model[:a0]))
    # has_values(model) && println(value.(model[:gab]))

    if solution.shifts != (0,0) && instancename != "compilation"
        wcpg_val = WCPG(solution)
        wcpg_eps_val = WCPG_eps(solution)
        write_flopoco_solution(solution, instancename, wordlength, wcpg_val, wcpg_eps_val)
        varepsilon = coefs_fit_specs_varepsilon(solution.coefficients_b[1]/(1<<solution.shifts[2]),
                                   solution.coefficients_b[2]/(1<<solution.shifts[2]),
                                   solution.coefficients_b[3]/(1<<solution.shifts[2]),
                                   solution.coefficients_a[1]/(1<<solution.shifts[1]),
                                   solution.coefficients_a[2]/(1<<solution.shifts[1]),
                                   get_specifications(fbands, abands))
        write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_shiftandadd",
                            sum(solution.coefficients_b[i] != 0 ? 1 : 0 for i in 1:3)-1+
                                sum(solution.coefficients_a[i] != 0 ? 1 : 0 for i in 1:2)+
                                length(solution.addergraph_b.constants)+length(solution.addergraph_a.constants),
                            varepsilon, wordlength)
    end

    return nothing
end


function float_to_fxp(number_init, wordlength)
    @assert number_init < 1
    return round(Int, number_init*2.0^(wordlength-1))#, RoundDown)
end


function float_to_fxp_both(number_init, wordlength)
    @assert number_init < 1
    return (round(Int, number_init*2.0^(wordlength-1), RoundDown), round(Int, number_init*2.0^(wordlength-1), RoundUp))
end


function truncate_all(b0_init::Float64, b1_init::Float64, b2_init::Float64,
                  a1_init::Float64, a2_init::Float64,
                  specifications::Vector{Tuple{Float64, Float64, Float64}},
                  instancename::String)
    smallestwordlength_found = false

    b0_saved, b1_saved, b2_saved, a1_saved, a2_saved = 0,0,0,0,0
    shiftb_saved, shifta_saved = 0,0
    wordlength_saved = 0

    shifta_init = 0
    shiftb_init = 0
    while max(abs(b0_init),abs(b1_init),abs(b2_init),abs(a1_init),abs(a2_init)) < 1/2
        b0_init,b1_init,b2_init,a1_init,a2_init = (b0_init,b1_init,b2_init,a1_init,a2_init).*2
        shifta_init += 1
        shiftb_init += 1
    end

    while max(abs(b0_init),abs(b1_init),abs(b2_init)) < 1/2
        b0_init,b1_init,b2_init = (b0_init,b1_init,b2_init).*2
        shiftb_init += 1
    end
    while max(abs(a1_init),abs(a2_init)) < 1/2
        a1_init,a2_init = (a1_init,a2_init).*2
        shifta_init += 1
    end

    for wordlength in 4:1:10
        shifta = shifta_init + wordlength-1
        shiftb = shiftb_init + wordlength-1
        b0_inf, b0_sup = float_to_fxp_both(b0_init, wordlength)
        b1_inf, b1_sup = float_to_fxp_both(b1_init, wordlength)
        b2_inf, b2_sup = float_to_fxp_both(b2_init, wordlength)
        a1_inf, a1_sup = float_to_fxp_both(a1_init, wordlength)
        a2_inf, a2_sup = float_to_fxp_both(a2_init, wordlength)

        varepsilon = Inf
        b0, b1, b2, a1, a2 = (b0_inf, b1_inf, b2_inf, a1_inf, a2_inf)

        for b0_tmp in (b0_inf, b0_sup)
            for b1_tmp in (b1_inf, b1_sup)
                for b2_tmp in (b2_inf, b2_sup)
                    for a1_tmp in (a1_inf, a1_sup)
                        for a2_tmp in (a2_inf, a2_sup)
                            if !is_stable(a1_tmp/(1<<shifta), a2_tmp/(1<<shifta))
                                continue
                            end
                            varepsilon_tmp = coefs_fit_specs_varepsilon(
                                                    b0_tmp/(1<<shiftb), b1_tmp/(1<<shiftb), b2_tmp/(1<<shiftb),
                                                    a1_tmp/(1<<shifta), a2_tmp/(1<<shifta),
                                                    specifications
                            )
                            if varepsilon_tmp < varepsilon
                                b0, b1, b2, a1, a2 = (b0_tmp, b1_tmp, b2_tmp, a1_tmp, a2_tmp)
                                varepsilon = varepsilon_tmp
                            end
                        end
                    end
                end
            end
        end

        if !smallestwordlength_found && varepsilon < 1e-7
            smallestwordlength_found = true
            b0_saved, b1_saved, b2_saved, a1_saved, a2_saved = (b0, b1, b2, a1, a2)
            shiftb_saved, shifta_saved = (shiftb, shifta)
            wordlength_saved = wordlength
        end

        println("MCM for $instancename -- $wordlength")
        @time begin
            model = Model(CPLEX.Optimizer)
            set_silent(model)
            solution_b = mcm(model, [b0, b1, b2], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true, avoid_internal_shifts = true)
            model = Model(CPLEX.Optimizer)
            set_silent(model)
            solution_a = mcm(model, [a1, a2], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true, avoid_internal_shifts = true)
        end # @time begin
        if instancename != "compilation"
            write_flopoco_solution([(a1, a2), (b0, b1, b2), solution_a, solution_b], instancename, wordlength, shifta, shiftb)
            write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd", sum(coef != 0 ? 1 : 0 for coef in [a1, a2, b0, b1, b2])-1+length(solution_b)+length(solution_a), varepsilon, wordlength)
            write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain", 0, varepsilon, wordlength)
            write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain0dsp", 0, varepsilon, wordlength)
        end
    end

    if !smallestwordlength_found
        for wordlength in 11:1:32
            shifta = shifta_init + wordlength-1
            shiftb = shiftb_init + wordlength-1
            b0_inf, b0_sup = float_to_fxp_both(b0_init, wordlength)
            b1_inf, b1_sup = float_to_fxp_both(b1_init, wordlength)
            b2_inf, b2_sup = float_to_fxp_both(b2_init, wordlength)
            a1_inf, a1_sup = float_to_fxp_both(a1_init, wordlength)
            a2_inf, a2_sup = float_to_fxp_both(a2_init, wordlength)

            varepsilon = Inf
            b0, b1, b2, a1, a2 = (b0_inf, b1_inf, b2_inf, a1_inf, a2_inf)

            for b0_tmp in (b0_inf, b0_sup)
                for b1_tmp in (b1_inf, b1_sup)
                    for b2_tmp in (b2_inf, b2_sup)
                        for a1_tmp in (a1_inf, a1_sup)
                            for a2_tmp in (a2_inf, a2_sup)
                                if !is_stable(a1_tmp/(1<<shifta), a2_tmp/(1<<shifta))
                                    continue
                                end
                                varepsilon_tmp = coefs_fit_specs_varepsilon(
                                                        b0_tmp/(1<<shiftb), b1_tmp/(1<<shiftb), b2_tmp/(1<<shiftb),
                                                        a1_tmp/(1<<shifta), a2_tmp/(1<<shifta),
                                                        specifications
                                )
                                if varepsilon_tmp < varepsilon
                                    b0, b1, b2, a1, a2 = (b0_tmp, b1_tmp, b2_tmp, a1_tmp, a2_tmp)
                                    varepsilon = varepsilon_tmp
                                end
                            end
                        end
                    end
                end
            end

            if !smallestwordlength_found && varepsilon < 1e-7
                smallestwordlength_found = true
                b0_saved, b1_saved, b2_saved, a1_saved, a2_saved = (b0, b1, b2, a1, a2)
                shiftb_saved, shifta_saved = (shiftb, shifta)
                wordlength_saved = wordlength
                if instancename != "compilation"
                    solution_b = rpag([b0, b1, b2])
                    solution_a = rpag([a1, a2])
                    write_flopoco_solution([(a1, a2), (b0, b1, b2), solution_a, solution_b], instancename, wordlength, shifta, shiftb)
                    write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd", sum(coef != 0 ? 1 : 0 for coef in [a1, a2, b0, b1, b2])-1+length(solution_b)+length(solution_a), varepsilon, wordlength)
                    write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain", 0, varepsilon, wordlength)
                    write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain0dsp", 0, varepsilon, wordlength)
                end
            end
        end
    end

    return b0_saved, b1_saved, b2_saved, a1_saved, a2_saved, shiftb_saved, shifta_saved, wordlength_saved
end


function truncate(b0_init::Float64, b1_init::Float64, b2_init::Float64,
                  a1_init::Float64, a2_init::Float64,
                  specifications::Vector{Tuple{Float64, Float64, Float64}},
                  instancename::String; write_files::Bool=true)
    smallestwordlength_found = false

    b0_saved, b1_saved, b2_saved, a1_saved, a2_saved = 0,0,0,0,0
    shiftb_saved, shifta_saved = 0,0
    wordlength_saved = 0

    shifta_init = 0
    shiftb_init = 0
    while max(abs(b0_init),abs(b1_init),abs(b2_init),abs(a1_init),abs(a2_init)) < 1/2
        b0_init,b1_init,b2_init,a1_init,a2_init = (b0_init,b1_init,b2_init,a1_init,a2_init).*2
        shifta_init += 1
        shiftb_init += 1
    end

    while max(abs(b0_init),abs(b1_init),abs(b2_init)) < 1/2
        b0_init,b1_init,b2_init = (b0_init,b1_init,b2_init).*2
        shiftb_init += 1
    end
    while max(abs(a1_init),abs(a2_init)) < 1/2
        a1_init,a2_init = (a1_init,a2_init).*2
        shifta_init += 1
    end

    for wordlength in 4:1:10
        shifta = shifta_init + wordlength-1
        shiftb = shiftb_init + wordlength-1
        b0 = float_to_fxp(b0_init, wordlength)
        b1 = float_to_fxp(b1_init, wordlength)
        b2 = float_to_fxp(b2_init, wordlength)
        a1 = float_to_fxp(a1_init, wordlength)
        a2 = float_to_fxp(a2_init, wordlength)

        if !is_stable(a1/(1<<shifta), a2/(1<<shifta))
            continue
        end

        varepsilon = coefs_fit_specs_varepsilon(
                                b0/(1<<shiftb), b1/(1<<shiftb), b2/(1<<shiftb),
                                a1/(1<<shifta), a2/(1<<shifta),
                                specifications
        )

        if !smallestwordlength_found && varepsilon < 1e-7
            smallestwordlength_found = true
            b0_saved, b1_saved, b2_saved, a1_saved, a2_saved = (b0, b1, b2, a1, a2)
            shiftb_saved, shifta_saved = (shiftb, shifta)
            wordlength_saved = wordlength
        end

        println("MCM for $instancename -- $wordlength --- varepsilon: $varepsilon")
        @time begin
            model = Model(CPLEX.Optimizer)
            set_silent(model)
            solution_b = mcm(model, [b0, b1, b2], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true, avoid_internal_shifts = true)
            model = Model(CPLEX.Optimizer)
            set_silent(model)
            solution_a = mcm(model, [a1, a2], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true, avoid_internal_shifts = true)
        end # @time begin
        if instancename != "compilation"
            write_files && write_flopoco_solution([(a1, a2), (b0, b1, b2), solution_a, solution_b], instancename, wordlength, shifta, shiftb)
            write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd", sum(coef != 0 ? 1 : 0 for coef in [a1, a2, b0, b1, b2])-1+length(solution_b)+length(solution_a), varepsilon, wordlength)
            write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain", 0, varepsilon, wordlength)
            write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain0dsp", 0, varepsilon, wordlength)
        end
    end

    if !smallestwordlength_found
        for wordlength in 11:1:32
            shifta = shifta_init + wordlength-1
            shiftb = shiftb_init + wordlength-1
            b0 = float_to_fxp(b0_init, wordlength)
            b1 = float_to_fxp(b1_init, wordlength)
            b2 = float_to_fxp(b2_init, wordlength)
            a1 = float_to_fxp(a1_init, wordlength)
            a2 = float_to_fxp(a2_init, wordlength)

            if !is_stable(a1/(1<<shifta), a2/(1<<shifta))
                continue
            end

            varepsilon = coefs_fit_specs_varepsilon(
                                    b0/(1<<shiftb), b1/(1<<shiftb), b2/(1<<shiftb),
                                    a1/(1<<shifta), a2/(1<<shifta),
                                    specifications)

            println("\t$instancename -- $wordlength --- varepsilon: $varepsilon")
            if !smallestwordlength_found && varepsilon < 1e-7
                smallestwordlength_found = true
                b0_saved, b1_saved, b2_saved, a1_saved, a2_saved = (b0, b1, b2, a1, a2)
                shiftb_saved, shifta_saved = (shiftb, shifta)
                wordlength_saved = wordlength
                if instancename != "compilation"
                    solution_b = rpag([b0, b1, b2])
                    solution_a = rpag([a1, a2])
                    write_files && write_flopoco_solution([(a1, a2), (b0, b1, b2), solution_a, solution_b], instancename, wordlength, shifta, shiftb)
                    write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd", sum(coef != 0 ? 1 : 0 for coef in [a1, a2, b0, b1, b2])-1+length(solution_b)+length(solution_a), varepsilon, wordlength)
                    write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain", 0, varepsilon, wordlength)
                    write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain0dsp", 0, varepsilon, wordlength)
                end
            end
        end
    end

    return b0_saved, b1_saved, b2_saved, a1_saved, a2_saved, shiftb_saved, shifta_saved, wordlength_saved
end




function benchmarking()
    verbose=true
    named_instances = [
        ("compilation", [8], (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        #("hp1", 4:10, (0.7, 1.0), (0.0, 0.3), 0.1, 0.1),
        ("lp4", 4:10, (0.0, 0.5), (0.9, 1.0), 0.1, 0.1),
        ("lp1x0", 4:10, (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        ("lp1x1", 4:10, (0.0, 0.3), (0.7, 1.0), 0.09, 0.09),
        ("lp1x2", 4:10, (0.0, 0.3), (0.7, 1.0), 0.08, 0.08),
        ("lp1x3", 4:10, (0.0, 0.3), (0.7, 1.0), 0.07, 0.07),
        ("lp1x4", 4:10, (0.0, 0.3), (0.7, 1.0), 0.06, 0.06),
        ("lp1x5", 4:10, (0.0, 0.3), (0.7, 1.0), 0.05, 0.05),
        ("lp1x6", 4:10, (0.0, 0.3), (0.7, 1.0), 0.04, 0.04),
        ("lp1x7", 4:10, (0.0, 0.3), (0.7, 1.0), 0.03, 0.03),
        ("lp1x8", 4:10, (0.0, 0.3), (0.7, 1.0), 0.02, 0.02),
        ("lp1x9", 4:10, (0.0, 0.3), (0.7, 1.0), 0.01, 0.01),
        ("lp2x0", 4:10, (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        ("lp2x1", 4:10, (0.0, 0.3), (0.65, 1.0), 0.1, 0.1),
        ("lp2x2", 4:10, (0.0, 0.3), (0.6, 1.0), 0.1, 0.1),
        ("lp2x3", 4:10, (0.0, 0.3), (0.55, 1.0), 0.1, 0.1),
        ("lp2x4", 4:10, (0.0, 0.3), (0.5, 1.0), 0.1, 0.1),
        ("lp2x5", 4:10, (0.0, 0.3), (0.45, 1.0), 0.1, 0.1),
        ("lp2x6", 4:10, (0.0, 0.3), (0.4, 1.0), 0.1, 0.1),
        ("lp2x7", 4:10, (0.0, 0.3), (0.35, 1.0), 0.1, 0.1),
        ("lp3x0", 4:10, (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        ("lp3x1", 4:10, (0.0, 0.35), (0.7, 1.0), 0.1, 0.1),
        ("lp3x2", 4:10, (0.0, 0.4), (0.7, 1.0), 0.1, 0.1),
        ("lp3x3", 4:10, (0.0, 0.45), (0.7, 1.0), 0.1, 0.1),
        ("lp3x4", 4:10, (0.0, 0.5), (0.7, 1.0), 0.1, 0.1),
        ("lp3x5", 4:10, (0.0, 0.55), (0.7, 1.0), 0.1, 0.1),
        ("lp3x6", 4:10, (0.0, 0.6), (0.7, 1.0), 0.1, 0.1),
        ("lp3x7", 4:10, (0.0, 0.65), (0.7, 1.0), 0.1, 0.1)
    ]

    for instance in named_instances
        instancename, wordlengths, passband, stopband, deltapass, deltastop = instance
        for wordlength in wordlengths
            one_step_method((instancename, wordlength, passband, stopband, deltapass, deltastop), verbose=verbose)
            verbose && println()
        end
    end

    verbose && println()

    design_method = "ellip"
    match_exactly = "both"
    for instance in named_instances
        instancename, wordlengths, passband, stopband, deltapass, deltastop = instance
        fbands = [passband, stopband]
        abands = [(1-deltapass, 1+deltapass), (0.0, deltastop)]
        specifications = get_specifications(fbands, abands)
        eval_string("$instancename = designfilt('lowpassiir','DesignMethod','$design_method','PassbandFrequency',$(passband[2]),'StopbandFrequency',$(stopband[1]),'MatchExactly','$match_exactly','PassbandRipple',$(20*log10(1/(1.0-deltapass))),'StopbandAttenuation',$(20*log10(1/deltastop)));
            coefs_str = num2hex($instancename.Coefficients);
            nb_sections = size($instancename.Coefficients, 1);")

        if @mget(nb_sections) != 1
            println("No second order IIR filter using $design_method ($match_exactly) for $instancename")
            continue
        end

        eval_string("coef = coefs_str(1,:);")
        b0_hex = @mget(coef)
        b0 = reinterpret(Float64, reverse(hex2bytes(b0_hex)))[1]
        eval_string("coef = coefs_str(2,:);")
        b1_hex = @mget(coef)
        b1 = reinterpret(Float64, reverse(hex2bytes(b1_hex)))[1]
        eval_string("coef = coefs_str(3,:);")
        b2_hex = @mget(coef)
        b2 = reinterpret(Float64, reverse(hex2bytes(b2_hex)))[1]
        eval_string("coef = coefs_str(5,:);")
        a1_hex = @mget(coef)
        a1 = reinterpret(Float64, reverse(hex2bytes(a1_hex)))[1]
        eval_string("coef = coefs_str(6,:);")
        a2_hex = @mget(coef)
        a2 = reinterpret(Float64, reverse(hex2bytes(a2_hex)))[1]

        # Check that they fit
        varepsilon = coefs_fit_specs_varepsilon(b0,b1,b2,a1,a2,specifications)
        if varepsilon != 0.0
            println("Matlab error for $instancename, varepsilon = $varepsilon")
        end
        if instancename != "compilation"
            write_table_epsilon("$(instancename)_dw00_cw00_fixiir", 0, varepsilon, 0)
        end

        # Save coefficients in a file
        if instancename != "compilation"
            write_flopoco_solution(b0_hex, b1_hex, b2_hex, a1_hex, a2_hex, instancename)
        end

        # Truncate coefficients
        b0,b1,b2,a1,a2,shiftb,shifta,wordlength = truncate(b0,b1,b2,a1,a2,specifications,instancename)
        verbose && println("wordlength: $(wordlength)")
        verbose && println("coefficients_b: $(b0), $(b1), $(b2)")
        verbose && println("coefficients_a: $(a1), $(a2)")
        verbose && println("shiftb: $shiftb")
        verbose && println("shifta: $shifta")
        verbose && println("eps: $(varepsilon)")
        verbose && println()
    end

    return nothing
end










function special_cases()
    verbose=true
    named_instances = [
        ("compilation", (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        ("lp1x5", (0.0, 0.3), (0.7, 1.0), 0.05, 0.05),
        ("lp2x3", (0.0, 0.3), (0.55, 1.0), 0.1, 0.1),
        ("lp3x3", (0.0, 0.45), (0.7, 1.0), 0.1, 0.1),
    ]

    design_method = "ellip"
    match_exactly = "both"
    for instance in named_instances
        instancename, passband, stopband, deltapass, deltastop = instance
        fbands = [passband, stopband]
        abands = [(1-deltapass, 1+deltapass), (0.0, deltastop)]
        specifications = get_specifications(fbands, abands)
        eval_string("$instancename = designfilt('lowpassiir','DesignMethod','$design_method','PassbandFrequency',$(passband[2]),'StopbandFrequency',$(stopband[1]),'MatchExactly','$match_exactly','PassbandRipple',$(20*log10(1/(1.0-deltapass))),'StopbandAttenuation',$(20*log10(1/deltastop)));
            coefs_str = num2hex($instancename.Coefficients);
            nb_sections = size($instancename.Coefficients, 1);")

        println("Before while: $deltapass")
        while @mget(nb_sections) != 1
            deltapass += 1e-5
            deltastop += 1e-5
            eval_string("$instancename = designfilt('lowpassiir','DesignMethod','$design_method','PassbandFrequency',$(passband[2]),'StopbandFrequency',$(stopband[1]),'MatchExactly','$match_exactly','PassbandRipple',$(20*log10(1/(1.0-deltapass))),'StopbandAttenuation',$(20*log10(1/deltastop)));
                coefs_str = num2hex($instancename.Coefficients);
                nb_sections = size($instancename.Coefficients, 1);")
        end
        println(deltapass)
        if @mget(nb_sections) != 1
            println("No second order IIR filter using $design_method ($match_exactly) for $instancename")
            continue
        end

        eval_string("coef = coefs_str(1,:);")
        b0_hex = @mget(coef)
        b0 = reinterpret(Float64, reverse(hex2bytes(b0_hex)))[1]
        eval_string("coef = coefs_str(2,:);")
        b1_hex = @mget(coef)
        b1 = reinterpret(Float64, reverse(hex2bytes(b1_hex)))[1]
        eval_string("coef = coefs_str(3,:);")
        b2_hex = @mget(coef)
        b2 = reinterpret(Float64, reverse(hex2bytes(b2_hex)))[1]
        eval_string("coef = coefs_str(5,:);")
        a1_hex = @mget(coef)
        a1 = reinterpret(Float64, reverse(hex2bytes(a1_hex)))[1]
        eval_string("coef = coefs_str(6,:);")
        a2_hex = @mget(coef)
        a2 = reinterpret(Float64, reverse(hex2bytes(a2_hex)))[1]

        # Check that they fit
        varepsilon = coefs_fit_specs_varepsilon(b0,b1,b2,a1,a2,specifications)
        if varepsilon != 0.0
            println("Matlab error for $instancename, varepsilon = $varepsilon")
        end
        if instancename != "compilation"
            write_table_epsilon("$(instancename)_dw00_cw00_fixiir", 0, varepsilon, 0)
        end

        # Save coefficients in a file
        if instancename != "compilation"
            write_flopoco_solution(b0_hex, b1_hex, b2_hex, a1_hex, a2_hex, instancename)
        end

        # Truncate coefficients
        b0,b1,b2,a1,a2,shiftb,shifta,wordlength = truncate(b0,b1,b2,a1,a2,specifications,instancename)
        varepsilon = coefs_fit_specs_varepsilon(
                                b0/(1<<shiftb), b1/(1<<shiftb), b2/(1<<shiftb),
                                a1/(1<<shifta), a2/(1<<shifta),
                                specifications
        )
        verbose && println("wordlength: $(wordlength)")
        verbose && println("coefficients_b: $(b0), $(b1), $(b2)")
        verbose && println("coefficients_a: $(a1), $(a2)")
        verbose && println("shiftb: $shiftb")
        verbose && println("shifta: $shifta")
        verbose && println("eps: $(varepsilon)")
    end

    return nothing
end



function main_hp0()
    verbose = true
    instancename, wordlengths, passband, stopband, deltapass, deltastop = ("hp0", 4:10, (0.041, 1.0), (0.0, 0.0), 0.016, 0.0)
    for wordlength in wordlengths
        one_step_method((instancename, wordlength, passband, stopband, deltapass, deltastop), verbose=verbose)
        verbose && println()
    end

    return nothing
end














function get_all_solutions()
    verbose = true
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

            @time solutions = all_coefficients(model, fbands, abands, wordlength, with_symmetry_breaking=false, verbose=verbose)
            if !isempty(solutions)
                open((@__DIR__)*"/solutions/$(instancename)_$(wordlength)_all_solutions.txt", "w") do io
                    write(io, "a1;a2;b0;b1;b2;shifta;shiftb\n")
                    for solution in solutions
                        write(io, "$(solution.coefficients_a[1]);$(solution.coefficients_a[2]);$(solution.coefficients_b[1]);$(solution.coefficients_b[2]);$(solution.coefficients_b[3]);$(solution.shifts[1]);$(solution.shifts[2])\n")
                    end
                end
            end
        end
    end

    return nothing
end





function WCPG(a::Vector{Float64}, b::Vector{Float64}, solution=nothing)::Float64
    W = zeros(1)
    if ccall((:WCPG_tf, "libwcpg"), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cuintmax_t, Cuintmax_t, Cint), W, b, a, length(b), length(a), 0) == 0
        @warn "WCPG computation errored"
        if !isnothing(solution)
            a1, a2 = solution[1]
            b0, b1, b2 = solution[2]
            shifta, shiftb = solution[3], solution[4]
            C = Matrix{BigInt}([b1*2^(shifta)-b0*a1 b2*2^(shifta)-b0*a2])
            A = Matrix{BigInt}([-a1 -a2; 2^(shifta) 0])
            B = Vector{BigInt}([1, 0])
            D = Vector{BigInt}([b0])
            N = 1001
            W[1] = (BigInt(2)^(shiftb)*sum((abs.(C*(A^k)*B))[1]*(BigInt(2)^(shifta))^(N-k) for k in 0:N)+abs(D[1])*(BigInt(2)^(shifta)*BigInt(2)^(shiftb))*(BigInt(2)^(shifta))^N)/(BigInt(2)^(shiftb)*BigInt(2)^(shiftb)*BigInt(2)^(shifta)*(BigInt(2)^(shifta))^N)
        end
    end
    return W[1]
end

function WCPG(solution::Iir2AdderGraphs)::Float64
    return WCPG(collect(solution.coefficients_a)./(2^solution.shifts[1]), collect(solution.coefficients_b)./(2^solution.shifts[2]),
        (solution.coefficients_a, solution.coefficients_b, solution.shifts[1], solution.shifts[2]))
end


function WCPG_eps(a::Vector{Float64}, b::Vector{Float64}, solution=nothing)::Float64
    W = zeros(1)
    if ccall((:WCPG_tf, "libwcpg"), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cuintmax_t, Cuintmax_t, Cint), W, [1.0], a, 1, length(a), 0) == 0
        @warn "WCPG_eps computation errored"
        if !isnothing(solution)
            a1, a2 = solution[1]
            b0, b1, b2 = 1,0,0
            shifta, shiftb = solution[3], 0
            C = Matrix{BigInt}([b1*2^(shifta)-b0*a1 b2*2^(shifta)-b0*a2])
            A = Matrix{BigInt}([-a1 -a2; 2^(shifta) 0])
            B = Vector{BigInt}([1, 0])
            D = Vector{BigInt}([b0])
            N = 1001
            W[1] = (BigInt(2)^(shiftb)*sum((abs.(C*(A^k)*B))[1]*(BigInt(2)^(shifta))^(N-k) for k in 0:N)+abs(D[1])*(BigInt(2)^(shifta)*BigInt(2)^(shiftb))*(BigInt(2)^(shifta))^N)/(BigInt(2)^(shiftb)*BigInt(2)^(shiftb)*BigInt(2)^(shifta)*(BigInt(2)^(shifta))^N)
        end
    end
    return W[1]
end

function WCPG_eps(solution::Iir2AdderGraphs)::Float64
    return WCPG_eps(collect(solution.coefficients_a)./(2^solution.shifts[1]), collect(solution.coefficients_b)./(2^solution.shifts[2]),
        (solution.coefficients_a, solution.coefficients_b, solution.shifts[1], solution.shifts[2]))
end

function get_poles_coordinates(a::Vector{Float64})
    poles_complex = roots(Polynomial(push!(reverse(a),1)))

    return [(real(pole_complex), imag(pole_complex)) for pole_complex in poles_complex]
end

function plot_all_wcpg()
    for filename in readdir((@__DIR__)*"/solutions/")
        if length(filename) <= 5 || filename[(end-3):end] != ".txt"
            continue
        end
        println("$filename")
        solutions = Vector{Tuple{Tuple{Int, Int}, Tuple{Int, Int, Int}, Int, Int}}()

        open((@__DIR__)*"/solutions/"*filename, "r") do f
            lines = readlines(f)
            for line in lines[2:end]
                splitted_line = parse.(Int, split(line, ";"))
                push!(solutions, ((splitted_line[1], splitted_line[2]), (splitted_line[3], splitted_line[4], splitted_line[5]), splitted_line[6], splitted_line[7]))
            end
        end
        if !isempty(solutions)
            all_wcpgs = Vector{Float64}()
            all_wcpgs_eps = Vector{Float64}()
            all_poles_pos = Vector{Tuple{Float64, Float64}}()
            new_poles = 0
            for solution in solutions
                tmp = length(all_poles_pos)
                append!(all_poles_pos, get_poles_coordinates(collect(solution[1])./(2^solution[3])))
                new_poles = length(all_poles_pos) - tmp
                push!(all_wcpgs, WCPG(collect(solution[1])./(2^solution[3]), collect(solution[2])./(2^solution[4]), solution))
                push!(all_wcpgs_eps, WCPG_eps(collect(solution[1])./(2^solution[3]), collect(solution[2])./(2^solution[4]), solution))
                for _ in 2:new_poles
                    push!(all_wcpgs, all_wcpgs[end])
                    push!(all_wcpgs_eps, all_wcpgs_eps[end])
                end
            end

            xy = [(cos(t), sin(t)) for t in 0:0.01:(2*pi)]
            p = scatter(first.(all_poles_pos), last.(all_poles_pos), zcolor = all_wcpgs, color = :Reds, xlims = (-1, 1), ylims = (-1, 1), framestyle = :origin, aspect_ratio = :equal)#, legend = false)
            plot!(p, xy, legend = false, color = :black)
            display(p)
            savefig((@__DIR__)*"/solutions/wcpg/svg/$(filename[1:(end-4)]).svg")
            savefig((@__DIR__)*"/solutions/wcpg/png/$(filename[1:(end-4)]).png")
            p = scatter(first.(all_poles_pos), last.(all_poles_pos), zcolor = all_wcpgs_eps, color = :Reds, xlims = (-1, 1), ylims = (-1, 1), framestyle = :origin, aspect_ratio = :equal)#, legend = false)
            plot!(p, xy, legend = false, color = :black)
            display(p)
            savefig((@__DIR__)*"/solutions/wcpg/svg/$(filename[1:(end-4)])_eps.svg")
            savefig((@__DIR__)*"/solutions/wcpg/png/$(filename[1:(end-4)])_eps.png")
        end
    end

    return nothing
end


function plot_best_wcpg()
    for filename in readdir((@__DIR__)*"/solutions/")
        if length(filename) <= 5 || filename[(end-3):end] != ".txt"
            continue
        end
        println("$filename")
        solutions = Vector{Tuple{Tuple{Int, Int}, Tuple{Int, Int, Int}, Int, Int, Int, Int}}()

        open((@__DIR__)*"/solutions/"*filename, "r") do f
            lines = readlines(f)
            for line in lines[2:end]
                splitted_line = parse.(Int, split(line, ";"))
                push!(solutions, ((splitted_line[1], splitted_line[2]), (splitted_line[3], splitted_line[4], splitted_line[5]), splitted_line[6], splitted_line[7], splitted_line[8], splitted_line[9]))
            end
        end
        println("\t$(length(solutions)) solutions")
        if !isempty(solutions)
            best_mcm = typemax(Int)
            for solution in solutions
                coefs_a = sort(collect(solution[1]))
                coefs_b = sort(collect(solution[2]))
                if solution[5]+solution[6]+sum(coefs_b[i] != 0 ? 1 : 0 for i in 1:3)-1+sum(coefs_a[i] != 0 ? 1 : 0 for i in 1:2) < best_mcm
                    best_mcm = solution[5]+solution[6]+sum(coefs_b[i] != 0 ? 1 : 0 for i in 1:3)-1+sum(coefs_a[i] != 0 ? 1 : 0 for i in 1:2)
                end
            end
            all_wcpgs = Vector{Float64}()
            all_wcpgs_eps = Vector{Float64}()
            all_poles_pos = Vector{Tuple{Float64, Float64}}()
            new_poles = 0
            for solution in solutions
                coefs_a = sort(collect(solution[1]))
                coefs_b = sort(collect(solution[2]))
                if solution[5]+solution[6]+sum(coefs_b[i] != 0 ? 1 : 0 for i in 1:3)-1+sum(coefs_a[i] != 0 ? 1 : 0 for i in 1:2) > best_mcm
                    continue
                end
                tmp = length(all_poles_pos)
                append!(all_poles_pos, get_poles_coordinates(collect(solution[1])./(2^solution[3])))
                new_poles = length(all_poles_pos) - tmp
                push!(all_wcpgs, WCPG(collect(solution[1])./(2^solution[3]), collect(solution[2])./(2^solution[4]), solution))
                push!(all_wcpgs_eps, WCPG_eps(collect(solution[1])./(2^solution[3]), collect(solution[2])./(2^solution[4]), solution))
                for _ in 2:new_poles
                    push!(all_wcpgs, all_wcpgs[end])
                    push!(all_wcpgs_eps, all_wcpgs_eps[end])
                end
            end

            xy = [(cos(t), sin(t)) for t in 0:0.01:(2*pi)]
            println("\t$(length(all_poles_pos)) poles")
            p = scatter(first.(all_poles_pos), last.(all_poles_pos), zcolor = all_wcpgs, color = :Reds, xlims = (-1, 1), ylims = (-1, 1), framestyle = :origin, aspect_ratio = :equal)#, legend = false)
            plot!(p, xy, legend = false, color = :black)
            display(p)
            savefig((@__DIR__)*"/solutions/wcpg/lowest_mcm/svg/$(filename[1:(end-4)]).svg")
            savefig((@__DIR__)*"/solutions/wcpg/lowest_mcm/png/$(filename[1:(end-4)]).png")
            p = scatter(first.(all_poles_pos), last.(all_poles_pos), zcolor = all_wcpgs_eps, color = :Reds, xlims = (-1, 1), ylims = (-1, 1), framestyle = :origin, aspect_ratio = :equal)#, legend = false)
            plot!(p, xy, legend = false, color = :black)
            display(p)
            savefig((@__DIR__)*"/solutions/wcpg/lowest_mcm/svg/$(filename[1:(end-4)])_eps.svg")
            savefig((@__DIR__)*"/solutions/wcpg/lowest_mcm/png/$(filename[1:(end-4)])_eps.png")
        end
    end

    return nothing
end




function compute_all_mcm()
    mcm_values = Dict{Tuple{Int, Int, Int}, Int}()
    for filename in readdir((@__DIR__)*"/solutions/")
        if length(filename) <= 5 || filename[(end-3):end] != ".txt"
            continue
        end
        println("$filename")
        solutions = Vector{Tuple{Tuple{Int, Int}, Tuple{Int, Int, Int}, Int, Int}}()

        open((@__DIR__)*"/solutions/$(filename)", "r") do f
            lines = readlines(f)
            for line in lines[2:end]
                splitted_line = parse.(Int, split(line, ";"))
                push!(solutions, ((splitted_line[1], splitted_line[2]), (splitted_line[3], splitted_line[4], splitted_line[5]), splitted_line[6], splitted_line[7]))
            end
        end
        if !isempty(solutions)
            open((@__DIR__)*"/solutions/$(filename)", "w") do io
                write(io, "a1;a2;b0;b1;b2;shifta;shiftb;mcma;mcmb\n")
                for solution in solutions
                    coefs_a = sort(collect(solution[1]))
                    coefs_b = sort(collect(solution[2]))
                    model = Model(CPLEX.Optimizer)
                    mcma = get!(mcm_values, (0, coefs_a[1], coefs_a[2]), length(mcm(model, coefs_a, avoid_internal_shifts=true)))
                    model = Model(CPLEX.Optimizer)
                    mcmb = get!(mcm_values, (coefs_b[1], coefs_b[2], coefs_b[3]), length(mcm(model, coefs_b, avoid_internal_shifts=true)))
                    write(io, "$(solution[1][1]);$(solution[1][2]);$(solution[2][1]);$(solution[2][2]);$(solution[2][3]);$(solution[3]);$(solution[4]);$(mcma);$(mcmb)\n")
                end
            end
        end
    end

    return nothing
end



function plot_all_mcm()
    for filename in readdir((@__DIR__)*"/solutions/")
        if length(filename) <= 5 || filename[(end-3):end] != ".txt"
            continue
        end
        println("$filename")
        solutions = Vector{Tuple{Tuple{Int, Int}, Tuple{Int, Int, Int}, Int, Int, Int, Int}}()

        open((@__DIR__)*"/solutions/$(filename)", "r") do f
            lines = readlines(f)
            for line in lines[2:end]
                splitted_line = parse.(Int, split(line, ";"))
                push!(solutions, ((splitted_line[1], splitted_line[2]), (splitted_line[3], splitted_line[4], splitted_line[5]), splitted_line[6], splitted_line[7], splitted_line[8], splitted_line[9]))
            end
        end
        if !isempty(solutions)
            all_mcm = Vector{Int}()
            all_poles_pos = Vector{Tuple{Float64, Float64}}()
            new_poles = 0
            for solution in solutions
                tmp = length(all_poles_pos)
                append!(all_poles_pos, get_poles_coordinates(collect(solution[1])./(2^solution[3])))
                new_poles = length(all_poles_pos) - tmp
                coefs_a = sort(collect(solution[1]))
                coefs_b = sort(collect(solution[2]))
                push!(all_mcm, solution[5]+solution[6]+sum(coefs_b[i] != 0 ? 1 : 0 for i in 1:3)-1+sum(coefs_a[i] != 0 ? 1 : 0 for i in 1:2))
                for _ in 2:new_poles
                    push!(all_mcm, all_mcm[end])
                end
            end

            xy = [(cos(t), sin(t)) for t in 0:0.01:(2*pi)]
            p = scatter(first.(all_poles_pos), last.(all_poles_pos), zcolor = all_mcm, color = :Reds, xlims = (-1, 1), ylims = (-1, 1), framestyle = :origin, aspect_ratio = :equal)#, legend = false)
            plot!(p, xy, legend = false, color = :black)
            display(p)
            savefig((@__DIR__)*"/solutions/mcm/svg/$(filename[1:(end-4)])_mcm.svg")
            savefig((@__DIR__)*"/solutions/mcm/png/$(filename[1:(end-4)])_mcm.png")
        end
    end

    return nothing
end





function plot_wcpg()
    verbose=true
    named_instances = [
        #("compilation", [8], (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        #("hp1", 4:10, (0.7, 1.0), (0.0, 0.3), 0.1, 0.1),
        ("hp0", 4:8, (0.041, 1.0), (0.0, 0.0), 0.016, 0.0),
        ("lp4", 4:10, (0.0, 0.5), (0.9, 1.0), 0.1, 0.1),
        ("lp1x0", 4:10, (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        ("lp1x1", 4:10, (0.0, 0.3), (0.7, 1.0), 0.09, 0.09),
        ("lp1x2", 4:10, (0.0, 0.3), (0.7, 1.0), 0.08, 0.08),
        ("lp1x3", 4:10, (0.0, 0.3), (0.7, 1.0), 0.07, 0.07),
        ("lp1x4", 4:10, (0.0, 0.3), (0.7, 1.0), 0.06, 0.06),
        ("lp1x5", 4:10, (0.0, 0.3), (0.7, 1.0), 0.05, 0.05),
        ("lp1x6", 4:10, (0.0, 0.3), (0.7, 1.0), 0.04, 0.04),
        ("lp1x7", 4:10, (0.0, 0.3), (0.7, 1.0), 0.03, 0.03),
        ("lp1x8", 4:10, (0.0, 0.3), (0.7, 1.0), 0.02, 0.02),
        ("lp1x9", 4:10, (0.0, 0.3), (0.7, 1.0), 0.01, 0.01),
        ("lp2x0", 4:10, (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        ("lp2x1", 4:10, (0.0, 0.3), (0.65, 1.0), 0.1, 0.1),
        ("lp2x2", 4:10, (0.0, 0.3), (0.6, 1.0), 0.1, 0.1),
        ("lp2x3", 4:10, (0.0, 0.3), (0.55, 1.0), 0.1, 0.1),
        ("lp2x4", 4:10, (0.0, 0.3), (0.5, 1.0), 0.1, 0.1),
        ("lp2x5", 4:10, (0.0, 0.3), (0.45, 1.0), 0.1, 0.1),
        ("lp2x6", 4:10, (0.0, 0.3), (0.4, 1.0), 0.1, 0.1),
        ("lp2x7", 4:10, (0.0, 0.3), (0.35, 1.0), 0.1, 0.1),
        ("lp3x0", 4:10, (0.0, 0.3), (0.7, 1.0), 0.1, 0.1),
        ("lp3x1", 4:10, (0.0, 0.35), (0.7, 1.0), 0.1, 0.1),
        ("lp3x2", 4:10, (0.0, 0.4), (0.7, 1.0), 0.1, 0.1),
        ("lp3x3", 4:10, (0.0, 0.45), (0.7, 1.0), 0.1, 0.1),
        ("lp3x4", 4:10, (0.0, 0.5), (0.7, 1.0), 0.1, 0.1),
        ("lp3x5", 4:10, (0.0, 0.55), (0.7, 1.0), 0.1, 0.1),
        ("lp3x6", 4:10, (0.0, 0.6), (0.7, 1.0), 0.1, 0.1),
        ("lp3x7", 4:10, (0.0, 0.65), (0.7, 1.0), 0.1, 0.1)
    ]

    for instance in named_instances
        instancename, wordlengths, passband, stopband, deltapass, deltastop = instance
        fbands = [passband, stopband]
        abands = [(1-deltapass, 1+deltapass), (0.0, deltastop)]
        for wordlength in wordlengths
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

            println(" --- $(instancename) --- $wordlength --- ")
            solution = design_second_order_iir(model, fbands, abands, wordlength, avoid_internal_shifts=true, use_big_m=true, presolve_time_sec=600, verbose=verbose)
            # has_values(model) && println(value.(model[:a0]))
            # has_values(model) && println(value.(model[:gab]))

            if solution.shifts != (0,0) && instancename != "compilation"
                nb_wcpg = 1001
                wcpg_n = zeros(nb_wcpg)
                wcpg_eps_n = zeros(nb_wcpg)
                a1, a2 = solution.coefficients_a
                b0, b1, b2 = solution.coefficients_b
                shifta, shiftb = solution.shifts
                C = Matrix{BigInt}([b1*2^(shifta)-b0*a1 b2*2^(shifta)-b0*a2])
                A = Matrix{BigInt}([-a1 -a2; 2^(shifta) 0])
                B = Vector{BigInt}([1, 0])
                D = Vector{BigInt}([b0])
                wcpg_n[0+1] = (BigInt(2)^(shiftb)*(abs.(C*B))[1]+abs(D[1])*(BigInt(2)^(shifta)*BigInt(2)^(shiftb)))/(BigInt(2)^(shiftb)*BigInt(2)^(shiftb)*BigInt(2)^(shifta))
                for N in 1:(nb_wcpg-1)
                    wcpg_n[N+1] = (BigInt(2)^(shiftb)*sum((abs.(C*(A^k)*B))[1]*(BigInt(2)^(shifta))^(N-k) for k in 0:N)+abs(D[1])*(BigInt(2)^(shifta)*BigInt(2)^(shiftb))*(BigInt(2)^(shifta))^N)/(BigInt(2)^(shiftb)*BigInt(2)^(shiftb)*BigInt(2)^(shifta)*(BigInt(2)^(shifta))^N)
                end
                a1, a2 = solution.coefficients_a
                b0, b1, b2 = (1, 0, 0)
                shifta, shiftb = solution.shifts[1], 0
                C = Matrix{BigInt}([b1*2^(shifta)-b0*a1 b2*2^(shifta)-b0*a2])
                A = Matrix{BigInt}([-a1 -a2; 2^(shifta) 0])
                B = Vector{BigInt}([1, 0])
                D = Vector{BigInt}([b0])
                wcpg_eps_n[0+1] = (BigInt(2)^(shiftb)*(abs.(C*B))[1]+abs(D[1])*(BigInt(2)^(shifta)*BigInt(2)^(shiftb)))/(BigInt(2)^(shiftb)*BigInt(2)^(shiftb)*BigInt(2)^(shifta))
                for N in 1:(nb_wcpg-1)
                    wcpg_eps_n[N+1] = (BigInt(2)^(shiftb)*sum((abs.(C*(A^k)*B))[1]*(BigInt(2)^(shifta))^(N-k) for k in 0:N)+abs(D[1])*(BigInt(2)^(shifta)*BigInt(2)^(shiftb))*(BigInt(2)^(shifta))^N)/(BigInt(2)^(shiftb)*BigInt(2)^(shiftb)*BigInt(2)^(shifta)*(BigInt(2)^(shifta))^N)
                end
                current_n = 10
                while current_n < nb_wcpg
                    p = plot(0:current_n, wcpg_n[1:(current_n+1)])
                    plot!(p, legend=:bottomright)
                    display(p)
                    savefig((@__DIR__)*"/solutions/wcpg/N/svg/$(instancename)_$(wordlength)_wcpg_N$(current_n).svg")
                    savefig((@__DIR__)*"/solutions/wcpg/N/png/$(instancename)_$(wordlength)_wcpg_N$(current_n).png")
                    p = plot(0:current_n, wcpg_eps_n[1:(current_n+1)])
                    plot!(p, legend=:bottomright)
                    display(p)
                    savefig((@__DIR__)*"/solutions/wcpg/N/svg/$(instancename)_$(wordlength)_wcpg_eps_N$(current_n).svg")
                    savefig((@__DIR__)*"/solutions/wcpg/N/png/$(instancename)_$(wordlength)_wcpg_eps_N$(current_n).png")
                    current_n *= 10
                end
            end
        end
    end
    return nothing
end


function maxdiff(b0_init,b1_init,b2_init,a1_init,a2_init, b0, b1, b2, a1, a2)
    varepsilon = 0.0
    a0 = 1.0
    a0_init = 1.0
    x = 0.0:0.01:1.0
    y1 = [abs((b0*exp(2*pi*im*val)+b1*exp(im*pi*val)+b2)/(a0*exp(2*pi*im*val)+a1*exp(im*pi*val)+a2)) for val in x];
    y2 = [abs((b0_init*exp(2*pi*im*val)+b1_init*exp(im*pi*val)+b2_init)/(a0_init*exp(2*pi*im*val)+a1_init*exp(im*pi*val)+a2_init)) for val in x];
    println(y1)
    println(y2)
    varepsilon = maximum(abs.(y1 .- y2))
    return varepsilon
end

function truncate(b0_init::Float64, b1_init::Float64, b2_init::Float64,
                  a1_init::Float64, a2_init::Float64,
                  instancename::String; write_files::Bool=true)
    smallestwordlength_found = false

    b0_saved, b1_saved, b2_saved, a1_saved, a2_saved = 0,0,0,0,0
    shiftb_saved, shifta_saved = 0,0
    wordlength_saved = 0

    shifta_init = 0
    shiftb_init = 0
    while max(abs(b0_init),abs(b1_init),abs(b2_init),abs(a1_init),abs(a2_init)) < 1/2
        b0_init,b1_init,b2_init,a1_init,a2_init = (b0_init,b1_init,b2_init,a1_init,a2_init).*2
        shifta_init += 1
        shiftb_init += 1
    end

    while max(abs(b0_init),abs(b1_init),abs(b2_init)) < 1/2
        b0_init,b1_init,b2_init = (b0_init,b1_init,b2_init).*2
        shiftb_init += 1
    end
    while max(abs(a1_init),abs(a2_init)) < 1/2
        a1_init,a2_init = (a1_init,a2_init).*2
        shifta_init += 1
    end

    while max(abs(b0_init),abs(b1_init),abs(b2_init),abs(a1_init),abs(a2_init)) >= 1
        b0_init,b1_init,b2_init,a1_init,a2_init = (b0_init,b1_init,b2_init,a1_init,a2_init)./2
        shifta_init -= 1
        shiftb_init -= 1
    end

    while max(abs(b0_init),abs(b1_init),abs(b2_init)) >= 1
        b0_init,b1_init,b2_init = (b0_init,b1_init,b2_init)./2
        shiftb_init -= 1
    end
    while max(abs(a1_init),abs(a2_init)) >= 1
        a1_init,a2_init = (a1_init,a2_init)./2
        shifta_init -= 1
    end
    #println("b0_init: $b0_init, b1_init: $b1_init, b2_init: $b2_init, a1_init: $a1_init, a2_init: $a2_init")

    for wordlength in 4:1:10
        shifta = shifta_init + wordlength-1
        shiftb = shiftb_init + wordlength-1
        b0 = float_to_fxp(b0_init, wordlength)
        b1 = float_to_fxp(b1_init, wordlength)
        b2 = float_to_fxp(b2_init, wordlength)
        a1 = float_to_fxp(a1_init, wordlength)
        a2 = float_to_fxp(a2_init, wordlength)

        if shifta >= 0
            if !is_stable(a1/(1<<shifta), a2/(1<<shifta))
                continue
            end
        else
            if !is_stable(Float64(a1*(1<<-shifta)), Float64(a2*(1<<-shifta)))
                continue
            end
        end

        if shifta >= 0
            a1_eps, a2_eps = a1/(1<<shifta), a2/(1<<shifta)
        else
            a1_eps, a2_eps = a1*(1<<-shifta), a2*(1<<-shifta)
        end
        if shiftb >= 0
            b0_eps, b1_eps, b2_eps = b0/(1<<shiftb), b1/(1<<shiftb), b2/(1<<shiftb)
        else
            b0_eps, b1_eps, b2_eps = b0*(1<<-shiftb), b1*(1<<-shiftb), b2*(1<<-shiftb)
        end
        if shifta_init >= 0
            a1_init_eps, a2_init_eps = a1/(1<<shifta_init), a2/(1<<shifta_init)
        else
            a1_init_eps, a2_init_eps = a1*(1<<-shifta_init), a2*(1<<-shifta_init)
        end
        if shiftb_init >= 0
            b0_init_eps, b1_init_eps, b2_init_eps = b0/(1<<shiftb_init), b1/(1<<shiftb_init), b2/(1<<shiftb_init)
        else
            b0_init_eps, b1_init_eps, b2_init_eps = b0*(1<<-shiftb_init), b1*(1<<-shiftb_init), b2*(1<<-shiftb_init)
        end
        varepsilon = maxdiff(b0_init_eps, b1_init_eps, b2_init_eps,a1_init_eps, a2_init_eps,b0_eps, b1_eps, b2_eps, a1_eps, a2_eps)

        if !smallestwordlength_found && varepsilon < 1e-7
            smallestwordlength_found = true
            b0_saved, b1_saved, b2_saved, a1_saved, a2_saved = (b0, b1, b2, a1, a2)
            shiftb_saved, shifta_saved = (shiftb, shifta)
            wordlength_saved = wordlength
            println("MCM for $instancename -- $wordlength --- varepsilon: $varepsilon")
            @time begin
                model = Model(CPLEX.Optimizer)
                set_silent(model)
                solution_b = mcm(model, [b0, b1, b2], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true, avoid_internal_shifts = true)
                model = Model(CPLEX.Optimizer)
                set_silent(model)
                solution_a = mcm(model, [a1, a2], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true, avoid_internal_shifts = true)
            end # @time begin
            if instancename != "compilation"
                write_files && write_flopoco_solution([(a1, a2), (b0, b1, b2), solution_a, solution_b], instancename, wordlength, shifta, shiftb)
                write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd", sum(coef != 0 ? 1 : 0 for coef in [a1, a2, b0, b1, b2])-1+length(solution_b)+length(solution_a), varepsilon, wordlength)
                write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain", 0, varepsilon, wordlength)
                write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain0dsp", 0, varepsilon, wordlength)
            end
        end
    end

    if !smallestwordlength_found
        for wordlength in 11:1:32
            shifta = shifta_init + wordlength-1
            shiftb = shiftb_init + wordlength-1
            b0 = float_to_fxp(b0_init, wordlength)
            b1 = float_to_fxp(b1_init, wordlength)
            b2 = float_to_fxp(b2_init, wordlength)
            a1 = float_to_fxp(a1_init, wordlength)
            a2 = float_to_fxp(a2_init, wordlength)

            if !is_stable(a1/(1<<shifta), a2/(1<<shifta))
                continue
            end

            if shifta >= 0
                a1_eps, a2_eps = a1/(1<<shifta), a2/(1<<shifta)
            else
                a1_eps, a2_eps = a1*(1<<-shifta), a2*(1<<-shifta)
            end
            if shiftb >= 0
                b0_eps, b1_eps, b2_eps = b0/(1<<shiftb), b1/(1<<shiftb), b2/(1<<shiftb)
            else
                b0_eps, b1_eps, b2_eps = b0*(1<<-shiftb), b1*(1<<-shiftb), b2*(1<<-shiftb)
            end
            if shifta_init >= 0
                a1_init_eps, a2_init_eps = a1/(1<<shifta_init), a2/(1<<shifta_init)
            else
                a1_init_eps, a2_init_eps = a1*(1<<-shifta_init), a2*(1<<-shifta_init)
            end
            if shiftb_init >= 0
                b0_init_eps, b1_init_eps, b2_init_eps = b0/(1<<shiftb_init), b1/(1<<shiftb_init), b2/(1<<shiftb_init)
            else
                b0_init_eps, b1_init_eps, b2_init_eps = b0*(1<<-shiftb_init), b1*(1<<-shiftb_init), b2*(1<<-shiftb_init)
            end
            varepsilon = maxdiff(b0_init_eps, b1_init_eps, b2_init_eps,a1_init_eps, a2_init_eps,b0_eps, b1_eps, b2_eps, a1_eps, a2_eps)

            println("\t$instancename -- $wordlength --- varepsilon: $varepsilon")
            if !smallestwordlength_found && varepsilon < 1e-7
                smallestwordlength_found = true
                b0_saved, b1_saved, b2_saved, a1_saved, a2_saved = (b0, b1, b2, a1, a2)
                shiftb_saved, shifta_saved = (shiftb, shifta)
                wordlength_saved = wordlength
                if instancename != "compilation"
                    solution_b = rpag([b0, b1, b2])
                    solution_a = rpag([a1, a2])
                    write_files && write_flopoco_solution([(a1, a2), (b0, b1, b2), solution_a, solution_b], instancename, wordlength, shifta, shiftb)
                    write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd", sum(coef != 0 ? 1 : 0 for coef in [a1, a2, b0, b1, b2])-1+length(solution_b)+length(solution_a), varepsilon, wordlength)
                    write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain", 0, varepsilon, wordlength)
                    write_files && write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain0dsp", 0, varepsilon, wordlength)
                end
            end
        end
    end

    return b0_saved, b1_saved, b2_saved, a1_saved, a2_saved, shiftb_saved, shifta_saved, wordlength_saved
end

function truncate_hp0()
    verbose = true
    instancename, passband, stopband, deltapass, deltastop = ("hp0", (0.041, 1.0), (0.0, 0.0), 0.016, 0.0)
    fbands = [passband, stopband]
    abands = [(1-deltapass, 1+deltapass), (0.0, deltastop)]
    specifications = get_specifications(fbands, abands)

    design_method = "ellip"
    match_exactly = "both"
    eval_string("""
        Fstop = 1e-09;          % Stopband Frequency
        Fpass = 0.041;          % Passband Frequency
        Astop = 290;            % Stopband Attenuation (dB)
        Apass = 0.14009803137;  % Passband Ripple (dB)
        $instancename = designfilt('highpassiir','DesignMethod','$design_method','StopbandFrequency',Fstop,'PassbandFrequency',Fpass,'MatchExactly','$match_exactly','StopbandAttenuation',Astop,'PassbandRipple',Apass);
        coefs_str = num2hex($instancename.Coefficients);
        nb_sections = size($instancename.Coefficients, 1);"""
    )

    if @mget(nb_sections) != 1
        println("No second order IIR filter using $design_method ($match_exactly) for $instancename")
    end

    eval_string("coef = coefs_str(1,:);")
    b0_hex = @mget(coef)
    b0 = reinterpret(Float64, reverse(hex2bytes(b0_hex)))[1]
    eval_string("coef = coefs_str(2,:);")
    b1_hex = @mget(coef)
    b1 = reinterpret(Float64, reverse(hex2bytes(b1_hex)))[1]
    eval_string("coef = coefs_str(3,:);")
    b2_hex = @mget(coef)
    b2 = reinterpret(Float64, reverse(hex2bytes(b2_hex)))[1]
    eval_string("coef = coefs_str(5,:);")
    a1_hex = @mget(coef)
    a1 = reinterpret(Float64, reverse(hex2bytes(a1_hex)))[1]
    eval_string("coef = coefs_str(6,:);")
    a2_hex = @mget(coef)
    a2 = reinterpret(Float64, reverse(hex2bytes(a2_hex)))[1]

    verbose && println("coefficients_b: $(b0), $(b1), $(b2)")
    verbose && println("coefficients_a: $(a1), $(a2)")

    # Check that they fit
    varepsilon = coefs_fit_specs_varepsilon(b0,b1,b2,a1,a2,specifications)
    if varepsilon != 0.0
        println("Matlab error for $instancename, varepsilon = $varepsilon")
    end
    write_table_epsilon("$(instancename)_dw00_cw00_fixiir", 0, varepsilon, 0)
    write_flopoco_solution(b0_hex, b1_hex, b2_hex, a1_hex, a2_hex, instancename)

    b0,b1,b2,a1,a2,shiftb,shifta,wordlength = truncate(b0,b1,b2,a1,a2,specifications,instancename, write_files=true)
    verbose && println("wordlength: $(wordlength)")
    verbose && println("coefficients_b: $(b0), $(b1), $(b2)")
    verbose && println("coefficients_a: $(a1), $(a2)")
    verbose && println("shiftb: $shiftb")
    verbose && println("shifta: $shifta")
    verbose && println("eps: $(varepsilon)")
end


function truncate_hp0()
    verbose = true
    instancename = "hp0"
    write_table_epsilon("$(instancename)_dw00_cw00_fixiir", 0, 0.0, 0)
    b0, b1, b2 = 101.8, -203.4, 101.6
    a1, a2 = -1.967, 0.968
    write_flopoco_solution(b0, b1, b2, a1, a2, instancename)

    b0,b1,b2,a1,a2,shiftb,shifta,wordlength = truncate(b0,b1,b2,a1,a2,instancename, write_files=true)
    verbose && println("wordlength: $(wordlength)")
    verbose && println("coefficients_b: $(b0), $(b1), $(b2)")
    verbose && println("coefficients_a: $(a1), $(a2)")
    verbose && println("shiftb: $shiftb")
    verbose && println("shifta: $shifta")
    #verbose && println("eps: $(varepsilon)")
end




function main()
    initialize_files()
    main_hp0()
    truncate_hp0()
    benchmarking()
    #special_cases()
    read_and_rewrite()
    write_special_cases()
    return nothing
end
