using JuMP
using CPLEX
using MATLAB
# include("../src/AdderGraphs/src/AdderGraphs.jl")
# using .AdderGraphs
# include("../src/jMCM/src/jMCM.jl")
# using .jMCM
# include("../src/jIIR2AG/src/jIIR2AG.jl")
# using .jIIR2AG
include("../src/includepkg.jl")


function using_matlab()
    lp1x4_5 = ("lp1x4", 5, (0.0, 0.3), (0.7, 1.0), 0.06, 0.06)
    instancename, wordlength, passband, stopband, deltapass, deltastop = lp1x4_5
    design_method = "ellip"
    match_exactly = "both"
    eval_string("$instancename = designfilt('lowpassiir','DesignMethod','$design_method','PassbandFrequency',$(passband[2]),'StopbandFrequency',$(stopband[1]),'MatchExactly','$match_exactly','PassbandRipple',$(20*log10(1/(1.0-deltapass))),'StopbandAttenuation',$(20*log10(1/deltastop)));
        coefs_str = num2hex($instancename.Coefficients);
        nb_sections = size($instancename.Coefficients, 1);")

    if @mget(nb_sections) != 1
        println("No second order IIR filter using $design_method ($match_exactly) for $instancename")
        return Iir2AdderGraphs()
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

    # Truncate coefficients
    b0_init,b1_init,b2_init,a1_init,a2_init = b0,b1,b2,a1,a2

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
    shifta = shifta_init + wordlength-1
    shiftb = shiftb_init + wordlength-1
    b0_truncated = round(Int, b0_init*2.0^(wordlength-1))
    b1_truncated = round(Int, b1_init*2.0^(wordlength-1))
    b2_truncated = round(Int, b2_init*2.0^(wordlength-1))
    a1_truncated = round(Int, a1_init*2.0^(wordlength-1))
    a2_truncated = round(Int, a2_init*2.0^(wordlength-1))

    model = Model(CPLEX.Optimizer)
    set_silent(model)
    solution_b = mcm(model, [b0_truncated, b1_truncated, b2_truncated], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true, avoid_internal_shifts = true)
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    solution_a = mcm(model, [a1_truncated, a2_truncated], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true, avoid_internal_shifts = true)

    return Iir2AdderGraphs((a1_truncated, a2_truncated), (b0_truncated, b1_truncated, b2_truncated), (shifta, shiftb), solution_a, solution_b)
end



function our_method()
    lp1x4_5 = ("lp1x4", 5, (0.0, 0.3), (0.7, 1.0), 0.06, 0.06)
    instancename, wordlength, passband, stopband, deltapass, deltastop = lp1x4_5
    fbands = [passband, stopband]
    abands = [(1-deltapass, 1+deltapass), (0.0, deltastop)]
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    solution = design_second_order_iir(model, fbands, abands, wordlength, avoid_internal_shifts=true)
    return solution
end
