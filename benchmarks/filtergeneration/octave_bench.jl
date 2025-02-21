using JuMP, jIIR2AG, CPLEX

function generate_passband(Fc, Fs, N)
    f1 = Fc/(2^(1/6))
    f2 = Fc*(2^(1/6))
    Qr = Fc/(f2-f1)
    Qd = Qr*(pi/(2*N))/sin(pi/(2*N))
    alpha = (1+sqrt(1+4*Qd^2))/(2*Qd)
    W1 = 2*Fc/(Fs*alpha)
    W2 = 2*Fc*alpha/Fs
    return (W1, W2)
end

function generate_passband(Fc, Fs, N)
    W1 = 2*(Fc/sqrt(2))/Fs
    W2 = 2*(Fc*sqrt(2))/Fs
    return (W1, W2)
end


function generate_specifications(Fc, Fs, N)
    specifications = Vector{Tuple{Float64, Float64, Float64}}()
    w1, w2 = generate_passband(Fc, Fs, N)
    println("w1: $w1 -- w2: $w2")
    dbatten = 10^(-3/20)
    onepercent = (w2-w1)/100
    for i in 0.0:onepercent:(w1-onepercent)
        push!(specifications, (i, 0.0, dbatten))
    end
    for i in (w1+onepercent):onepercent:(w2-onepercent)
        push!(specifications, (i, dbatten, 1.0))
    end
    for i in (w2+onepercent):onepercent:1.0
        push!(specifications, (i, 0.0, dbatten))
    end
    push!(specifications, (1.0, 0.0, dbatten))
    println(length(specifications))
    return specifications
end

wordlength = 6
model = Model(CPLEX.Optimizer)
verbose = true
design_second_order_iir(model, generate_specifications(1000, 44100, 2), wordlength, avoid_internal_shifts=true, use_big_m=true, presolve_time_sec=600, verbose=verbose)
