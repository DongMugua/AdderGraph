function wcpg(solution::Iir2AdderGraphs)::Float64
    W = zeros(1)
    wcpg_success = 0
    if !isempty(Libc.find_library("libwcpg"))
        b = collect(solution.coefficients_b)./(2^solution.shifts[2])
        a = collect(solution.coefficients_a)./(2^solution.shifts[1])
        wcpg_success = ccall((:WCPG_tf, "libwcpg"), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cuintmax_t, Cuintmax_t, Cint), W, collect(b), collect(a), length(b), length(a), 0)
    else
        @warn "libwcpg not found"
    end
    if wcpg_success == 0
        @warn "WCPG computation errored"
        a1, a2 = solution.coefficients_a
        b0, b1, b2 = solution.coefficients_b
        shifta, shiftb = solution.shifts
        C = Matrix{BigInt}([b1*2^(shifta)-b0*a1 b2*2^(shifta)-b0*a2])
        A = Matrix{BigInt}([-a1 -a2; 2^(shifta) 0])
        B = Vector{BigInt}([1, 0])
        D = Vector{BigInt}([b0])
        N = 1001
        W[1] = (BigInt(2)^(shiftb)*sum((abs.(C*(A^k)*B))[1]*(BigInt(2)^(shifta))^(N-k) for k in 0:N)+abs(D[1])*(BigInt(2)^(shifta)*BigInt(2)^(shiftb))*(BigInt(2)^(shifta))^N)/(BigInt(2)^(shiftb)*BigInt(2)^(shiftb)*BigInt(2)^(shifta)*(BigInt(2)^(shifta))^N)
    end
    return W[1]
end


function wcpg_eps(solution::Iir2AdderGraphs)::Float64
    W = zeros(1)
    wcpg_success = 0
    if !isempty(Libc.find_library("libwcpg"))
        b = zeros(length(solution.coefficients_b))
        b[1] = 1.0
        a = collect(solution.coefficients_a)./(2^solution.shifts[1])
        wcpg_success = ccall((:WCPG_tf, "libwcpg"), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cuintmax_t, Cuintmax_t, Cint), W, collect(b), collect(a), length(b), length(a), 0)
    else
        @warn "libwcpg not found"
    end
    if wcpg_success == 0
        @warn "WCPG_eps computation errored"
        a1, a2 = solution.coefficients_a
        b0, b1, b2 = 1,0,0
        shifta, shiftb = solution.shifts[1], 0
        C = Matrix{BigInt}([b1*2^(shifta)-b0*a1 b2*2^(shifta)-b0*a2])
        A = Matrix{BigInt}([-a1 -a2; 2^(shifta) 0])
        B = Vector{BigInt}([1, 0])
        D = Vector{BigInt}([b0])
        N = 1001
        W[1] = (BigInt(2)^(shiftb)*sum((abs.(C*(A^k)*B))[1]*(BigInt(2)^(shifta))^(N-k) for k in 0:N)+abs(D[1])*(BigInt(2)^(shifta)*BigInt(2)^(shiftb))*(BigInt(2)^(shifta))^N)/(BigInt(2)^(shiftb)*BigInt(2)^(shiftb)*BigInt(2)^(shifta)*(BigInt(2)^(shifta))^N)
    end
    return W[1]
end
