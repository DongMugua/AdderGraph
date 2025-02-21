"""
    get_gb_secondorder_iir(model::Model,
                           specifications::Vector{Tuple{Float64, Float64, Float64}})

Gives the global scaling factor for `b` according to the `specifications`.
"""
function get_gb_secondorder_iir(model::Model,
                                specifications::Vector{Tuple{Float64, Float64, Float64}},
                                ;verbose::Bool=false)::Int
    init_silent = get_optimizer_attribute(model, MOI.Silent())
    empty!(model)
    if !init_silent
        set_silent(model)
    end
    @variable(model, b0)
    @variable(model, b1)
    @variable(model, b2)
    @constraint(model, [i in 1:length(specifications)], b0^2+b1^2+b2^2+2*cos(2*pi*specifications[i][1])*b0*b2+2*cos(pi*specifications[i][1])*b1*b2+2*cos(pi*specifications[i][1])*b0*b1 <= specifications[i][3]^2*16)
    @objective(model, Min, b0)
    optimize!(model)
    if !has_values(model)
        model[:gb_obtained] = false
        return 0
    end
    b0 = -value(b0)
    verbose && println("\tMax b0: $b0")
    @objective(model, Min, b1)
    optimize!(model)
    b1 = -value(b1)
    verbose && println("\tMax b1: $b1")
    verbose && println("\tMax b2: $b0")
    gb = Int(ceil(log2(max(b0, b1))))
    if !init_silent
        unset_silent(model)
    end
    empty!(model)
    model[:gb_obtained] = true
    return gb
end



"""
    iir_design!(model::Model,
                specifications::Vector{Tuple{Float64, Float64, Float64}},
                wordlength::Int,
                ;verbose::Bool = false)

Modify `model` to add iir specific constraints according to `specifications`.
"""
function iir_design!(model::Model,
                     specifications::Vector{Tuple{Float64, Float64, Float64}},
                     wordlength::Int,
                     ;with_symmetry_breaking::Bool=true,
                     use_big_m::Bool=true,
                     verbose::Bool=false,
                     maximum_g::Int=-1)

    minimum_value = -1 << (wordlength-1)
    maximum_value = 1 << (wordlength-1)-1
    minimum_value_coeffs = minimum_value
    maximum_value_coeffs = maximum_value
    if maximum_g < 0
        maximum_g = wordlength
    end
    init_ga, init_gb = model[:scaling]

    max_abs_value_coeffs = max(abs(minimum_value_coeffs), abs(maximum_value_coeffs))
    minimum_value_quad = minimum_value_coeffs*maximum_value_coeffs
    maximum_value_quad = max_abs_value_coeffs*max_abs_value_coeffs
    maximum_value_g = 1 << maximum_g
    maximum_value_gsq = maximum_value_g^2
    minimum_value_quad_g = minimum_value_coeffs*maximum_value_g
    maximum_value_quad_g = maximum_value_coeffs*maximum_value_g

    verbose && use_big_m && println("Use big M for filter design")
    verbose && !use_big_m && println("Use indicator constraints for filter design")
    verbose && with_symmetry_breaking && println("Use symmetry breaking constraint")

    @variable(model, minimum_value_coeffs <= a[1:2] <= maximum_value_coeffs, Int)
    @variable(model, 0 <= apos[1:2] <= max_abs_value_coeffs, Int)
    @variable(model, asign[1:2], Bin)
    @variable(model, ta[1:2, 0:(wordlength-1)], Bin)
    @variable(model, 0 <= asquare[1:2] <= maximum_value_quad, Int)
    @variable(model, minimum_value_quad <= aquad[m in 1:(2-1), (m+1):2] <= maximum_value_quad, Int)
    @variable(model, 0 <= uaquadsq[m in 1:2, m:2, 0:(wordlength-1)] <= max_abs_value_coeffs, Int)

    if use_big_m
        @constraint(model, [m in 1:2], a[m] >= minimum_value_coeffs*(1-asign[m]))
        @constraint(model, [m in 1:2], a[m] <= maximum_value_coeffs*asign[m])
    else
        @constraint(model, [m in 1:2], asign[m] => {a[m] >= 0})
        @constraint(model, [m in 1:2], !asign[m] => {a[m] <= 0})
    end
    if use_big_m
        @constraint(model, [m in 1:2], apos[m] >= a[m])
        @constraint(model, [m in 1:2], apos[m] <= a[m] + (maximum_value_coeffs-minimum_value_coeffs)*(1-asign[m]))
        @constraint(model, [m in 1:2], apos[m] >= -a[m])
        @constraint(model, [m in 1:2], apos[m] <= -a[m] + (maximum_value_coeffs-minimum_value_coeffs)*asign[m])
    else
        @constraint(model, [m in 1:2], asign[m] => {apos[m] == a[m]})
        @constraint(model, [m in 1:2], !asign[m] => {apos[m] == -a[m]})
    end
    @constraint(model, [m in 1:2], apos[m] == sum(2^j*ta[m,j] for j in 0:(wordlength-1)))

    if use_big_m
        @constraint(model, [m in 1:2, n in m:2, j in 0:(wordlength-1)], uaquadsq[m,n,j] <= max_abs_value_coeffs*ta[n,j])
        @constraint(model, [m in 1:2, n in m:2, j in 0:(wordlength-1)], uaquadsq[m,n,j] <= apos[m])
        @constraint(model, [m in 1:2, n in m:2, j in 0:(wordlength-1)], uaquadsq[m,n,j] >= apos[m] - max_abs_value_coeffs*(1-ta[n,j]))
    else
        @constraint(model, [m in 1:2, n in m:2, j in 0:(wordlength-1)], ta[n,j] => {uaquadsq[m,n,j] == apos[m]})
        @constraint(model, [m in 1:2, n in m:2, j in 0:(wordlength-1)], !ta[n,j] => {uaquadsq[m,n,j] == 0})
    end
    @constraint(model, [m in 1:2], asquare[m] == sum(2^j*uaquadsq[m,m,j] for j in 0:(wordlength-1)))

    if use_big_m
        @constraint(model, [m in 1:(2-1), n in (m+1):2], aquad[m,n] >= sum(2^j*uaquadsq[m,n,j] for j in 0:(wordlength-1)) - (maximum_value_quad-minimum_value_quad)*(asign[m]+asign[n]))
        @constraint(model, [m in 1:(2-1), n in (m+1):2], aquad[m,n] >= sum(2^j*uaquadsq[m,n,j] for j in 0:(wordlength-1)) - (maximum_value_quad-minimum_value_quad)*(2-asign[m]-asign[n]))
        @constraint(model, [m in 1:(2-1), n in (m+1):2], aquad[m,n] <= sum(2^j*uaquadsq[m,n,j] for j in 0:(wordlength-1)))
        @constraint(model, [m in 1:(2-1), n in (m+1):2], aquad[m,n] >= -sum(2^j*uaquadsq[m,n,j] for j in 0:(wordlength-1)))
        @constraint(model, [m in 1:(2-1), n in (m+1):2], aquad[m,n] <= -sum(2^j*uaquadsq[m,n,j] for j in 0:(wordlength-1)) + (maximum_value_quad-minimum_value_quad)*(1-asign[m]+asign[n]))
        @constraint(model, [m in 1:(2-1), n in (m+1):2], aquad[m,n] <= -sum(2^j*uaquadsq[m,n,j] for j in 0:(wordlength-1)) + (maximum_value_quad-minimum_value_quad)*(1+asign[m]-asign[n]))
    else
        @variable(model, samesigna[m in 1:2, m:2], Bin)
        @constraint(model, [m in 1:(2-1), n in (m+1):2], samesigna[m,n] => {asign[m] == asign[n]})
        @constraint(model, [m in 1:(2-1), n in (m+1):2], !samesigna[m,n] => {asign[m] == 1-asign[n]})
        @constraint(model, [m in 1:(2-1), n in (m+1):2], samesigna[m,n] => {aquad[m,n] == sum(2^j*uaquadsq[m,n,j] for j in 0:(wordlength-1))})
        @constraint(model, [m in 1:(2-1), n in (m+1):2], !samesigna[m,n] => {aquad[m,n] == -sum(2^j*uaquadsq[m,n,j] for j in 0:(wordlength-1))})
    end

    @variable(model, minimum_value_coeffs <= b[0:2] <= maximum_value_coeffs, Int)
    @variable(model, 0 <= bpos[0:2] <= max_abs_value_coeffs, Int)
    @variable(model, bsign[0:2], Bin)
    @variable(model, tb[0:2, 0:(wordlength-1)], Bin)
    @variable(model, 0 <= bsquare[0:2] <= maximum_value_quad, Int)
    @variable(model, minimum_value_quad <= bquad[m in 0:(2-1), (m+1):2] <= maximum_value_quad, Int)
    @variable(model, 0 <= ubquadsq[m in 0:2, m:2, 0:(wordlength-1)] <= max_abs_value_coeffs, Int)

    if use_big_m
        @constraint(model, [m in 0:2], b[m] >= minimum_value_coeffs*(1-bsign[m]))
        @constraint(model, [m in 0:2], b[m] <= maximum_value_coeffs*bsign[m])
    else
        @constraint(model, [m in 0:2], bsign[m] => {b[m] >= 0})
        @constraint(model, [m in 0:2], !bsign[m] => {b[m] <= 0})
    end
    if use_big_m
        @constraint(model, [m in 0:2], bpos[m] >= b[m])
        @constraint(model, [m in 0:2], bpos[m] <= b[m] + (maximum_value_coeffs-minimum_value_coeffs)*(1-bsign[m]))
        @constraint(model, [m in 0:2], bpos[m] >= -b[m])
        @constraint(model, [m in 0:2], bpos[m] <= -b[m] + (maximum_value_coeffs-minimum_value_coeffs)*bsign[m])
    else
        @constraint(model, [m in 0:2], bsign[m] => {bpos[m] == b[m]})
        @constraint(model, [m in 0:2], !bsign[m] => {bpos[m] == -b[m]})
    end
    @constraint(model, [m in 0:2], bpos[m] == sum(2^j*tb[m,j] for j in 0:(wordlength-1)))

    if use_big_m
        @constraint(model, [m in 0:2, n in m:2, j in 0:(wordlength-1)], ubquadsq[m,n,j] <= max_abs_value_coeffs*tb[n,j])
        @constraint(model, [m in 0:2, n in m:2, j in 0:(wordlength-1)], ubquadsq[m,n,j] <= bpos[m])
        @constraint(model, [m in 0:2, n in m:2, j in 0:(wordlength-1)], ubquadsq[m,n,j] >= bpos[m] - max_abs_value_coeffs*(1-tb[n,j]))
    else
        @constraint(model, [m in 0:2, n in m:2, j in 0:(wordlength-1)], tb[n,j] => {ubquadsq[m,n,j] == bpos[m]})
        @constraint(model, [m in 0:2, n in m:2, j in 0:(wordlength-1)], !tb[n,j] => {ubquadsq[m,n,j] == 0})
    end
    @constraint(model, [m in 0:2], bsquare[m] == sum(2^j*ubquadsq[m, m, j] for j in 0:(wordlength-1)))

    if use_big_m
        @constraint(model, [m in 0:(2-1), n in (m+1):2], bquad[m,n] >= sum(2^j*ubquadsq[m,n,j] for j in 0:(wordlength-1)) - (maximum_value_quad-minimum_value_quad)*(bsign[m]+bsign[n]))
        @constraint(model, [m in 0:(2-1), n in (m+1):2], bquad[m,n] >= sum(2^j*ubquadsq[m,n,j] for j in 0:(wordlength-1)) - (maximum_value_quad-minimum_value_quad)*(2-bsign[m]-bsign[n]))
        @constraint(model, [m in 0:(2-1), n in (m+1):2], bquad[m,n] <= sum(2^j*ubquadsq[m,n,j] for j in 0:(wordlength-1)))
        @constraint(model, [m in 0:(2-1), n in (m+1):2], bquad[m,n] >= -sum(2^j*ubquadsq[m,n,j] for j in 0:(wordlength-1)))
        @constraint(model, [m in 0:(2-1), n in (m+1):2], bquad[m,n] <= -sum(2^j*ubquadsq[m,n,j] for j in 0:(wordlength-1)) + (maximum_value_quad-minimum_value_quad)*(1-bsign[m]+bsign[n]))
        @constraint(model, [m in 0:(2-1), n in (m+1):2], bquad[m,n] <= -sum(2^j*ubquadsq[m,n,j] for j in 0:(wordlength-1)) + (maximum_value_quad-minimum_value_quad)*(1+bsign[m]-bsign[n]))
    else
        @variable(model, samesignb[m in 0:2, m:2], Bin)
        @constraint(model, [m in 0:(2-1), n in (m+1):2], samesignb[m,n] => {bsign[m] == bsign[n]})
        @constraint(model, [m in 0:(2-1), n in (m+1):2], !samesignb[m,n] => {bsign[m] == 1-bsign[n]})
        @constraint(model, [m in 0:(2-1), n in (m+1):2], samesignb[m,n] => {bquad[m,n] == sum(2^j*ubquadsq[m,n,j] for j in 0:(wordlength-1))})
        @constraint(model, [m in 0:(2-1), n in (m+1):2], !samesignb[m,n] => {bquad[m,n] == -sum(2^j*ubquadsq[m,n,j] for j in 0:(wordlength-1))})
    end

    @variable(model, a0b[0:maximum_g], Bin)
    @constraint(model, sum(a0b[j] for j in 0:maximum_g) == 1)
    @variable(model, 1 <= a0 <= maximum_value_g, Int)
    @constraint(model, a0 == sum(2^j * a0b[j] for j in 0:maximum_g))
    @variable(model, 1 <= a0square <= maximum_value_gsq, Int)
    @constraint(model, a0square == sum(2^(2j) * a0b[j] for j in 0:maximum_g))
    @variable(model, minimum_value_quad_g <= a01 <= maximum_value_quad_g, Int)
    @variable(model, minimum_value_quad_g <= a02 <= maximum_value_quad_g, Int)
    @variable(model, 0 <= ua0quad[m in 1:2, 0:(wordlength-1)] <= max_abs_value_coeffs, Int)

    if use_big_m
        @constraint(model, [m in 1:2, j in 0:(wordlength-1)], ua0quad[m,j] <= max_abs_value_coeffs*a0b[j])
        @constraint(model, [m in 1:2, j in 0:(wordlength-1)], ua0quad[m,j] <= apos[m])
        @constraint(model, [m in 1:2, j in 0:(wordlength-1)], ua0quad[m,j] >= apos[m] - max_abs_value_coeffs*(1-a0b[j]))
    else
        @constraint(model, [m in 1:2, j in 0:(wordlength-1)], a0b[j] => {ua0quad[m,j] == apos[m]})
        @constraint(model, [m in 1:2, j in 0:(wordlength-1)], !a0b[j] => {ua0quad[m,j] == 0})
    end

    if use_big_m
        @constraint(model, a01 <= sum(2^j*ua0quad[1,j] for j in 0:(wordlength-1)) + (maximum_value_quad_g)*(1-asign[1]))
        @constraint(model, a01 >= sum(2^j*ua0quad[1,j] for j in 0:(wordlength-1)) - (maximum_value_quad_g-minimum_value_quad_g)*(1-asign[1]))
        @constraint(model, a01 <= -sum(2^j*ua0quad[1,j] for j in 0:(wordlength-1)) + (maximum_value_quad_g-minimum_value_quad_g)*(asign[1]))
        @constraint(model, a01 >= -sum(2^j*ua0quad[1,j] for j in 0:(wordlength-1)) - (-minimum_value_quad_g)*(asign[1]))
        @constraint(model, a02 <= sum(2^j*ua0quad[2,j] for j in 0:(wordlength-1)) + (maximum_value_quad_g)*(1-asign[2]))
        @constraint(model, a02 >= sum(2^j*ua0quad[2,j] for j in 0:(wordlength-1)) - (maximum_value_quad_g-minimum_value_quad_g)*(1-asign[2]))
        @constraint(model, a02 <= -sum(2^j*ua0quad[2,j] for j in 0:(wordlength-1)) + (maximum_value_quad_g-minimum_value_quad_g)*(asign[2]))
        @constraint(model, a02 >= -sum(2^j*ua0quad[2,j] for j in 0:(wordlength-1)) - (-minimum_value_quad_g)*(asign[2]))
    else
        @constraint(model, asign[1] => {a01 == sum(2^j*ua0quad[1,j] for j in 0:(wordlength-1))})
        @constraint(model, !asign[1] => {a01 == -sum(2^j*ua0quad[1,j] for j in 0:(wordlength-1))})
        @constraint(model, asign[2] => {a02 == sum(2^j*ua0quad[2,j] for j in 0:(wordlength-1))})
        @constraint(model, !asign[2] => {a02 == -sum(2^j*ua0quad[2,j] for j in 0:(wordlength-1))})
    end

    # Frequency
    @constraint(model, [i in 1:length(specifications)],
        specifications[i][2]^2*((2^(2*wordlength-2-2*init_ga))*a0square+asquare[1]+asquare[2]
            +2*cos(2*pi*specifications[i][1])*a02*(2^(wordlength-1-init_ga))+2*cos(pi*specifications[i][1])*aquad[1,2]+2*cos(pi*specifications[i][1])*a01*(2^(wordlength-1-init_ga)))
            <= (bsquare[0]+bsquare[1]+bsquare[2]
                +2*cos(2*pi*specifications[i][1])*bquad[0,2]+2*cos(pi*specifications[i][1])*bquad[1,2]+2*cos(pi*specifications[i][1])*bquad[0,1]))
    @constraint(model, [i in 1:length(specifications)],
        (bsquare[0]+bsquare[1]+bsquare[2]
            +2*cos(2*pi*specifications[i][1])*bquad[0,2]+2*cos(pi*specifications[i][1])*bquad[1,2]+2*cos(pi*specifications[i][1])*bquad[0,1])
            <= specifications[i][3]^2*((2^(2*wordlength-2-2*init_ga))*a0square+asquare[1]+asquare[2]
                +2*cos(2*pi*specifications[i][1])*a02*(2^(wordlength-1-init_ga))+2*cos(pi*specifications[i][1])*aquad[1,2]+2*cos(pi*specifications[i][1])*a01*(2^(wordlength-1-init_ga))))

    # Symmetry breaking
    if with_symmetry_breaking
        @constraint(model, b[0] >= bpos[2])
    end

    # Poles
    @constraint(model, a[2] >= apos[1] - a0*2^(wordlength-1-init_ga) + 1)
    @constraint(model, a[2] <= a0*2^(wordlength-1-init_ga) - 1)
    @constraint(model, apos[1] <= 2*a0*2^(wordlength-1-init_ga) - 1)

    # For link with MCM:
    model[:ToLink] = ([apos[1], apos[2]], [bpos[0], bpos[1], bpos[2]])

    return model
end
