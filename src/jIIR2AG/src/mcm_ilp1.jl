function mcm_ilp_odd!(model::Model, C::Tuple{Vector{VariableRef}, Vector{VariableRef}},
                      wordlength::Int, NA::Int;
                      avoid_internal_shifts::Bool=false,
                      only_positive::Bool=false,
                      no_output_shifts::Bool=false)
    Smin = -wordlength
    Smax = wordlength
    maximum_value = 1 << (wordlength-1)

    @variable(model, 1 <= ca[0:NA] <= maximum_value, Int)

    if NA > 0
        @variable(model, 1 <= ca_no_shift[1:NA] <= maximum_value, Int)
        @variable(model, 1 <= cai[1:NA, 1:2] <= maximum_value, Int)
        @variable(model, 1 <= ca_right_sh[1:NA] <= maximum_value, Int)
        @variable(model, -maximum_value <= cai_right_shsg[1:NA] <= maximum_value, Int)
        @variable(model, -maximum_value <= cai_left_sg[1:NA] <= maximum_value, Int)

        @variable(model, Phiai[1:NA, 1:2], Bin)
        @variable(model, caik[a in 1:NA, 1:2, 0:(a-1)], Bin)
        @variable(model, phias[1:NA, 0:Smax], Bin)
    end

    @variable(model, oaasgsh[0:NA, 1:2, 0:Smax, 0:1], Bin)
    @variable(model, oabsgsh[0:NA, 0:2, 0:Smax, 0:1], Bin)
    if only_positive
        @constraint(model, [l in 0:NA, m in 1:2, s in 0:Smax], oaasgsh[l,m,s,1] == 0)
        @constraint(model, [l in 0:NA, m in 0:2, s in 0:Smax], oabsgsh[l,m,s,1] == 0)
    end
    if no_output_shifts
        @constraint(model, [l in 0:NA, m in 1:2, s in 1:Smax, phi in 0:1], oaasgsh[l,m,s,phi] == 0)
        @constraint(model, [l in 0:NA, m in 0:2, s in 1:Smax, phi in 0:1], oabsgsh[l,m,s,phi] == 0)
    end

    if NA > 0
        # Glue
        @variable(model, mcma[0:NA], Bin)
        @variable(model, mcmb[0:NA], Bin)
        @variable(model, samemcm[k in 0:NA, (k+1):NA], Bin)
        @variable(model, mcmbin[k in 1:NA, (k+1):NA], Bin)
        fix(mcma[0], 1, force=true)
        fix(mcmb[0], 1, force=true)
        for i in 1:NA
            fix(samemcm[0,i], 1, force=true)
        end
        @constraint(model, [i in 1:NA], mcma[i] + mcmb[i] == 1)
        @constraint(model, [i in 1:(NA-1)], mcma[i+1] <= mcma[i])
        @constraint(model, [k in 1:NA, l in (k+1):NA], 2*samemcm[k,l] >= mcma[k]+mcma[l]-1)
        @constraint(model, [k in 1:NA, l in (k+1):NA], 2*samemcm[k,l] >= mcmb[k]+mcmb[l]-1)
        @constraint(model, [k in 1:NA, l in (k+1):NA], mcma[k]+mcma[l]+1 == 2*mcmbin[k,l]+samemcm[k,l])
    end

    @variable(model, azero[m in 1:2], Bin)
    @variable(model, bzero[m in 0:2], Bin)

    maximum_value = 1 << (wordlength)

    @constraint(model, [m in 1:2], C[1][m] <= maximum_value*(1-azero[m]))
    @constraint(model, [m in 1:2], C[1][m] >= -maximum_value*(1-azero[m]))
    @constraint(model, [m in 0:2], C[2][m+1] <= maximum_value*(1-bzero[m]))
    @constraint(model, [m in 0:2], C[2][m+1] >= -maximum_value*(1-bzero[m]))

    # C1
    fix(ca[0], 1, force=true)

    if NA > 0
        # C2 - Modified
        @constraint(model, [l in 1:NA], ca_no_shift[l] == cai_left_sg[l] + cai_right_shsg[l])
        # C3a - C3b - add one constraint
        @constraint(model, [l in 1:NA, i in 1:2, k in 0:(l-1)], cai[l,i] <= ca[k] + (1-caik[l,i,k])*maximum_value)
        @constraint(model, [l in 1:NA, i in 1:2, k in 0:(l-1)], cai[l,i] >= ca[k] - (1-caik[l,i,k])*maximum_value)
        @constraint(model, [l in 1:NA, i in 1:2, k in 0:(l-1)], caik[l,i,k] <= samemcm[k,l])
        @constraint(model, [l in 1:NA, i in 1:2], sum(caik[l,i,k] for k in 0:(l-1)) == 1)
        # C4a - C4b - Modified
        @constraint(model, [l in 1:NA, s in 0:Smax], ca_right_sh[l] <= 2^s*cai[l,1] + (1-phias[l,s])*maximum_value)
        @constraint(model, [l in 1:NA, s in 0:Smax], ca_right_sh[l] >= 2^s*cai[l,1] - (1-phias[l,s])*(maximum_value*(2^s)))
        @constraint(model, [l in 1:NA], sum(phias[l,s] for s in 0:Smax) == 1)
        # C5a - C5b - C5c - Modified
        @constraint(model, [l in 1:NA], cai_right_shsg[l] <= ca_right_sh[l] + Phiai[l,1]*maximum_value)
        @constraint(model, [l in 1:NA], cai_right_shsg[l] >= ca_right_sh[l] - Phiai[l,1]*(2*maximum_value))
        @constraint(model, [l in 1:NA], cai_right_shsg[l] <= -ca_right_sh[l] + (1-Phiai[l,1])*(2*maximum_value))
        @constraint(model, [l in 1:NA], cai_right_shsg[l] >= -ca_right_sh[l] - (1-Phiai[l,1])*maximum_value)
        @constraint(model, [l in 1:NA], cai_left_sg[l] <= cai[l,2] + Phiai[l,2]*maximum_value)
        @constraint(model, [l in 1:NA], cai_left_sg[l] >= cai[l,2] - Phiai[l,2]*(2*maximum_value))
        @constraint(model, [l in 1:NA], cai_left_sg[l] <= -cai[l,2] + (1-Phiai[l,2])*(2*maximum_value))
        @constraint(model, [l in 1:NA], cai_left_sg[l] >= -cai[l,2] - (1-Phiai[l,2])*maximum_value)
        @constraint(model, [l in 1:NA], Phiai[l,1] + Phiai[l,2] <= 1)
    end

    # C6a - C6b - Modified
    @constraint(model, [l in 0:NA, m in 1:2, s in 0:Smax, phi in 0:1],
        2^s*ca[l] <= (1-2*phi)*C[1][m] + (1 - oaasgsh[l,m,s,phi])*maximum_value*2^s)
    @constraint(model, [l in 0:NA, m in 1:2, s in 0:Smax, phi in 0:1],
        2^s*ca[l] >= (1-2*phi)*C[1][m] - (1 - oaasgsh[l,m,s,phi])*maximum_value*2^s)
    @constraint(model, [l in 0:NA, m in 0:2, s in 0:Smax, phi in 0:1],
        2^s*ca[l] <= (1-2*phi)*C[2][m+1] + (1 - oabsgsh[l,m,s,phi])*maximum_value*2^s)
    @constraint(model, [l in 0:NA, m in 0:2, s in 0:Smax, phi in 0:1],
        2^s*ca[l] >= (1-2*phi)*C[2][m+1] - (1 - oabsgsh[l,m,s,phi])*maximum_value*2^s)
    if NA > 0
        @constraint(model, [l in 0:NA, m in 1:2, s in 0:Smax, phi in 0:1], oaasgsh[l,m,s,phi] <= mcma[l])
        @constraint(model, [l in 0:NA, m in 0:2, s in 0:Smax, phi in 0:1], oabsgsh[l,m,s,phi] <= mcmb[l])
    end
    @constraint(model, [m in 1:2], sum(oaasgsh[l,m,s,phi] for l in 0:NA, s in 0:Smax, phi in 0:1) + azero[m] == 1)
    @constraint(model, [m in 0:2], sum(oabsgsh[l,m,s,phi] for l in 0:NA, s in 0:Smax, phi in 0:1) + bzero[m] == 1)

    if NA > 0
        # Odd
        @variable(model, 0 <= force_odd[1:NA] <= div(maximum_value, 2), Int)
        @constraint(model, [a in 1:NA], ca[a] == 2*force_odd[a]+1)
        @variable(model, Psias[1:NA, Smin:0], Bin)
        @constraint(model, [a in 1:NA, s in Smin:0], ca_no_shift[a] >= 2^(-s)*ca[a] + (Psias[a,s] - 1)*(maximum_value*(2^(-s))))
        @constraint(model, [a in 1:NA, s in Smin:0], ca_no_shift[a] <= 2^(-s)*ca[a] + (1 - Psias[a,s])*(maximum_value*(2^(-s))))
        @constraint(model, [a in 1:NA], sum(Psias[a,s] for s in Smin:0) == 1)
        @constraint(model, [a in 1:NA], phias[a,0] <= sum(Psias[a,s] for s in Smin:-1))
    end

    if !avoid_internal_shifts || NA == 0
        @objective(model, Max, sum(azero[m] for m in 1:2) + sum(bzero[m] for m in 0:2))
    else
        @objective(model, Max, NA*sum(azero[m] for m in 1:2) + NA*sum(bzero[m] for m in 0:2) - sum(Psias[a,s] for a in 1:NA, s in Smin:-1))
    end

    model[:ToLinkZeros] = [[azero[m] for m in 1:2]; [bzero[m] for m in 0:2]]

    return model
end


function mcm_ilp_odd_increment!(model::Model,
                                wordlength::Int,
                                specifications::Vector{Tuple{Float64, Float64, Float64}},
                                bounds_a::Vector{Tuple{Int, Int}},
                                bounds_b::Vector{Tuple{Int, Int}},
                                ;verbose::Bool = false,
                                use_big_m::Bool=true,
                                nb_adders_lb::Int = 0,
                                nb_total_adders_ub::Int = typemax(Int),
                                avoid_internal_shifts::Bool = false,
                                only_positive::Bool=false,
                                no_output_shifts::Bool=false)
    if nb_adders_lb > nb_total_adders_ub
        @warn "Parameters are inconsistent"
    end
    NA = max(0, nb_adders_lb)-1
    scaling = model[:scaling]
    empty!(model)
    solved_once = false
    timelimit = time_limit_sec(model)
    total_solve_time = 0.0
    while (!solved_once) || (termination_status(model) in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED] && timelimit > 0 && nb_total_adders_ub >= NA+1)
        NA += 1
        empty!(model)
        model[:scaling] = scaling
        iir_design!(model, specifications, wordlength, use_big_m=use_big_m, verbose=verbose)
        fix_to_bounds!(model, bounds_a, bounds_b)

        C = model[:ToLink]
        mcm_ilp_odd!(model, C, wordlength, NA, avoid_internal_shifts=avoid_internal_shifts, only_positive=only_positive, no_output_shifts=no_output_shifts)
        if nb_total_adders_ub < sum(length.(C))+-1+NA # Number of structural adders + multiplier adders
            verbose && println("Add constraint on the total number of adders")
            @constraint(model, sum(1-model[:ToLinkZeros][m] for m in 1:length(model[:ToLinkZeros]))-1+NA <= nb_total_adders_ub)
        end
        set_silent(model)
        set_time_limit_sec(model, timelimit)
        optimize!(model)
        solved_once = true

        current_solve_time = solve_time(model)
        total_solve_time += current_solve_time
        timelimit -= current_solve_time
        verbose && println("$(termination_status(model)) for NA=$(NA)")
    end
    model[:NA] = NA
    if !solved_once
        model[:NA] += 1
    else
        verbose && println("\tSolving time: $(total_solve_time)")
    end
    model[:scaling] = scaling

    return model
end
