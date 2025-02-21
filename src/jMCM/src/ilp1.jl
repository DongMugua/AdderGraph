#using JuMP

"""
    model_mcm_formulation_1_odd_big_m!(model::Model, C::Vector{Int},
                                       wordlength::Int, S::Tuple{Int, Int},
                                       NA::Int, _adder_depth_max::Int)::Model

Formulation 1 from "M. Kumm -- Optimal Constant Multiplication Using Integer
Linear Programming" with small modifications: no left shifts and only positive
right shifts
"""
function model_mcm_formulation_1_odd_big_m!(model::Model, C::Vector{Int},
                                            wordlength::Int, S::Tuple{Int, Int},
                                            NA::Int, adder_depth_max::Int,
                                            avoid_internal_shifts::Bool)::Model
    Smin, Smax = S
    NO = length(C)
    maximum_target = maximum(C)
    maximum_value = 2^wordlength

    @variable(model, 1 <= ca[0:NA] <= maximum_value-1, Int)
    @variable(model, 1 <= ca_no_shift[1:NA] <= maximum_value-1, Int)
    @variable(model, 1 <= cai[1:NA, 1:2] <= maximum_value-1, Int)
    @variable(model, 1 <= ca_right_sh[1:NA] <= maximum_value-1, Int)
    @variable(model, -maximum_value <= cai_right_shsg[1:NA] <= maximum_value-1, Int)
    @variable(model, -maximum_value <= cai_left_sg[1:NA] <= maximum_value-1, Int)

    @variable(model, Phiai[1:NA, 1:2], Bin)
    @variable(model, caik[a in 1:NA, 1:2, 0:(a-1)], Bin)
    @variable(model, phias[1:NA, 0:Smax], Bin)
    @variable(model, oaj[0:NA, 1:NO], Bin)

    if adder_depth_max != 0
        # Add adder depth
        @variable(model, 0 <= ada[0:NA] <= NA)
        fix(ada[0], 0, force=true)
        @variable(model, 1 <= max_ad <= NA, Int)
        @constraint(model, [a in 1:NA, i in 1:2, k in 0:(a-1)], ada[a] >= ada[k]+1 - (1-caik[a,i,k])*NA)
        @constraint(model, [a in 1:NA], max_ad >= ada[a])

        @constraint(model, max_ad <= adder_depth_max)
    end

    # C1
    fix(ca[0], 1, force=true)
    # C2 - Modified
    @constraint(model, [a in 1:NA], ca_no_shift[a] == cai_left_sg[a] + cai_right_shsg[a])
    # C3a - C3b
    @constraint(model, [a in 1:NA, i in 1:2, k in 0:(a-1)], cai[a,i] <= ca[k] + (1-caik[a,i,k])*maximum_value)
    @constraint(model, [a in 1:NA, i in 1:2, k in 0:(a-1)], cai[a,i] >= ca[k] - (1-caik[a,i,k])*maximum_value)
    @constraint(model, [a in 1:NA, i in 1:2], sum(caik[a,i,k] for k in 0:(a-1)) == 1)
    # C4a - C4b - Modified
    @constraint(model, [a in 1:NA, s in 0:Smax], ca_right_sh[a] <= 2^s*cai[a,1] + (1-phias[a,s])*maximum_value)
    @constraint(model, [a in 1:NA, s in 0:Smax], ca_right_sh[a] >= 2^s*cai[a,1] - (1-phias[a,s])*(maximum_value*(2^s)))
    @constraint(model, [a in 1:NA], sum(phias[a,s] for s in 0:Smax) == 1)
    # C5a - C5b - C5c - Modified
    @constraint(model, [a in 1:NA], cai_right_shsg[a] <= ca_right_sh[a] + Phiai[a,1]*maximum_value)
    @constraint(model, [a in 1:NA], cai_right_shsg[a] >= ca_right_sh[a] - Phiai[a,1]*(2*maximum_value))
    @constraint(model, [a in 1:NA], cai_right_shsg[a] <= -ca_right_sh[a] + (1-Phiai[a,1])*(2*maximum_value))
    @constraint(model, [a in 1:NA], cai_right_shsg[a] >= -ca_right_sh[a] - (1-Phiai[a,1])*maximum_value)
    @constraint(model, [a in 1:NA], cai_left_sg[a] <= cai[a,2] + Phiai[a,2]*maximum_value)
    @constraint(model, [a in 1:NA], cai_left_sg[a] >= cai[a,2] - Phiai[a,2]*(2*maximum_value))
    @constraint(model, [a in 1:NA], cai_left_sg[a] <= -cai[a,2] + (1-Phiai[a,2])*(2*maximum_value))
    @constraint(model, [a in 1:NA], cai_left_sg[a] >= -cai[a,2] - (1-Phiai[a,2])*maximum_value)
    @constraint(model, [a in 1:NA], Phiai[a,1] + Phiai[a,2] <= 1)
    # C6a - C6b
    @constraint(model, [a in 0:NA, j in 1:NO], ca[a] <= C[j] + (1-oaj[a,j])*maximum_value)
    @constraint(model, [a in 0:NA, j in 1:NO], ca[a] >= C[j] - (1-oaj[a,j])*maximum_target)
    @constraint(model, [j in 1:NO], sum(oaj[a,j] for a in 0:NA) == 1)

    # Odd
    @variable(model, 0 <= force_odd[1:NA] <= div(maximum_value, 2), Int)
    @constraint(model, [a in 1:NA], ca[a] == 2*force_odd[a]+1)
    @variable(model, Psias[1:NA, Smin:0], Bin)
    @constraint(model, [a in 1:NA, s in Smin:0], ca_no_shift[a] >= 2^(-s)*ca[a] + (Psias[a,s] - 1)*(maximum_value*(2^(-s))))
    @constraint(model, [a in 1:NA, s in Smin:0], ca_no_shift[a] <= 2^(-s)*ca[a] + (1 - Psias[a,s])*(maximum_value*(2^(-s))))
    @constraint(model, [a in 1:NA], sum(Psias[a,s] for s in Smin:0) == 1)
    @constraint(model, [a in 1:NA], phias[a,0] == sum(Psias[a,s] for s in Smin:-1))

    if !avoid_internal_shifts
        # Mock objective
        @objective(model, Min, NA)
    else
        # Objective that permit to minimize the number of inside shifts
        @objective(model, Min, sum(Psias[a,s] for a in 1:NA, s in Smin:-1))
    end

    return model
end




"""
    optimize_increment!(model::Model, model_mcm_forumlation!::Function,
                        C::Vector{Int}, wordlength::Int, S::Tuple{Int, Int},
                        verbose::Bool)::Model

Increment NA until a solution is found for the coefficients in `C`.
"""
function optimize_increment!(model::Model, model_mcm_forumlation!::Function,
                             C::Vector{Int}, wordlength::Int, S::Tuple{Int, Int},
                             ;verbose::Bool = false, adder_depth_max::Int = 0,
                             nb_adders_lb::Int = 0, avoid_internal_shifts::Bool = false)::Model
    NA = max(length(C), nb_adders_lb)
    model_mcm_forumlation!(model, C, wordlength, S, NA, adder_depth_max, avoid_internal_shifts)
    timelimit = time_limit_sec(model)
    total_solve_time = 0.0
    optimize!(model)
    current_solve_time = solve_time(model)
    total_solve_time += current_solve_time
    while termination_status(model) in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        timelimit -= current_solve_time
        verbose && println("$(termination_status(model)) for NA = $NA in $current_solve_time seconds")
        if timelimit <= 0.0
            break
        end
        NA += 1
        empty!(model)
        model_mcm_forumlation!(model, C, wordlength, S, NA, adder_depth_max, avoid_internal_shifts)
        set_time_limit_sec(model, timelimit)
        optimize!(model)
        current_solve_time = solve_time(model)
        total_solve_time += current_solve_time
    end
    model[:NA] = NA
    verbose && println("\n$(termination_status(model)) for NA = $NA in $current_solve_time seconds\nTotal time: $total_solve_time seconds\n\n\n")

    return model
end
