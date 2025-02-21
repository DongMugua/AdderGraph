function presolve_a!(model::Model; presolve_time_sec_for_each_coef::Float64=Inf, verbose::Bool=false)
    model[:infeasibility_proven] = false
    set_silent(model)
    tmp_time_limit = time_limit_sec(model)
    if isnothing(tmp_time_limit)
        tmp_time_limit = Inf
    end
    if presolve_time_sec_for_each_coef != Inf
        set_time_limit_sec(model, presolve_time_sec_for_each_coef/2)
    end
    time_for_presolve_a = 0.0
    model[:bounds_a] = Vector{Tuple{Int, Int}}()
    for m in 1:2
        @objective(model, Min, model[:a][m])
        optimize!(model)
        time_for_presolve_a += solve_time(model)
        if !has_values(model)
            if termination_status(model) == MOI.INFEASIBLE
                model[:infeasibility_proven] = true
            end
            model[:success_presolve] = false
            unset_silent(model)
            tmp_time_limit -= time_for_presolve_a
            set_time_limit_sec(model, tmp_time_limit)
            return model
        end
        a_lower_bound = 0
        if termination_status(model) == MOI.TIME_LIMIT
            a_lower_bound = round(Int, objective_bound(model))
            set_lower_bound(model[:a][m], a_lower_bound)
            #a_lower_bound = round(Int, lower_bound(model[:a][m]))
        else
            a_lower_bound = round(Int, value(model[:a][m]))
            set_lower_bound(model[:a][m], a_lower_bound)
        end
        @objective(model, Max, model[:a][m])
        optimize!(model)
        time_for_presolve_a += solve_time(model)
        if !has_values(model)
            if termination_status(model) == MOI.INFEASIBLE
                model[:infeasibility_proven] = true
            end
            model[:success_presolve] = false
            unset_silent(model)
            tmp_time_limit -= time_for_presolve_a
            set_time_limit_sec(model, tmp_time_limit)
            return model
        end
        a_upper_bound = 0
        if termination_status(model) == MOI.TIME_LIMIT
            a_upper_bound = round(Int, objective_bound(model))
            set_upper_bound(model[:a][m], a_upper_bound)
            #a_upper_bound = round(Int, upper_bound(model[:a][m]))
        else
            a_upper_bound = round(Int, value(model[:a][m]))
            set_upper_bound(model[:a][m], a_upper_bound)
        end
        push!(model[:bounds_a], (a_lower_bound, a_upper_bound))
        if 0 < a_lower_bound
            set_lower_bound(model[:asquare][m], a_lower_bound^2)
        elseif 0 > a_upper_bound
            set_lower_bound(model[:asquare][m], a_upper_bound^2)
        end
        set_upper_bound(model[:asquare][m], max(a_lower_bound^2, a_upper_bound^2))
        for n in 1:(m-1)
            if 0 < a_lower_bound && 0 < lower_bound(model[:a][n])
                set_lower_bound(model[:aquad][n,m], a_lower_bound*lower_bound(model[:a][n]))
                set_upper_bound(model[:aquad][n,m], a_upper_bound*upper_bound(model[:a][n]))
            elseif 0 > a_upper_bound && 0 > upper_bound(model[:a][n])
                set_lower_bound(model[:aquad][n,m], a_upper_bound*upper_bound(model[:a][n]))
                set_upper_bound(model[:aquad][n,m], a_lower_bound*lower_bound(model[:a][n]))
            elseif 0 < a_lower_bound && 0 > upper_bound(model[:a][n])
                set_lower_bound(model[:aquad][n,m], a_upper_bound*lower_bound(model[:a][n]))
                set_upper_bound(model[:aquad][n,m], a_lower_bound*upper_bound(model[:a][n]))
            elseif 0 > a_upper_bound && 0 < lower_bound(model[:a][n])
                set_lower_bound(model[:aquad][n,m], a_lower_bound*upper_bound(model[:a][n]))
                set_upper_bound(model[:aquad][n,m], a_upper_bound*lower_bound(model[:a][n]))
            end
        end
        verbose && println("\tmodel[:a][$m] in ", (a_lower_bound, a_upper_bound))
    end
    unset_silent(model)
    tmp_time_limit -= time_for_presolve_a
    if !isinf(tmp_time_limit)
        set_time_limit_sec(model, tmp_time_limit)
    end

    model[:success_presolve] = true
    return model
end

function presolve_b!(model::Model; presolve_time_sec_for_each_coef::Float64=Inf, verbose::Bool=false)
    model[:infeasibility_proven] = false
    set_silent(model)
    tmp_time_limit = time_limit_sec(model)
    if presolve_time_sec_for_each_coef != Inf
        set_time_limit_sec(model, presolve_time_sec_for_each_coef/2)
    end
    time_for_presolve_b = 0.0
    model[:bounds_b] = Vector{Tuple{Int, Int}}()
    for m in 0:2
        @objective(model, Min, model[:b][m])
        optimize!(model)
        time_for_presolve_b += solve_time(model)
        if !has_values(model)
            if termination_status(model) == MOI.INFEASIBLE
                model[:infeasibility_proven] = true
            end
            model[:success_presolve] = false
            unset_silent(model)
            tmp_time_limit -= time_for_presolve_b
            set_time_limit_sec(model, tmp_time_limit)
            return model
        end
        b_lower_bound = 0
        if termination_status(model) == MOI.TIME_LIMIT
            b_lower_bound = round(Int, objective_bound(model))
            set_lower_bound(model[:b][m], b_lower_bound)
            #b_lower_bound = round(Int, lower_bound(model[:b][m]))
        else
            b_lower_bound = round(Int, value(model[:b][m]))
            set_lower_bound(model[:b][m], b_lower_bound)
        end
        @objective(model, Max, model[:b][m])
        optimize!(model)
        time_for_presolve_b += solve_time(model)
        if !has_values(model)
            if termination_status(model) == MOI.INFEASIBLE
                model[:infeasibility_proven] = true
            end
            model[:success_presolve] = false
            unset_silent(model)
            tmp_time_limit -= time_for_presolve_a
            set_time_limit_sec(model, tmp_time_limit)
            return model
        end
        b_upper_bound = 0
        if termination_status(model) == MOI.TIME_LIMIT
            b_upper_bound = round(Int, objective_bound(model))
            set_upper_bound(model[:b][m], b_upper_bound)
            #b_upper_bound = round(Int, upper_bound(model[:b][m]))
        else
            b_upper_bound = round(Int, value(model[:b][m]))
            set_upper_bound(model[:b][m], b_upper_bound)
        end
        push!(model[:bounds_b], (b_lower_bound, b_upper_bound))
        if 0 < b_lower_bound
            set_lower_bound(model[:bsquare][m], b_lower_bound^2)
        elseif 0 > b_upper_bound
            set_lower_bound(model[:bsquare][m], b_upper_bound^2)
        end
        set_upper_bound(model[:bsquare][m], max(b_lower_bound^2, b_upper_bound^2))
        for n in 0:(m-1)
            if 0 < b_lower_bound && 0 < lower_bound(model[:b][n])
                set_lower_bound(model[:bquad][n,m], b_lower_bound*lower_bound(model[:b][n]))
                set_upper_bound(model[:bquad][n,m], b_upper_bound*upper_bound(model[:b][n]))
            elseif 0 > b_upper_bound && 0 > upper_bound(model[:b][n])
                set_lower_bound(model[:bquad][n,m], b_upper_bound*upper_bound(model[:b][n]))
                set_upper_bound(model[:bquad][n,m], b_lower_bound*lower_bound(model[:b][n]))
            elseif 0 < b_lower_bound && 0 > upper_bound(model[:b][n])
                set_lower_bound(model[:bquad][n,m], b_upper_bound*lower_bound(model[:b][n]))
                set_upper_bound(model[:bquad][n,m], b_lower_bound*upper_bound(model[:b][n]))
            elseif 0 > b_upper_bound && 0 < lower_bound(model[:b][n])
                set_lower_bound(model[:bquad][n,m], b_lower_bound*upper_bound(model[:b][n]))
                set_upper_bound(model[:bquad][n,m], b_upper_bound*lower_bound(model[:b][n]))
            end
        end
        verbose && println("\tmodel[:b][$m] in ", (b_lower_bound, b_upper_bound))
    end
    unset_silent(model)
    tmp_time_limit -= time_for_presolve_b
    set_time_limit_sec(model, tmp_time_limit)

    model[:success_presolve] = true
    return model
end
