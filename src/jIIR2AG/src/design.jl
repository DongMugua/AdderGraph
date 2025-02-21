function all_coefficients(
        model::Model,
        specifications_init::Vector{Tuple{Float64, Float64, Float64}},
        wordlength::Int;
        with_symmetry_breaking::Bool=true,
        presolve::Int=1, # default: 1 => presolve a then b -- 2 => b then a, 0 => no presolve (incompatible with ilp1)
        presolve_time_sec::Union{Float64, Int}=Inf,
        use_big_m::Bool=true,
        nb_max_solutions::Union{Float64, Int}=Inf,
        verbose::Bool=false,
    )
    presolve_order = presolve
    presolve_time_sec = Float64(presolve_time_sec)
    solutions = Vector{Iir2AdderGraphs}()

    verbose && println("Finding bounds for b_k")
    gb = get_gb_secondorder_iir(model, specifications_init, verbose=verbose)
    if !model[:gb_obtained]
        return solutions
    end
    scaling_init = (1, gb) # stability -- a_k in [-2,2]

    scalings = union([(scaling_a, scaling_init[2]) for scaling_a in scaling_init[1]:-1:(scaling_init[1]-wordlength)], [(scaling_init[1], scaling_b) for scaling_b in (scaling_init[2]-1):-1:(scaling_init[2]-wordlength)])
    for scaling in scalings
        empty!(model)
        verbose && println("Scaling (a, b): $(scaling)")
        specifications = scale_specifications(specifications_init, scaling)

        # Call iir_design!
        verbose && println("Generation of IIR design part")
        model[:scaling] = scaling
        iir_design!(model, specifications, wordlength, with_symmetry_breaking=with_symmetry_breaking, use_big_m=use_big_m, verbose=verbose)

        # Presolve? a_k? b_k?
        verbose && println("Presolving")
        first_presolve! = presolve_a!
        second_presolve! = presolve_b!
        if presolve_order == 2
            first_presolve! = presolve_b!
            second_presolve! = presolve_a!
        end
        presolve_time_sec_for_each_coef = presolve_time_sec/5

        timelimit = time_limit_sec(model)
        set_time_limit_sec(model, timelimit)

        first_presolve!(model, presolve_time_sec_for_each_coef=presolve_time_sec_for_each_coef, verbose=verbose)
        if !model[:success_presolve]
            continue
        end
        second_presolve!(model, presolve_time_sec_for_each_coef=presolve_time_sec_for_each_coef, verbose=verbose)
        @assert model[:success_presolve]
        bounds_a, bounds_b = model[:bounds_a], model[:bounds_b]
        verbose && println(model[:bounds_a])
        verbose && println(model[:bounds_b])

        set_silent(model)
        nb_solution = 0
        verbosemultiplier = 1

        a1_lower_bound = round(Int, lower_bound(model[:a][1]))
        a1_upper_bound = round(Int, upper_bound(model[:a][1]))
        a2_lower_bound = round(Int, lower_bound(model[:a][2]))
        a2_upper_bound = round(Int, upper_bound(model[:a][2]))
        b0_lower_bound = round(Int, lower_bound(model[:b][0]))
        b0_upper_bound = round(Int, upper_bound(model[:b][0]))
        b1_lower_bound = round(Int, lower_bound(model[:b][1]))
        b1_upper_bound = round(Int, upper_bound(model[:b][1]))
        b2_lower_bound = round(Int, lower_bound(model[:b][2]))
        b2_upper_bound = round(Int, upper_bound(model[:b][2]))
        current_a1_lower_bound = round(Int, lower_bound(model[:a][1]))
        current_a1_upper_bound = round(Int, upper_bound(model[:a][1]))
        for a1 in current_a1_lower_bound:current_a1_upper_bound
            fix(model[:a][1], a1, force=true)
            @objective(model, Min, model[:a][2])
            optimize!(model)
            if termination_status(model) != MOI.OPTIMAL
                continue
            end
            current_a2_lower_bound = round(Int, value(model[:a][2]))
            @objective(model, Max, model[:a][2])
            optimize!(model)
            current_a2_upper_bound = round(Int, value(model[:a][2]))
            for a2 in current_a2_lower_bound:current_a2_upper_bound
                fix(model[:a][2], a2, force=true)
                @objective(model, Min, model[:b][0])
                optimize!(model)
                if termination_status(model) != MOI.OPTIMAL
                    continue
                end
                current_b0_lower_bound = round(Int, value(model[:b][0]))
                @objective(model, Max, model[:b][0])
                optimize!(model)
                current_b0_upper_bound = round(Int, value(model[:b][0]))
                for b0 in current_b0_lower_bound:current_b0_upper_bound
                    fix(model[:b][0], b0, force = true)
                    @objective(model, Min, model[:b][1])
                    optimize!(model)
                    if termination_status(model) != MOI.OPTIMAL
                        continue
                    end
                    current_b1_lower_bound = round(Int, value(model[:b][1]))
                    @objective(model, Max, model[:b][1])
                    optimize!(model)
                    current_b1_upper_bound = round(Int, value(model[:b][1]))
                    for b1 in current_b1_lower_bound:current_b1_upper_bound
                        fix(model[:b][1], b1, force = true)
                        @objective(model, Min, model[:b][2])
                        optimize!(model)
                        if termination_status(model) != MOI.OPTIMAL
                            continue
                        end
                        current_b2_lower_bound = round(Int, value(model[:b][2]))
                        @objective(model, Max, model[:b][2])
                        optimize!(model)
                        current_b2_upper_bound = round(Int, value(model[:b][2]))
                        set_lower_bound(model[:b][2], b2_lower_bound)
                        set_upper_bound(model[:b][2], b2_upper_bound)
                        for b2 in current_b2_lower_bound:current_b2_upper_bound
                            for a0 in 0:wordlength
                                if is_stable(a1/(1<<(wordlength-1-scaling[1]+a0)), a2/(1<<(wordlength-1-scaling[1]+a0))) &&
                                    coefs_fit_specs(a1/(1<<(wordlength-1-scaling[1]+a0)), a2/(1<<(wordlength-1-scaling[1]+a0)),
                                        b0/(1<<(wordlength-1-scaling[2]+a0)), b1/(1<<(wordlength-1-scaling[2]+a0)), b2/(1<<(wordlength-1-scaling[2]+a0)),
                                        specifications_init)
                                    push!(solutions, Iir2AdderGraphs((a1, a2), (b0, b1, b2), (wordlength-1-scaling[1]+a0, wordlength-1-scaling[2]+a0)))
                                    nb_solution += 1
                                    if nb_solution < 10*verbosemultiplier && nb_solution % (1*verbosemultiplier) == 0
                                        verbose && println("nb_solution: $nb_solution")
                                    elseif nb_solution == 10*verbosemultiplier
                                        verbose && println("nb_solution: $nb_solution")
                                        verbosemultiplier *= 10
                                    end
                                    if nb_solution >= nb_max_solutions
                                        return solutions
                                    end
                                end
                            end
                        end
                    end
                    unfix(model[:b][1])
                    set_lower_bound(model[:b][1], b1_lower_bound)
                    set_upper_bound(model[:b][1], b1_upper_bound)
                end
                unfix(model[:b][0])
                set_lower_bound(model[:b][0], b0_lower_bound)
                set_upper_bound(model[:b][0], b0_upper_bound)
            end
            unfix(model[:a][2])
            set_lower_bound(model[:a][2], a2_lower_bound)
            set_upper_bound(model[:a][2], a2_upper_bound)
        end
        unfix(model[:a][1])
        set_lower_bound(model[:a][1], a1_lower_bound)
        set_upper_bound(model[:a][1], a1_upper_bound)
    end

    return solutions
end

all_coefficients(optimizer::DataType, args...; kwargs...) = all_coefficients(Model(optimizer), args...; kwargs...)

all_coefficients(model::Model,
                 fbands::Vector{Tuple{Float64, Float64}},
                 dbands::Union{Vector{Float64},Vector{Int}},
                 delta::Float64,
                 wordlength::Int;
                 size_of_grid::Int=480,
                 kwargs...
    ) = all_coefficients(model, get_specifications(fbands, dbands, delta, size_of_grid), wordlength; kwargs...)

all_coefficients(model::Model,
                 fbands::Vector{Tuple{Float64, Float64}},
                 dbands::Union{Vector{Float64},Vector{Int}},
                 deltas::Vector{Float64},
                 wordlength::Int;
                 size_of_grid::Int=480,
                 kwargs...
    ) = all_coefficients(model, get_specifications(fbands, dbands, deltas, size_of_grid), wordlength; kwargs...)

all_coefficients(model::Model,
                 fbands::Vector{Tuple{Float64, Float64}},
                 abands::Vector{Tuple{Float64, Float64}},
                 wordlength::Int;
                 size_of_grid::Int=480,
                 kwargs...
    ) = all_coefficients(model, get_specifications(fbands, abands, size_of_grid), wordlength; kwargs...)

all_coefficients(model::Model,
                 transfer_function::Function,
                 error_margin_percent::Union{Float64, Int},
                 error_minimum::Float64,
                 wordlength::Int;
                 size_of_grid::Int=480,
                 kwargs...
    ) = all_coefficients(model, get_specifications(transfer_function, error_margin_percent, error_minimum, size_of_grid), wordlength; kwargs...)



function design_second_order_iir(
        model::Model,
        specifications::Vector{Tuple{Float64, Float64, Float64}},
        wordlength::Int;
        ilp::Int=1, # default: 1 => ilp1 -- 2 => ilp2
        use_big_m::Bool=true,
        avoid_internal_shifts::Bool=false, # Only implemented for ILP1 yet
        with_symmetry_breaking::Bool=true, # Only implemented for ILP1 yet
        max_depth::Int=0, # default: 0 => no bound (incompatible with ilp2)
        presolve::Int=1, # default: 1 => presolve a then b -- 2 => b then a, 0 => no presolve (incompatible with ilp1)
        presolve_time_sec::Union{Float64, Int}=Inf,
        nb_adders_lb::Int = 0, # 0: automatically compute
        verbose::Bool=false, debug::Bool=true,
    )
    @assert ilp in [1,2]
    @assert max_depth >= 0
    @assert presolve in [0,1,2]
    @assert ilp != 1 || presolve != 0
    @assert ilp != 2 || max_depth != 0
    presolve_time_sec = Float64(presolve_time_sec)

    tmp_time = @timed begin
    if ilp == 1
        if use_big_m && wordlength >= 10
            @warn "Default integrality tolerance might lead to numerical instability"
        end
        solution = iir_ilp1!(model, wordlength, specifications, nb_adders_lb=nb_adders_lb, presolve_order=presolve, use_big_m=use_big_m, avoid_internal_shifts=avoid_internal_shifts, presolve_time_sec=presolve_time_sec, with_symmetry_breaking=with_symmetry_breaking, verbose=verbose, debug=debug)
    else #ilp = 2
        solution = iir_ilp2!(model, wordlength, specifications, presolve_order=presolve, presolve_time_sec=presolve_time_sec, max_depth=max_depth)
        if termination_status(model) == MOI.OPTIMAL

        end
    end
    end #@timed
    verbose && println("Solving time: $(tmp_time[2])")

    return solution
end

design_second_order_iir(optimizer::DataType, args...; kwargs...) = design_second_order_iir(Model(optimizer), args...; kwargs...)

design_second_order_iir(model::Model,
                        fbands::Vector{Tuple{Float64, Float64}},
                        abands::Vector{Tuple{Float64, Float64}},
                        wordlength::Int;
                        size_of_grid::Int=480,
                        kwargs...
    ) = design_second_order_iir(model, get_specifications(fbands, abands, size_of_grid), wordlength; kwargs...)

design_second_order_iir(model::Model,
                        fbands::Vector{Tuple{Float64, Float64}},
                        dbands::Union{Vector{Float64},Vector{Int}},
                        deltas::Vector{Float64},
                        wordlength::Int;
                        size_of_grid::Int=480,
                        kwargs...
    ) = design_second_order_iir(model, get_specifications(fbands, dbands, deltas, size_of_grid), wordlength; kwargs...)

design_second_order_iir(model::Model,
                        fbands::Vector{Tuple{Float64, Float64}},
                        dbands::Union{Vector{Float64},Vector{Int}},
                        delta::Float64,
                        wordlength::Int;
                        size_of_grid::Int=480,
                        kwargs...
    ) = design_second_order_iir(model, get_specifications(fbands, dbands, delta, size_of_grid), wordlength; kwargs...)

design_second_order_iir(model::Model,
                        transfer_function::Function,
                        error_margin_percent::Union{Float64, Int},
                        error_minimum::Float64,
                        wordlength::Int;
                        size_of_grid::Int=480,
                        kwargs...
    ) = design_second_order_iir(model, get_specifications(transfer_function, error_margin_percent, error_minimum, size_of_grid), wordlength; kwargs...)
