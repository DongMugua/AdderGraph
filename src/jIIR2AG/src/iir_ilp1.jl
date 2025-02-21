function iir_ilp1!(model::Model,
                   wordlength::Int,
                   specifications_init::Vector{Tuple{Float64, Float64, Float64}},
                   ;presolve_order::Int=1, # 1 => a,b | 2 => b,a
                   presolve_time_sec::Float64=Inf,
                   use_big_m::Bool=true,
                   nb_adders_lb::Int=0,
                   with_symmetry_breaking::Bool=true,
                   avoid_internal_shifts::Bool = false,
                   verbose::Bool = false,
                   debug::Bool = false,
    )
    solution = Iir2AdderGraphs((0,0), (0,0,0), (0,0), AdderGraph(Vector{AdderNode}()), AdderGraph(Vector{AdderNode}()))
    # Call get_gb_secondorder_iir
    verbose && println("Finding bounds for b_k")
    gb = get_gb_secondorder_iir(model, specifications_init, verbose=verbose)
    if !model[:gb_obtained]
        return solution
    end
    infeasibility_proven = true
    scaling_init = (1, gb) # stability -- a_k in [-2,2]
    best_A = typemax(Int)

    scalings = union([(scaling_a, scaling_init[2]) for scaling_a in scaling_init[1]:-1:(scaling_init[1]-wordlength)], [(scaling_init[1], scaling_b) for scaling_b in (scaling_init[2]-1):-1:(scaling_init[2]-wordlength)])
    for scaling in scalings
        empty!(model)
        model[:presolve_done] = false
        model[:success_presolve] = false
        model[:scaling] = scaling
        presolve_done, success_presolve = model[:presolve_done], model[:success_presolve]
        current_AM = nb_adders_lb
        verbose && println("Scaling (a, b): $(scaling)")
        specifications = scale_specifications(specifications_init, scaling)

        while (!presolve_done || success_presolve) && current_AM < best_A
            if !presolve_done
                # Call iir_design!
                verbose && println("Generation of IIR design part")
                iir_design!(model, specifications, wordlength, use_big_m=use_big_m, with_symmetry_breaking=with_symmetry_breaking, verbose=verbose)

                presolve_done = true
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
                    if !model[:infeasibility_proven]
                        infeasibility_proven = false
                    end
                    continue
                else
                    infeasibility_proven = false
                end
                second_presolve!(model, presolve_time_sec_for_each_coef=presolve_time_sec_for_each_coef, verbose=verbose)
                if !model[:success_presolve]
                    if !model[:infeasibility_proven]
                        infeasibility_proven = false
                    end
                    continue
                else
                    infeasibility_proven = false
                end
                @assert model[:success_presolve]
                success_presolve = model[:success_presolve]
                verbose && println(model[:bounds_a])
                verbose && println(model[:bounds_b])
            end
            bounds_a, bounds_b = model[:bounds_a], model[:bounds_b]

            # Add MCM ILP1 to model
            verbose && println("Adding MCM to the model and solving...")
            mcm_ilp_odd_increment!(model, wordlength, specifications,
                                   model[:bounds_a], model[:bounds_b], nb_adders_lb=current_AM,
                                   use_big_m=use_big_m, avoid_internal_shifts=avoid_internal_shifts,
                                   verbose=verbose, only_positive=true, no_output_shifts=false,
                                   nb_total_adders_ub=best_A-1)

            current_AM = model[:NA]+1
            if termination_status(model) == MOI.OPTIMAL
                NA = model[:NA]
                best_A = NA+sum(1-round(Int, value(model[:azero][m])) for m in 1:2)+sum(1-round(Int, value(model[:bzero][m])) for m in 0:2)-1
                verbose && println("Current solution -- A: $(best_A)")
                solution.coefficients_a = (round(Int, value(model[:a][1])), round(Int, value(model[:a][2])))
                solution.coefficients_b = (round(Int, value(model[:b][0])), round(Int, value(model[:b][1])), round(Int, value(model[:b][2])))
                verbose && println("\tcoefficients b: ", solution.coefficients_b)
                verbose && println("\tcoefficients a: ", solution.coefficients_a)
                solution.shifts = (wordlength-1-model[:scaling][1]+round(Int, log2(value(model[:a0]))), wordlength-1-model[:scaling][2]+round(Int, log2(value(model[:a0]))))
                verbose && println("\tshifts (a,b): ", solution.shifts)
                solution.addergraph_a = AdderGraph(Vector{AdderNode}(), collect(solution.coefficients_a))
                for i in 1:NA
                    if round(Int, value(model[:mcma][i])) == 1
                        node_shift = 0
                        for s in -wordlength:0
                            if round(Int, value(model[:Psias][i,s])) == 1
                                node_shift = s
                                break
                            end
                        end
                        input_shift = 0
                        for s in 0:wordlength
                            if round(Int, value(model[:phias][i,s])) == 1
                                input_shift = s
                                break
                            end
                        end
                        subtraction = [value(model[:cai_right_shsg][i]) < 0, value(model[:cai_left_sg][i]) < 0]
                        push_node!(solution.addergraph_a,
                            AdderNode(round(Int, value(model[:ca][i])),
                                      [InputEdge(get_addernodes_by_value(solution.addergraph_a, round(Int, value(model[:cai][i,1])))[end], input_shift+node_shift, subtraction[1]),
                                      InputEdge(get_addernodes_by_value(solution.addergraph_a, round(Int, value(model[:cai][i,2])))[end], node_shift, subtraction[2])]
                            ))
                    end
                end
                solution.addergraph_b = AdderGraph(Vector{AdderNode}(), collect(solution.coefficients_b))
                for i in 1:NA
                    if round(Int, value(model[:mcmb][i])) == 1
                        node_shift = 0
                        for s in -wordlength:0
                            if round(Int, value(model[:Psias][i,s])) == 1
                                node_shift = s
                                break
                            end
                        end
                        input_shift = 0
                        for s in 0:wordlength
                            if round(Int, value(model[:phias][i,s])) == 1
                                input_shift = s
                                break
                            end
                        end
                        subtraction = [value(model[:cai_right_shsg][i]) < 0, value(model[:cai_left_sg][i]) < 0]
                        push_node!(solution.addergraph_b,
                            AdderNode(round(Int, value(model[:ca][i])),
                                      [InputEdge(get_addernodes_by_value(solution.addergraph_b, round(Int, value(model[:cai][i,1])))[end], input_shift+node_shift, subtraction[1]),
                                      InputEdge(get_addernodes_by_value(solution.addergraph_b, round(Int, value(model[:cai][i,2])))[end], node_shift, subtraction[2])]
                            ))
                    end
                end
            end
            empty!(model)
            model[:scaling] = scaling
            model[:bounds_a], model[:bounds_b] = bounds_a, bounds_b
            model[:presolve_done], model[:success_presolve] = presolve_done, success_presolve
        end
    end

    if infeasibility_proven
        @warn "No second-order IIR filter for these specifications"
    end

    return solution
end
