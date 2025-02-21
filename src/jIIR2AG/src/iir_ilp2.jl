function iir_ilp2!(model::Model,
                   wordlength::Int,
                   specifications::Vector{Tuple{Float64, Float64, Float64}},
                   ;presolve_order::Int=1, # 1 => a,b | 2 => b,a | 0 no presolve
                   max_depth::Int=3,
                   verbose::Bool = false)
    # Call get_gb_secondorder_iir
    verbose && println("Finding bounds for b_k")
    gb = get_gb_secondorder_iir(model, specifications)
    if !model[:gb_obtained]
        return solution
    end
    empty!(model)
    scaling = (1, gb) # stability -- a_k in [-2,2]
    specifications = scale_specifications(specifications, scaling)

    # Call iir_design!
    verbose && println("Generation of IIR design part")
    model[:scaling] = scaling
    iir_design!(model, specifications, wordlength)

    # Presolve? a_k? b_k?
    if presolve_order != 0
        verbose && println("Presolving")
        first_presolve! = presolve_a!
        second_presolve! = presolve_b!
        if presolve_order == 2
            first_presolve! = presolve_b!
            second_presolve! = presolve_a!
        end
        if first_presolve!(model) == -1
            return model
        end
        second_presolve!(model)
    end

    # Add MCM ILP2 to model
    verbose && println("Adding MCM to the model and solving...")
    mcm_ilp_2!(model, model[:ToLink], wordlength, max_depth)

    return model
end
