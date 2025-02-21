#using GLPK

#include("../utils.jl")
#include("mcm/ilp1.jl")
#include("mcm/ilp2.jl")
#include("mcm/nlp.jl")

"""
    mcm(model::Model,
        C::Vector{Int},
        ;wordlength::Int = 0,
        use_nlp::Bool = false,
        ilp::Int = 1, # 1 : ILP1 -- 2 : ILP2
        use_big_m::Bool = true, # Only for ILP1
        avoid_internal_shifts::Bool = false, # Only implemented for ILP1 yet
        adder_depth_max::Int = 0, # 0 to deactivate this bound
        nb_adders_lb::Int = 0,
        verbose::Bool = false,
    )

Add to `model`, an empty model with a solver attached, an mcm modelization.
The model is then optimized and the solution is transposed into an AdderGraph.
"""
function mcm(model::Model,
             C::Vector{Int},
             ;wordlength::Int = 0,
             use_nlp::Bool = false,
             ilp::Int = 1, # 1 : ILP1 -- 2 : ILP2
             use_big_m::Bool = true, # Only for ILP1
             avoid_internal_shifts::Bool = false, # Only implemented for ILP1 yet
             adder_depth_max::Int = 0, # 0 to deactivate this bound
             nb_adders_lb::Int = 0, # 0: automatically compute
             verbose::Bool = false,
    )
    if wordlength == 0
        return mcm(model, C, wordlength = maximum(get_min_wordlength.(C)), nb_adders_lb=nb_adders_lb, use_nlp=use_nlp, ilp=ilp, use_big_m=use_big_m, adder_depth_max=adder_depth_max, verbose=verbose)
    end
    if isempty(C)
        return AdderGraph()
    end
    oddabsC = filter!(x -> x > 1, unique!(odd.(abs.(C))))
    if isempty(oddabsC)
        return AdderGraph(C)
    end
    if (1 << wordlength) - 1 < maximum(oddabsC)
        return AdderGraph()
    end
    !verbose && set_silent(model)
    addergraph = AdderGraph(C)
    model_mcm_forumlation! = model_mcm_formulation_1_odd_big_m!

    if use_nlp
        error("Not implemented yet.")
        model_mcm_forumlation! = model_mcm_nlp!
    else
        if ilp == 1
            if use_big_m
                model_mcm_forumlation! = model_mcm_formulation_1_odd_big_m!
            else
                error("Not implemented yet.")
                if adder_depth_max == 0
                    model_mcm_forumlation! = model_mcm_formulation_1_odd_ind!
                end
            end
            optimize_increment!(model, model_mcm_forumlation!, oddabsC, wordlength, (-wordlength,wordlength), adder_depth_max=adder_depth_max, avoid_internal_shifts=avoid_internal_shifts, nb_adders_lb=nb_adders_lb, verbose=verbose)

            if termination_status(model) == MOI.OPTIMAL
                for i in 1:model[:NA]
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
                    push_node!(addergraph,
                        AdderNode(round(Int, value(model[:ca][i])),
                                  [InputEdge(get_addernodes_by_value(addergraph, round(Int, value(model[:cai][i,1])))[end], input_shift+node_shift, subtraction[1]),
                                  InputEdge(get_addernodes_by_value(addergraph, round(Int, value(model[:cai][i,2])))[end], node_shift, subtraction[2])]
                        )
                    )
                end
            end

        elseif ilp == 2
            if adder_depth_max == 0
                # Default for ILP1 should be modified for ILP2
                adder_depth_max = 3
            end
            model_mcm_formulation_2!(model, C, construct_hypertree(wordlength, adder_depth_max), adder_depth_max)
            optimize!(model)
            nodes_by_depth = Vector{Vector{HypernodeID}}()
            for depth in 1:adder_depth_max
                push!(nodes_by_depth, Vector{HypernodeID}())
                # Read solution to compute nodes_by_depth
                println(value.(model[:a]))
            end


        else
            error("Wrong value for parameter ilp.")
        end
    end

    return addergraph
end
