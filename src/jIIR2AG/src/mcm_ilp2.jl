function mcm_ilp_addergraph!(model::Model, C::Tuple{Vector{VariableRef}, Vector{VariableRef}}, AG::CompleteAdderGraph)
    max_depth = maximum(keys(AG.depth_adders_relation))
    S = Vector{Vector{Int}}()
    for depth in 1:max_depth
        push!(S, Vector{Int}())
        S[depth] = collect(get_values_at_depth(AG, depth))
    end

    T = Vector{Dict{Int, Vector{Tuple{Int, Int}}}}()
    for depth in 1:max_depth
        push!(T, Dict{Int, Vector{Tuple{Int, Int}}}())
        for adder_id in get_adder_id_at_depth(AG, depth)
            for edge_id in get_edges(AG, adder_id)
                input_adders = get_input_adders(AG, edge_id)
                if !(adder_id in input_adders)
                    value_1 = min(get_adder_value(AG, input_adders[1]), get_adder_value(AG, input_adders[2]))
                    value_2 = max(get_adder_value(AG, input_adders[1]), get_adder_value(AG, input_adders[2]))
                    push!(get!(T[depth], get_adder_value(AG, adder_id), Vector{Tuple{Int, Int}}()), (value_1, value_2))
                end
            end
        end
    end

    @variable(model, a_a[s in 1:max_depth, S[s]], Bin)
    @variable(model, r_a[s in 1:max_depth, S[s]], Bin)
    @variable(model, a_b[s in 1:max_depth, S[s]], Bin)
    @variable(model, r_b[s in 1:max_depth, S[s]], Bin)

    @variable(model, 0 <= x_a[s in 1:(max_depth-1), [Set([S[s][i],S[s][j]]) for i in 1:length(S[s]) for j in i:length(S[s])]] <= 1)
    @variable(model, 0 <= x_b[s in 1:(max_depth-1), [Set([S[s][i],S[s][j]]) for i in 1:length(S[s]) for j in i:length(S[s])]] <= 1)

    @variable(model, a_value[1:2, 0:AG.maximum_value], Bin)
    @variable(model, b_value[0:2, 0:AG.maximum_value], Bin)

    # C3a / C3b
    @constraint(model, [m in 1:2], C[1][m] == sum(w*a_value[m, w] for w in 0:AG.maximum_value))
    @constraint(model, [m in 1:2], sum(a_value[m, w] for w in 0:AG.maximum_value) == 1)
    @constraint(model, [m in 0:2], C[2][m+1] == sum(w*b_value[m, w] for w in 0:AG.maximum_value))
    @constraint(model, [m in 0:2], sum(b_value[m, w] for w in 0:AG.maximum_value) == 1)
    # C4
    @constraint(model, [w in 1:AG.maximum_value], 2*r_a[max_depth,odd(w)] + 2*a_a[max_depth,odd(w)] >= sum(a_value[m, w] for m in 1:2))
    @constraint(model, [w in 1:AG.maximum_value], 2*r_b[max_depth,odd(w)] + 2*a_b[max_depth,odd(w)] >= sum(b_value[m, w] for m in 0:2))
    # C5
    @constraint(model, [w in setdiff(S[1],1)], r_a[1,w] == 0)
    @constraint(model, [w in setdiff(S[1],1)], r_b[1,w] == 0)
    @constraint(model, [s in 2:max_depth, w in setdiff(S[s],S[s-1])], r_a[s,w] == 0)
    @constraint(model, [s in 2:max_depth, w in setdiff(S[s],S[s-1])], r_b[s,w] == 0)
    # C6
    @constraint(model, [s in 2:max_depth, w in S[s-1]], r_a[s,w] - a_a[s-1,w] - r_a[s-1,w] <= 0)
    @constraint(model, [s in 2:max_depth, w in S[s-1]], r_b[s,w] - a_b[s-1,w] - r_b[s-1,w] <= 0)
    # C7
    @constraint(model, [s in 2:max_depth, w in S[s]], a_a[s,w] - sum(x_a[s-1,Set([u,v])] for (u,v) in get(T[s], w, Vector{Tuple{Int, Int}}())) <= 0)
    @constraint(model, [s in 2:max_depth, w in S[s]], a_b[s,w] - sum(x_b[s-1,Set([u,v])] for (u,v) in get(T[s], w, Vector{Tuple{Int, Int}}())) <= 0)
    # C8
    @constraint(model, [s in 1:(max_depth-1), i in 1:length(S[s]), j in i:length(S[s])], x_a[s,Set([S[s][i],S[s][j]])] - r_a[s,S[s][i]] - a_a[s,S[s][i]] <= 0)
    @constraint(model, [s in 1:(max_depth-1), i in 1:length(S[s]), j in i:length(S[s])], x_a[s,Set([S[s][i],S[s][j]])] - r_a[s,S[s][j]] - a_a[s,S[s][j]] <= 0)
    @constraint(model, [s in 1:(max_depth-1), i in 1:length(S[s]), j in i:length(S[s])], x_b[s,Set([S[s][i],S[s][j]])] - r_b[s,S[s][i]] - a_b[s,S[s][i]] <= 0)
    @constraint(model, [s in 1:(max_depth-1), i in 1:length(S[s]), j in i:length(S[s])], x_b[s,Set([S[s][i],S[s][j]])] - r_b[s,S[s][j]] - a_b[s,S[s][j]] <= 0)

    # Objective
    @objective(model, Min, sum(sum(a_a[s,w]+a_b[s,w] for w in S[s]) for s in 1:max_depth)-a_value[1,0]-a_value[2,0]-b_value[0,0]-b_value[1,0]-b_value[2,0])

    return model
end


function mcm_ilp_2!(model::Model, C::Tuple{Vector{VariableRef}, Vector{VariableRef}}, wordlength::Int, max_depth::Int)
    @time AG = construct_addergraph(wordlength, max_depth)
    mcm_ilp_addergraph!(model, C, AG)
    unset_silent(model)
    optimize!(model)
    S = Vector{Vector{Int}}()
    for depth in 1:max_depth
        push!(S, Vector{Int}())
        S[depth] = collect(get_values_at_depth(AG, depth))
    end
    model[:NA] = sum(sum(round(Int, value(model[:a_a][s,w]))+round(Int, value(model[:a_b][s,w])) for w in S[s]) for s in 1:max_depth)

    return model
end
