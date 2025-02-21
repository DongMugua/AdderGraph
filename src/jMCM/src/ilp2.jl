#using JuMP

# include("utils_ilp2.jl")

"""
    model_mcm_formulation_2!(model::Model, C::Vector{Int},
                             h::Hypergraph, ad::Int)::Model

Formulation 2 from "M. Kumm -- Multiple Constant Multiplication Optimizations
for Field Programmable Gate Arrays"
"""
function model_mcm_formulation_2!(model::Model, C::Vector{Int},
                                  h::Hypergraph, ad::Int)::Model
    s = 0
    S = Vector{Vector{Int}}()
    for id in 2:h.nodes_number
        value = get_node_value(h, id)
        if value == 1
            push!(S, Vector{Int}())
            s += 1
        end
        push!(S[s], value)
    end

    @variable(model, a[s in 1:ad, S[s]], Bin)
    @variable(model, r[s in 1:ad, S[s]], Bin)
    s = 0
    T = Vector{Dict{Int, Vector{Tuple{Int, Int}}}}()
    for node_id in 2:h.nodes_number
        value = get_node_value(h, node_id)
        if value == 1
            if s != 0
                push!(T, Dict{Int, Vector{Tuple{Int, Int}}}(deepcopy(T[s])))
            else
                push!(T, Dict{Int, Vector{Tuple{Int, Int}}}())
            end
            s += 1
        end
        for arc_id in get_arcs(h, node_id)
            input_nodes = get_input_nodes(h, arc_id)
            node_value_max = get_node_value(h, maximum(input_nodes))
            if !(node_id in input_nodes) && (s == 1 || node_value_max == 1 || node_value_max in keys(T[s-1]))
                push!(get!(T[s], value, Vector{Tuple{Int, Int}}()), (get_node_value(h, input_nodes[1]), node_value_max))
            end
        end
    end
    @variable(model, 0 <= x[s in 1:(ad-1), [Set([S[s][i],S[s][j]]) for i in 1:length(S[s]) for j in i:length(S[s])]] <= 1)

    # C1
    @constraint(model, [w in C], r[ad,w] + a[ad,w] == 1)
    # C2
    @constraint(model, [w in setdiff(S[1],1)], r[1,w] == 0)
    @constraint(model, [s in 2:ad, w in setdiff(S[s],S[s-1])], r[s,w] == 0)
    # C3
    @constraint(model, [s in 2:ad, w in S[s-1]], r[s,w] - a[s-1,w] - r[s-1,w] <= 0)
    # C4
    @constraint(model, [s in 2:ad, w in S[s]], a[s,w] - sum(x[s-1,Set([u,v])] for (u,v) in get(T[s], w, Vector{Tuple{Int, Int}}())) <= 0)
    # C5
    @constraint(model, [s in 1:(ad-1), i in 1:length(S[s]), j in i:length(S[s])], x[s,Set([S[s][i],S[s][j]])] - r[s,S[s][i]] - a[s,S[s][i]] <= 0)
    @constraint(model, [s in 1:(ad-1), i in 1:length(S[s]), j in i:length(S[s])], x[s,Set([S[s][i],S[s][j]])] - r[s,S[s][j]] - a[s,S[s][j]] <= 0)

    # Objective
    @objective(model, Min, sum(sum(a[s,w] for w in S[s]) for s in 1:ad))

    # Stores hypergraph in model as h is necessary to read the solution
    model[:h] = h

    return model
end
