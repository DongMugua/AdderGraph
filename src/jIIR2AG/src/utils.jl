import JuMP.num_constraints
function num_constraints(model::Model)
    return sum(num_constraints(model, F, S) for (F, S) in list_of_constraint_types(model))
end

# function odd(number::Int)
#     if number == 0
#         return 0
#     end
#     while mod(number, 2) == 0
#         number = div(number, 2)
#     end
#     return number
# end

import Base.isempty
function isempty(solution::Iir2AdderGraphs)
    return solution.coefficients_a == (0,0) && solution.coefficients_b == (0,0,0)
end


"""
    fix_to_bounds!(model::Model,
                   bounds_a::Vector{Tuple{Int, Int}},
                   bounds_b::Vector{Tuple{Int, Int}})

Fix variables a_k and b_k to given bounds as in [presolve_a](@ref) and
[presolve_b](@ref).
"""
function fix_to_bounds!(model::Model, bounds_a::Vector{Tuple{Int, Int}}, bounds_b::Vector{Tuple{Int, Int}})
    for m in 1:2
        a_lower_bound, a_upper_bound = bounds_a[m]
        set_lower_bound(model[:a][m], a_lower_bound)
        set_upper_bound(model[:a][m], a_upper_bound)
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
    end

    for m in 0:2
        b_lower_bound, b_upper_bound = bounds_b[m+1]
        set_lower_bound(model[:b][m], b_lower_bound)
        set_upper_bound(model[:b][m], b_upper_bound)
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
    end

    return model
end



"""
    get_specifications(fbands::Vector{Tuple{Float64, Float64}},
                       abands::Vector{Tuple{Float64, Float64}})

Return specifications on a uniform grid from `fbands` and `abands`.
"""
function get_specifications(fbands::Vector{Tuple{Float64, Float64}},
                            abands::Vector{Tuple{Float64, Float64}},
                            size_of_grid::Int=480)
    specifications = Vector{Tuple{Float64, Float64, Float64}}()
    total_band_size = 0
    for band in fbands
        total_band_size += band[2]-band[1]
    end
    for i in 1:length(fbands)
        band = fbands[i]
        size_of_local_grid = ceil(Int, size_of_grid*(band[2]-band[1])/total_band_size)
        for omega in 0:size_of_local_grid
            push!(specifications, (band[1]+(band[2]-band[1])*omega/max(size_of_local_grid,1), abands[i][1], abands[i][2]))
        end
    end
    return specifications
end

get_specifications(fbands::Vector{Tuple{Float64, Float64}},
                   dbands::Union{Vector{Float64},Vector{Int}},
                   delta::Vector{Float64},
                   size_of_grid::Int) = get_specifications(fbands, [(max(0, dbands[i]-delta[i]), dbands[i]+delta[i]) for i in 1:length(fbands)])

get_specifications(fbands::Vector{Tuple{Float64, Float64}},
                   dbands::Union{Vector{Float64},Vector{Int}},
                   delta::Float64,
                   size_of_grid::Int) = get_specifications(fbands, dbands, fill(delta, length(dbands)))


function get_specifications(transfer_function::Function,
                            error_margin_percent::Union{Float64, Int},
                            error_minimum::Float64,
                            size_of_grid::Int)
    specifications::Vector{Tuple{Float64, Float64, Float64}} = [(x/size_of_grid, (transfer_function(x/size_of_grid)*(1-error_margin_percent/100) < error_minimum ? 0.0 : transfer_function(x/size_of_grid)*(1-error_margin_percent/100)), max(transfer_function(x/size_of_grid)*(1+error_margin_percent/100),error_minimum)) for x in 0:size_of_grid if !isnan(transfer_function(x/size_of_grid))]
    return specifications
end




"""
    scale_specifications(specifications_init::Vector{Tuple{Float64, Float64, Float64}},
                         scaling::Tuple{Int, Int})

Scale the specifications by the scaling factors 2^ga and 2^gb where
(ga, gb) = `scaling`.
"""
function scale_specifications(specifications_init::Vector{Tuple{Float64, Float64, Float64}}, scaling::Tuple{Int, Int})
    specifications = Vector{Tuple{Float64, Float64, Float64}}(undef, length(specifications_init))
    g_a = scaling[1]
    g_b = scaling[2]
    for i in 1:length(specifications_init)
        specifications[i] = (specifications_init[i][1], max(0.0, specifications_init[i][2])*2.0^(g_a-g_b), specifications_init[i][3]*2.0^(g_a-g_b))
    end
    return specifications
end


function coefs_fit_specs(
        a1, a2, b0, b1, b2,
        specifications::Vector{Tuple{Float64, Float64, Float64}},
    )
    for (val, Bminus, Bplus) in specifications
        value = abs((b0*exp(2*im*val*pi)+b1*exp(im*val*pi)+b2)/(exp(2*im*val*pi)+a1*exp(im*val*pi)+a2))
        if value > Bplus || value < Bminus
            return false
        end
    end

    return true
end


function is_stable(a1::Float64, a2::Float64)
    return prod(abs.(roots(Polynomial([a2,a1,1]))) .< 1)
end
