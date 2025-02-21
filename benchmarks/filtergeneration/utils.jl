# function odd(number::Int)
#     if number == 0
#         return 0
#     end
#     while mod(number, 2) == 0
#         number = div(number, 2)
#     end
#     return number
# end

function coefs_fit_specs_varepsilon(
        b0, b1, b2, a1, a2,
        specifications::Vector{Tuple{Float64, Float64, Float64}},
    )
    epsilon = 0.0
    for (val, Bminus, Bplus) in specifications
        value = abs((b0*exp(2*im*val*pi)+b1*exp(im*val*pi)+b2)/(exp(2*im*val*pi)+a1*exp(im*val*pi)+a2))
        if value > Bplus || value < Bminus
            epsilon = max(epsilon, value-Bplus, Bminus-value)
        end
    end

    return epsilon
end


function get_specifications(fbands::Vector{Tuple{Float64, Float64}},
                            abands::Vector{Tuple{Float64, Float64}})
    specifications = Vector{Tuple{Float64, Float64, Float64}}()
    for i in 0:48000
        current_band = 0
        for i_bands in 1:length(fbands)
            if i/48000 <= fbands[i_bands][2] && i/48000 >= fbands[i_bands][1]
                current_band = i_bands
                break
            end
        end
        if current_band != 0
            push!(specifications, (i/48000, abands[current_band][1], abands[current_band][2]))
        end
    end
    return specifications
end

get_specifications(fbands::Vector{Tuple{Float64, Float64}},
                   dbands::Union{Vector{Float64},Vector{Int}},
                   delta::Vector{Float64}) = get_specifications(fbands, [(max(0, dbands[i]-delta[i]), dbands[i]+delta[i]) for i in 1:length(fbands)])

get_specifications(fbands::Vector{Tuple{Float64, Float64}},
                   dbands::Union{Vector{Float64},Vector{Int}},
                   delta::Float64) = get_specifications(fbands, dbands, fill(delta, length(dbands)))


function get_specifications(transfer_function::Function,
                            error_margin_percent::Union{Float64, Int},
                            error_minimum::Float64=0.001)
    size_of_grid = 48000
    specifications::Vector{Tuple{Float64, Float64, Float64}} = [(x/size_of_grid, (transfer_function(x/size_of_grid)*(1-error_margin_percent/100) < error_minimum ? 0.0 : transfer_function(x/size_of_grid)*(1-error_margin_percent/100)), max(transfer_function(x/size_of_grid)*(1+error_margin_percent/100),error_minimum)) for x in 0:size_of_grid]
    return specifications
end


function is_stable(a1::Float64, a2::Float64)
    return prod(abs.(roots(Polynomial([a2,a1,1]))) .< 1)
end
