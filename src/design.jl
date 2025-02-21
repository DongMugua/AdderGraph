function design2vhdl(
        model::Model,
        specifications::Vector{Tuple{Float64, Float64, Float64}},
        wordlength::Int;
        instance_name::String = "",
        design_iir_parameters::Vector{Pair} = Vector{Pair}(),
        vhdl_filename::String="",
        verbose::Bool=false, debug::Bool=true,
    )
    if isempty(vhdl_filename)
        vhdl_filename = tempname()
        @warn "No filename provided, usage: design2vhdl(..., vhdl_filename=\"filename.vhd\")\nFile stored at $(vhdl_filename)"
    end

    solution = design_second_order_iir(model, specifications, wordlength, design_iir_parameters..., verbose=verbose)
    # has_values(model) && println(value.(model[:a0]))
    # has_values(model) && println(value.(model[:gab]))

    if !isempty(solution)
        wcpg_val = wcpg(solution)
        wcpg_eps_val = wcpg_eps(solution)
        lsbin = -16
        lsbout = -16
        write_vhdl(vhdl_filename, solution, wcpg_val=wcpg_val, wcpg_eps_val=wcpg_eps_val, lsbin=lsbin, lsbout=lsbout)

    end

    return nothing
end

design2vhdl(optimizer::DataType, args...; kwargs...) = design2vhdl(Model(optimizer), args...; kwargs...)

design2vhdl(model::Model,
            fbands::Vector{Tuple{Float64, Float64}},
            abands::Vector{Tuple{Float64, Float64}},
            wordlength::Int;
            size_of_grid::Int=480,
            kwargs...
    ) = design2vhdl(model, get_specifications(fbands, abands, size_of_grid), wordlength; kwargs...)

design2vhdl(model::Model,
            fbands::Vector{Tuple{Float64, Float64}},
            dbands::Union{Vector{Float64},Vector{Int}},
            deltas::Vector{Float64},
            wordlength::Int;
            size_of_grid::Int=480,
            kwargs...
    ) = design2vhdl(model, get_specifications(fbands, dbands, deltas, size_of_grid), wordlength; kwargs...)

design2vhdl(model::Model,
            fbands::Vector{Tuple{Float64, Float64}},
            dbands::Union{Vector{Float64},Vector{Int}},
            delta::Float64,
            wordlength::Int;
            size_of_grid::Int=480,
            kwargs...
    ) = design2vhdl(model, get_specifications(fbands, dbands, delta, size_of_grid), wordlength; kwargs...)

design2vhdl(model::Model,
            transfer_function::Function,
            error_margin_percent::Union{Float64, Int},
            error_minimum::Float64,
            wordlength::Int;
            size_of_grid::Int=480,
            kwargs...
    ) = design2vhdl(model, get_specifications(transfer_function, error_margin_percent, error_minimum, size_of_grid), wordlength; kwargs...)
