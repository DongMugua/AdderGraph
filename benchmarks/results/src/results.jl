include("plots.jl")


function percentbetterdelayandluts(
        filename::String;
        acceptable_error::Union{Float64, Int}=-1,
        wordlengths_range::Tuple{Int, Int}=(4,10),
    )
    lines = Vector{String}()

    open(filename, "r") do f
        lines = readlines(f)
    end

    dict_luts_shiftandadd = Dict{String, Int}()
    dict_luts_truncatedplain = Dict{String, Int}()
    dict_luts_truncatedplain0dsp = Dict{String, Int}()
    dict_luts_truncatedshiftandadd = Dict{String, Int}()
    dict_dsps_truncatedplain = Dict{String, Int}()
    dict_times_shiftandadd = Dict{String, Float64}()
    dict_times_truncatedplain = Dict{String, Float64}()
    dict_times_truncatedplain0dsp = Dict{String, Float64}()
    dict_times_truncatedshiftandadd = Dict{String, Float64}()
    dict_eps_shiftandadd = Dict{String, Float64}()
    dict_eps_truncatedplain = Dict{String, Float64}()
    dict_eps_truncatedplain0dsp = Dict{String, Float64}()
    dict_eps_truncatedshiftandadd = Dict{String, Float64}()

    dict_eps_kcm = Dict{String, Float64}()
    dict_kcmlut = Dict{String, Int}()
    dict_kcmtime = Dict{String, Float64}()
    list_filtername = Vector{String}()
    list_filters = Vector{String}()

    for line in lines
        if occursin("fixiir", line)
            filtername = split(line, ";")[1]
            push!(list_filtername, filtername)
            push!(list_filters, split(line, "_")[1])
            dict_kcmlut[filtername] = parse(Int, split(line, ";")[2])
            dict_kcmtime[filtername] = parse(Float64, split(line, ";")[4][1:(end-2)])
            dict_eps_kcm[filtername] = parse(Float64, split(line, ";")[5])
        else
            filtername = split(line, ";")[1]
            push!(list_filtername, filtername)
            push!(list_filters, split(line, "_")[1])
            if occursin("_shiftandadd", line)
                dict_luts_shiftandadd[filtername] = parse(Int, split(line, ";")[2])
                dict_times_shiftandadd[filtername] = parse(Float64, split(line, ";")[4][1:(end-2)])
                dict_eps_shiftandadd[filtername] = parse(Float64, split(line, ";")[5])
            elseif occursin("_truncatedplain", line) && !occursin("_truncatedplain0dsp", line)
                dict_luts_truncatedplain[filtername] = parse(Int, split(line, ";")[2])
                dict_times_truncatedplain[filtername] = parse(Float64, split(line, ";")[4][1:(end-2)])
                dict_dsps_truncatedplain[filtername] = parse(Int, split(line, ";")[3])
                dict_eps_truncatedplain[filtername] = parse(Float64, split(line, ";")[5])
            elseif occursin("_truncatedplain0dsp", line)
                dict_luts_truncatedplain0dsp[filtername] = parse(Int, split(line, ";")[2])
                dict_times_truncatedplain0dsp[filtername] = parse(Float64, split(line, ";")[4][1:(end-2)])
                dict_eps_truncatedplain0dsp[filtername] = parse(Float64, split(line, ";")[5])
            elseif occursin("_truncatedshiftandadd", line)
                dict_luts_truncatedshiftandadd[filtername] = parse(Int, split(line, ";")[2])
                dict_times_truncatedshiftandadd[filtername] = parse(Float64, split(line, ";")[4][1:(end-2)])
                dict_eps_truncatedshiftandadd[filtername] = parse(Float64, split(line, ";")[5])
            else
                error(line)
            end
        end
    end
    unique!(list_filters)

    speedup_percent_fixiir = Vector{Float64}()
    speedup_percent_truncatedplain0dsp = Vector{Float64}()
    speedup_percent_truncatedplain = Vector{Float64}()
    speedup_percent_truncatedshiftandadd = Vector{Float64}()
    luts_decrease_percent_fixiir = Vector{Float64}()
    luts_decrease_percent_truncatedplain0dsp = Vector{Float64}()
    luts_decrease_percent_truncatedplain = Vector{Float64}()
    luts_decrease_percent_truncatedshiftandadd = Vector{Float64}()
    for filtername in list_filtername
        if occursin("_fixiir", filtername)
            for wordlength in wordlengths_range[1]:wordlengths_range[2]
                if (acceptable_error < 0 || dict_eps_kcm[filtername] <= acceptable_error) && get(dict_times_shiftandadd, replace(filtername, "_cw00_fixiir" => "_cw0$(wordlength)_shiftandadd"), nothing) != nothing
                    push!(speedup_percent_fixiir, (dict_kcmtime[filtername]-dict_times_shiftandadd[replace(filtername, "_cw00_fixiir" => "_cw0$(wordlength)_shiftandadd")])*100/dict_kcmtime[filtername])
                    push!(luts_decrease_percent_fixiir, (dict_kcmlut[filtername]-dict_luts_shiftandadd[replace(filtername, "_cw00_fixiir" => "_cw0$(wordlength)_shiftandadd")])*100/dict_kcmlut[filtername])
                    break
                end
            end
        elseif occursin("_truncatedplain", filtername) && !occursin("_truncatedplain0dsp", filtername)
            if (acceptable_error < 0 || dict_eps_truncatedplain[filtername] <= acceptable_error) && get(dict_times_shiftandadd, replace(filtername, "_truncatedplain" => "_shiftandadd"), nothing) != nothing
                push!(speedup_percent_truncatedplain, (dict_times_truncatedplain[filtername]-dict_times_shiftandadd[replace(filtername, "_truncatedplain" => "_shiftandadd")])*100/dict_times_truncatedplain[filtername])
                push!(luts_decrease_percent_truncatedplain, (dict_luts_truncatedplain[filtername]-dict_luts_shiftandadd[replace(filtername, "_truncatedplain" => "_shiftandadd")])*100/dict_luts_truncatedplain[filtername])
            end
        elseif occursin("_truncatedplain0dsp", filtername)
            if (acceptable_error < 0 || dict_eps_truncatedplain0dsp[filtername] <= acceptable_error) && get(dict_times_shiftandadd, replace(filtername, "_truncatedplain0dsp" => "_shiftandadd"), nothing) != nothing
                push!(speedup_percent_truncatedplain0dsp, (dict_times_truncatedplain0dsp[filtername]-dict_times_shiftandadd[replace(filtername, "_truncatedplain0dsp" => "_shiftandadd")])*100/dict_times_truncatedplain0dsp[filtername])
                push!(luts_decrease_percent_truncatedplain0dsp, (dict_luts_truncatedplain0dsp[filtername]-dict_luts_shiftandadd[replace(filtername, "_truncatedplain0dsp" => "_shiftandadd")])*100/dict_luts_truncatedplain0dsp[filtername])
            end
        elseif occursin("_truncatedshiftandadd", filtername)
            if (acceptable_error < 0 || dict_eps_truncatedshiftandadd[filtername] <= acceptable_error) && get(dict_times_shiftandadd, replace(filtername, "_truncatedshiftandadd" => "_shiftandadd"), nothing) != nothing
                push!(speedup_percent_truncatedshiftandadd, (dict_times_truncatedshiftandadd[filtername]-dict_times_shiftandadd[replace(filtername, "_truncatedshiftandadd" => "_shiftandadd")])*100/dict_times_truncatedshiftandadd[filtername])
                push!(luts_decrease_percent_truncatedshiftandadd, (dict_luts_truncatedshiftandadd[filtername]-dict_luts_shiftandadd[replace(filtername, "_truncatedshiftandadd" => "_shiftandadd")])*100/dict_luts_truncatedshiftandadd[filtername])
            end
        end
    end

    return ((sum(speedup_percent_fixiir)/length(speedup_percent_fixiir),sum(speedup_percent_truncatedplain)/length(speedup_percent_truncatedplain),sum(speedup_percent_truncatedplain0dsp)/length(speedup_percent_truncatedplain0dsp),sum(speedup_percent_truncatedshiftandadd)/length(speedup_percent_truncatedshiftandadd)),
    (sum(luts_decrease_percent_fixiir)/length(luts_decrease_percent_fixiir),sum(luts_decrease_percent_truncatedplain)/length(luts_decrease_percent_truncatedplain),sum(luts_decrease_percent_truncatedplain0dsp)/length(luts_decrease_percent_truncatedplain0dsp),sum(luts_decrease_percent_truncatedshiftandadd)/length(luts_decrease_percent_truncatedshiftandadd)))
end


function results()
    filename = (@__DIR__)*"/../all_results.csv"
    println(percentbetterdelayandluts(filename))
    main_plots(filename)
    return nothing
end
