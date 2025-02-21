function main()
    filename = "asic_results.csv"
    lines = Vector{String}()
    open(filename, "r") do f
        lines = readlines(f)
    end

    area_shiftandadd = Vector{Float64}()
    area_truncatedshiftandadd = Vector{Float64}()
    area_truncatedplain = Vector{Float64}()
    power_shiftandadd = Vector{Float64}()
    power_truncatedshiftandadd = Vector{Float64}()
    power_truncatedplain = Vector{Float64}()
    delay_shiftandadd = Vector{Float64}()
    delay_truncatedshiftandadd = Vector{Float64}()
    delay_truncatedplain = Vector{Float64}()
    for line in lines[2:end]
        area, power, delay = parse.(Float64, split(line, ",")[2:4])
        if occursin("_shiftandadd", line)
            push!(area_shiftandadd, area)
            push!(power_shiftandadd, power)
            push!(delay_shiftandadd, delay)
        end
        if occursin("_truncatedshiftandadd", line)
            push!(area_truncatedshiftandadd, area)
            push!(power_truncatedshiftandadd, power)
            push!(delay_truncatedshiftandadd, delay)
        end
        if occursin("_truncatedplain", line)
            push!(area_truncatedplain, area)
            push!(power_truncatedplain, power)
            push!(delay_truncatedplain, delay)
        end
    end
    all_positions = Vector{String}([
    "1",
    "1+1*\\xtickspaces",
    "1+2*\\xtickspaces",
    "1+3*\\xtickspaces",
    "1+4*\\xtickspaces",
    "1+5*\\xtickspaces",
    "1+6*\\xtickspaces+1*\\xtickblockspaces",
    "1+7*\\xtickspaces+1*\\xtickblockspaces",
    "1+8*\\xtickspaces+1*\\xtickblockspaces",
    "1+9*\\xtickspaces+1*\\xtickblockspaces",
    "1+10*\\xtickspaces+2*\\xtickblockspaces",
    "1+11*\\xtickspaces+2*\\xtickblockspaces",
    "1+12*\\xtickspaces+2*\\xtickblockspaces",
    "1+13*\\xtickspaces+2*\\xtickblockspaces",
    "1+14*\\xtickspaces+3*\\xtickblockspaces",
    "1+15*\\xtickspaces+4*\\xtickblockspaces"
    ])
    all_positions_truncated = Vector{String}([
    "1",
    "1+1*\\xtickspaces",
    "1+2*\\xtickspaces",
    "1+3*\\xtickspaces",
    "1+4*\\xtickspaces",
    "1+6*\\xtickspaces+1*\\xtickblockspaces",
    "1+7*\\xtickspaces+1*\\xtickblockspaces",
    "1+8*\\xtickspaces+1*\\xtickblockspaces",
    "1+10*\\xtickspaces+2*\\xtickblockspaces",
    "1+11*\\xtickspaces+2*\\xtickblockspaces",
    "1+12*\\xtickspaces+2*\\xtickblockspaces",
    "1+14*\\xtickspaces+3*\\xtickblockspaces",
    "1+15*\\xtickspaces+4*\\xtickblockspaces"
    ])
    area_string_shiftandadd = "{"
    area_string_truncatedshiftandadd = "{"
    area_string_truncatedplain = "{"
    power_string_shiftandadd = "{"
    power_string_truncatedshiftandadd = "{"
    power_string_truncatedplain = "{"
    delay_string_shiftandadd = "{"
    delay_string_truncatedshiftandadd = "{"
    delay_string_truncatedplain = "{"
    for i in 1:length(all_positions)
        area_string_shiftandadd *= "($(all_positions[i]), $(area_shiftandadd[i]))"
        power_string_shiftandadd *= "($(all_positions[i]), $(power_shiftandadd[i]))"
        delay_string_shiftandadd *= "($(all_positions[i]), $(delay_shiftandadd[i]))"
    end
    for i in 1:length(all_positions_truncated)
        area_string_truncatedshiftandadd *= "($(all_positions_truncated[i]), $(area_truncatedshiftandadd[i]))"
        power_string_truncatedshiftandadd *= "($(all_positions_truncated[i]), $(power_truncatedshiftandadd[i]))"
        delay_string_truncatedshiftandadd *= "($(all_positions_truncated[i]), $(delay_truncatedshiftandadd[i]))"

        area_string_truncatedplain *= "($(all_positions_truncated[i]), $(area_truncatedplain[i]))"
        power_string_truncatedplain *= "($(all_positions_truncated[i]), $(power_truncatedplain[i]))"
        delay_string_truncatedplain *= "($(all_positions_truncated[i]), $(delay_truncatedplain[i]))"
    end
    area_string_shiftandadd *= "};"
    area_string_truncatedshiftandadd *= "};"
    area_string_truncatedplain *= "};"
    power_string_shiftandadd *= "};"
    power_string_truncatedshiftandadd *= "};"
    power_string_truncatedplain *= "};"
    delay_string_shiftandadd *= "};"
    delay_string_truncatedshiftandadd *= "};"
    delay_string_truncatedplain *= "};"
    println(area_string_shiftandadd)
    println(area_string_truncatedshiftandadd)
    println(area_string_truncatedplain)
    println(power_string_shiftandadd)
    println(power_string_truncatedshiftandadd)
    println(power_string_truncatedplain)
    println(delay_string_shiftandadd)
    println(delay_string_truncatedshiftandadd)
    println(delay_string_truncatedplain)

    return nothing
end
