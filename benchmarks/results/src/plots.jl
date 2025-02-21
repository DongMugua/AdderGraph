# function odd(number::Int)
#     if number == 0
#         return 0
#     end
#     while mod(number, 2) == 0
#         number = div(number, 2)
#     end
#     return number
# end

function plots_for_each_filter(
        filename::String;
        acceptable_error::Union{Float64, Int}=-1,
        with_eps::Bool=false,
        wordlengths_range::Tuple{Int, Int}=(4,10),
        color_fixiir::String="blue",
        color_shiftandadd::String="green",
        color_truncatedshiftandadd::String="yellow",
        color_truncatedplain::String="red",
        color_truncatedplain0dsp::String="orange",
        fill_fixiir::String="blue!30!white",
        fill_shiftandadd::String="green!30!white",
        fill_truncatedshiftandadd::String="yellow!30!white",
        fill_truncatedplain::String="red!30!white",
        fill_truncatedplain0dsp::String="orange!30!white",
        name_fixiir::String="fixiir",
        name_shiftandadd::String="shiftandadd",
        name_truncatedplain::String="truncatedplain",
        name_truncatedplain0dsp::String="truncatedplain0dsp",
        name_truncatedshiftandadd::String="truncatedshiftandadd",
        luts_or_times::Int=1, #1 luts, 2 times, 3 both
        printnbadders::Bool=true,
    )
    lines = Vector{String}()
    open(filename, "r") do f
        lines = readlines(f)
    end

    wordlength_str_tick = string(collect((wordlengths_range[1]):(wordlengths_range[2])).-3)[2:(end-1)]
    wordlength_str = string(collect((wordlengths_range[1]):(wordlengths_range[2])))[2:(end-1)]

    figure_header_str = """
    \\begin{figure}$(luts_or_times == 3 ? "\n\\vspace{-3cm}" : "")
    \\centering
    """
    tikz_header_luts_str = """
    \\begin{tikzpicture}
        \\begin{axis}[
                ybar,
                bar width=0.15,
                legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},$(luts_or_times == 3 ? "\n            every axis legend/.code={\\let\\addlegendentry\\relax}," : "")
                ymin=0,
                ytick pos=left,
                xtick pos=bottom,
                width=\\textwidth,
                ylabel={\\#LUTs},
                xtick={0, $(wordlength_str_tick)},
                xticklabels={KCM, $(wordlength_str)},
                point meta=explicit symbolic
            ]
    """
    tikz_header_times_str = """
    \\begin{tikzpicture}
        \\begin{axis}[
                ybar,
                bar width=0.15,
                legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},
                ymin=0,
                ytick pos=left,
                xtick pos=bottom,
                width=\\textwidth,
                ylabel={Time in ns},
                xtick={0, $(wordlength_str_tick)},
                xticklabels={KCM, $(wordlength_str)},
                point meta=explicit symbolic
            ]
    """

    eps_axis_str = """
        \\end{axis}
        \\begin{axis}[
                %ymode=log,
                axis lines*=right,
                ymin=0.0,
                ytick pos=right,
                xmajorticks=false,
                width=\\textwidth,
                ylabel={error},
                xtick={0, $(wordlength_str_tick)},
                xticklabels={KCM, $(wordlength_str)},
                point meta=explicit symbolic,
            ]
    """

    tikz_footer_luts_str = tikz_footer_times_str = """
        \\end{axis}
    \\end{tikzpicture}
    """

    function figure_footer_str(caption_str)
        return "\\caption{$caption_str}\n\\end{figure}\n\n\\clearpage\n\n"
    end

    plots_str = ""

    luts_shiftandadd = zeros(Int, wordlengths_range[2])
    times_shiftandadd = zeros(Float64, wordlengths_range[2])
    eps_shiftandadd = zeros(Float64, wordlengths_range[2])
    luts_truncatedplain = zeros(Int, wordlengths_range[2])
    times_truncatedplain = zeros(Float64, wordlengths_range[2])
    eps_truncatedplain = zeros(Float64, wordlengths_range[2])
    luts_truncatedplain0dsp = zeros(Int, wordlengths_range[2])
    times_truncatedplain0dsp = zeros(Float64, wordlengths_range[2])
    eps_truncatedplain0dsp = zeros(Float64, wordlengths_range[2])
    luts_truncatedshiftandadd = zeros(Int, wordlengths_range[2])
    times_truncatedshiftandadd = zeros(Float64, wordlengths_range[2])
    eps_truncatedshiftandadd = zeros(Float64, wordlengths_range[2])
    dsps_truncatedplain = zeros(Int, wordlengths_range[2])

    adders_shiftandadd = zeros(Int, wordlengths_range[2])
    adders_truncatedshiftandadd = zeros(Int, wordlengths_range[2])

    line = lines[1]
    kcmlut = parse(Int, split(line, ";")[2])
    kcmtime = parse(Float64, split(line, ";")[4][1:(end-2)])
    eps_kcm = parse(Float64, split(line, ";")[5])
    filtername = split(line, "_")[1]*"\\_"*split(line, "_")[2]
    push!(lines, line) # To output the last filter

    for line in lines[2:end]
        if occursin("fixiir", line) || filtername != split(line, "_")[1]*"\\_"*split(line, "_")[2]
            initialize_xlabels_str = ""
            initialize_xlabels_eps_str = ""
            luts_shiftandadd_str = ""
            lutsdsp_truncatedplain_str = ""
            luts_truncatedplain0dsp_str = ""
            luts_truncatedshiftandadd_str = ""
            times_shiftandadd_str = ""
            times_truncatedplain_str = ""
            times_truncatedplain0dsp_str = ""
            times_truncatedshiftandadd_str = ""
            eps_shiftandadd_str = ""
            eps_truncatedplain_str = ""
            eps_truncatedplain0dsp_str = ""
            eps_truncatedshiftandadd_str = ""
            for wordlength in wordlengths_range[1]:wordlengths_range[2]
                initialize_xlabels_str *= "($(wordlength-3), 1) "
                initialize_xlabels_eps_str *= "($(wordlength-3), 0) "
                if luts_shiftandadd[wordlength] != 0
                    eps_shiftandadd_str *= "($(wordlength-3), $(eps_shiftandadd[wordlength])) "
                    if acceptable_error < 0 || eps_shiftandadd[wordlength] <= acceptable_error
                        luts_shiftandadd_str *= "($(wordlength-3), $(luts_shiftandadd[wordlength]))"
                        times_shiftandadd_str *= "($(wordlength-3), $(times_shiftandadd[wordlength]))"
                        if printnbadders
                            luts_shiftandadd_str *= "[\$$(adders_shiftandadd[wordlength])\$]"
                            times_shiftandadd_str *= "[\$$(adders_shiftandadd[wordlength])\$]"
                        end
                        luts_shiftandadd_str *= "\n            "
                        times_shiftandadd_str *= "\n            "
                    end
                end
                eps_truncatedplain0dsp_str *= "($(wordlength-3), $(eps_truncatedplain0dsp[wordlength])) "
                if acceptable_error < 0 || eps_truncatedplain0dsp[wordlength] <= acceptable_error
                    luts_truncatedplain0dsp_str *= "($(wordlength-3), $(luts_truncatedplain0dsp[wordlength]))\n            "
                    times_truncatedplain0dsp_str *= "($(wordlength-3), $(times_truncatedplain0dsp[wordlength]))\n            "
                end
                eps_truncatedplain_str *= "($(wordlength-3), $(eps_truncatedplain[wordlength])) "
                if acceptable_error < 0 || eps_truncatedplain[wordlength] <= acceptable_error
                    lutsdsp_truncatedplain_str *= "($(wordlength-3), $(luts_truncatedplain[wordlength]))"
                    times_truncatedplain_str *= "($(wordlength-3), $(times_truncatedplain[wordlength]))"
                    if dsps_truncatedplain[wordlength] != 0
                        lutsdsp_truncatedplain_str *= "[\$+$(dsps_truncatedplain[wordlength])\$ DSPs]"
                    end
                    lutsdsp_truncatedplain_str *= "\n            "
                end
                eps_truncatedshiftandadd_str *= "($(wordlength-3), $(eps_truncatedshiftandadd[wordlength])) "
                if acceptable_error < 0 || eps_truncatedshiftandadd[wordlength] <= acceptable_error
                    luts_truncatedshiftandadd_str *= "($(wordlength-3), $(luts_truncatedshiftandadd[wordlength]))"
                    times_truncatedshiftandadd_str *= "($(wordlength-3), $(times_truncatedshiftandadd[wordlength]))"
                    if printnbadders
                        luts_truncatedshiftandadd_str *= "[\$$(adders_truncatedshiftandadd[wordlength])\$]"
                        times_truncatedshiftandadd_str *= "[\$$(adders_truncatedshiftandadd[wordlength])\$]"
                    end
                    luts_truncatedshiftandadd_str *= "\n            "
                    times_truncatedshiftandadd_str *= "\n            "
                end
            end

            plots_str *= figure_header_str
            if luts_or_times in [1,3]
                plots_str *= tikz_header_luts_str
                plots_str *= "        \\addplot [white, only marks, draw opacity=0, fill opacity=0] coordinates {(0, 1) $(initialize_xlabels_str)};\n"
                if (acceptable_error < 0 || eps_kcm <= acceptable_error) && kcmlut != 0
                    plots_str *= "        \\addplot[legend entry=$(name_fixiir), $(color_fixiir), fill=$(fill_fixiir)] coordinates {(0, $kcmlut)};\n"
                end
                if luts_shiftandadd_str != ""
                    plots_str *= "        \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_shiftandadd), $(color_shiftandadd), fill=$(fill_shiftandadd)] coordinates {\n            $luts_shiftandadd_str};\n"
                end
                if lutsdsp_truncatedplain_str != ""
                    plots_str *= "        \\addplot[every node near coord/.append style={font=\\footnotesize, rotate=90, anchor=west, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_truncatedplain), $(color_truncatedplain), fill=$(fill_truncatedplain)] coordinates {\n            $lutsdsp_truncatedplain_str};\n"
                end
                if luts_truncatedplain0dsp_str != ""
                    plots_str *= "        \\addplot[legend entry=$(name_truncatedplain0dsp), $(color_truncatedplain0dsp), fill=$(fill_truncatedplain0dsp)] coordinates {\n            $luts_truncatedplain0dsp_str};\n"
                end
                if luts_truncatedshiftandadd_str != ""
                    plots_str *= "        \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_truncatedshiftandadd), $(color_truncatedshiftandadd), fill=$(fill_truncatedshiftandadd)] coordinates {\n            $luts_truncatedshiftandadd_str};\n"
                end

                if with_eps
                    plots_str *= eps_axis_str
                    plots_str *= "        \\addplot[draw=white,ultra thick,draw opacity=0] coordinates {(0, 0) $initialize_xlabels_eps_str};\n"
                    plots_str *= """
                            %\\addplot[draw=$(color_fixiir),ultra thick,mark=*,mark options={solid,draw=$(color_fixiir),fill=$(fill_fixiir)}] coordinates {(0, $eps_kcm)};
                            %\\addplot[draw=$(color_shiftandadd),ultra thick,mark=*,mark options={solid,draw=$(color_shiftandadd),fill=$(fill_shiftandadd)}] coordinates {$eps_shiftandadd_str};
                            \\addplot[draw=$(color_truncatedshiftandadd)!70!$(color_truncatedplain),ultra thick,mark=*,mark options={solid,draw=$(color_truncatedshiftandadd)!70!$(color_truncatedplain),fill=$(fill_truncatedshiftandadd)!70!$(fill_truncatedplain)}] coordinates {$eps_truncatedshiftandadd_str};
                    """
                end

                plots_str *= tikz_footer_luts_str
            end

            if luts_or_times in [2,3]
                plots_str *= tikz_header_times_str
                plots_str *= "        \\addplot [white, only marks, draw opacity=0, fill opacity=0] coordinates {(0, 1) $(initialize_xlabels_str)};\n"
                if (acceptable_error < 0 || eps_kcm <= acceptable_error) && kcmlut != 0
                    plots_str *= "        \\addplot[legend entry=$(name_fixiir), $(color_fixiir), fill=$(fill_fixiir)] coordinates {(0, $kcmtime)};\n"
                end
                if times_shiftandadd_str != ""
                    plots_str *= "        \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_shiftandadd), $(color_shiftandadd), fill=$(fill_shiftandadd)] coordinates {\n            $times_shiftandadd_str};\n"
                end
                if times_truncatedplain_str != ""
                    plots_str *= "        \\addplot[legend entry=$(name_truncatedplain), $(color_truncatedplain), fill=$(fill_truncatedplain)] coordinates {\n            $times_truncatedplain_str};\n"
                end
                if times_truncatedplain0dsp_str != ""
                    plots_str *= "        \\addplot[legend entry=$(name_truncatedplain0dsp), $(color_truncatedplain0dsp), fill=$(fill_truncatedplain0dsp)] coordinates {\n            $times_truncatedplain0dsp_str};\n"
                end
                if times_truncatedshiftandadd_str != ""
                    plots_str *= "        \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_truncatedshiftandadd), $(color_truncatedshiftandadd), fill=$(fill_truncatedshiftandadd)] coordinates {\n            $times_truncatedshiftandadd_str};\n"
                end

                if with_eps
                    plots_str *= eps_axis_str
                    plots_str *= "        \\addplot[draw=white,ultra thick,draw opacity=0] coordinates {(0, 0) $initialize_xlabels_eps_str};\n"
                    plots_str *= """
                            %\\addplot[draw=$(color_fixiir),ultra thick,mark=*,mark options={solid,draw=$(color_fixiir),fill=$(fill_fixiir)}] coordinates {(0, $eps_kcm)};
                            %\\addplot[draw=$(color_shiftandadd),ultra thick,mark=*,mark options={solid,draw=$(color_shiftandadd),fill=$(fill_shiftandadd)}] coordinates {$eps_shiftandadd_str};
                            \\addplot[draw=$(color_truncatedshiftandadd)!70!$(color_truncatedplain),ultra thick,mark=*,mark options={solid,draw=$(color_truncatedshiftandadd)!70!$(color_truncatedplain),fill=$(fill_truncatedshiftandadd)!70!$(fill_truncatedplain)}] coordinates {$eps_truncatedshiftandadd_str};
                    """
                end

                plots_str *= tikz_footer_times_str
            end

            plots_str *= figure_footer_str(filtername)

            luts_shiftandadd = zeros(Int, wordlengths_range[2])
            times_shiftandadd = zeros(Float64, wordlengths_range[2])
            eps_shiftandadd = zeros(Float64, wordlengths_range[2])
            luts_truncatedplain = zeros(Int, wordlengths_range[2])
            times_truncatedplain = zeros(Float64, wordlengths_range[2])
            eps_truncatedplain = zeros(Float64, wordlengths_range[2])
            luts_truncatedplain0dsp = zeros(Int, wordlengths_range[2])
            times_truncatedplain0dsp = zeros(Float64, wordlengths_range[2])
            eps_truncatedplain0dsp = zeros(Float64, wordlengths_range[2])
            luts_truncatedshiftandadd = zeros(Int, wordlengths_range[2])
            times_truncatedshiftandadd = zeros(Float64, wordlengths_range[2])
            eps_truncatedshiftandadd = zeros(Float64, wordlengths_range[2])
            dsps_truncatedplain = zeros(Int, wordlengths_range[2])
            kcmlut = 0
            kcmtime = 0.0
            eps_kcm = 0.0
            if occursin("_fixiir", line)
                kcmlut = parse(Int, split(line, ";")[2])
                kcmtime = parse(Float64, split(line, ";")[4][1:(end-2)])
                eps_kcm = parse(Float64, split(line, ";")[5])
            end
            filtername = split(line, "_")[1]*"\\_"*split(line, "_")[2]
            adders_shiftandadd = zeros(Int, wordlengths_range[2])
            adders_truncatedshiftandadd = zeros(Int, wordlengths_range[2])
        else
            wordlength = parse(Int, (split(line, "_")[3])[3:4])
            if wordlength <= wordlengths_range[2]
                if occursin("_shiftandadd", line)
                    luts_shiftandadd[wordlength] = parse(Int, split(line, ";")[2])
                    times_shiftandadd[wordlength] = parse(Float64, split(line, ";")[4][1:(end-2)])
                    eps_shiftandadd[wordlength] = parse(Float64, split(line, ";")[5])
                    adders_shiftandadd[wordlength] = parse(Int, split(line, ";")[6])
                elseif occursin("_truncatedplain", line) && !occursin("_truncatedplain0dsp", line)
                    luts_truncatedplain[wordlength] = parse(Int, split(line, ";")[2])
                    dsps_truncatedplain[wordlength] = parse(Int, split(line, ";")[3])
                    times_truncatedplain[wordlength] = parse(Float64, split(line, ";")[4][1:(end-2)])
                    eps_truncatedplain[wordlength] = parse(Float64, split(line, ";")[5])
                elseif occursin("_truncatedplain0dsp", line)
                    luts_truncatedplain0dsp[wordlength] = parse(Int, split(line, ";")[2])
                    times_truncatedplain0dsp[wordlength] = parse(Float64, split(line, ";")[4][1:(end-2)])
                    eps_truncatedplain0dsp[wordlength] = parse(Float64, split(line, ";")[5])
                elseif occursin("_truncatedshiftandadd", line)
                    luts_truncatedshiftandadd[wordlength] = parse(Int, split(line, ";")[2])
                    times_truncatedshiftandadd[wordlength] = parse(Float64, split(line, ";")[4][1:(end-2)])
                    eps_truncatedshiftandadd[wordlength] = parse(Float64, split(line, ";")[5])
                    adders_truncatedshiftandadd[wordlength] = parse(Int, split(line, ";")[6])
                else
                    error(line)
                end
            end
        end
    end

    return plots_str
end






function plots_fixiir(
        filename::String;
        acceptable_error::Union{Float64, Int}=-1,
        with_eps::Bool=false,
        wordlengths_range::Tuple{Int, Int}=(4,10),
        color_fixiir::String="blue",
        color_shiftandadd::String="green",
        color_truncatedshiftandadd::String="yellow",
        color_truncatedplain::String="red",
        color_truncatedplain0dsp::String="orange",
        fill_fixiir::String="blue!30!white",
        fill_shiftandadd::String="green!30!white",
        fill_truncatedshiftandadd::String="yellow!30!white",
        fill_truncatedplain::String="red!30!white",
        fill_truncatedplain0dsp::String="orange!30!white",
        name_fixiir::String="fixiir",
        name_shiftandadd::String="shiftandadd",
        name_truncatedplain::String="truncatedplain",
        name_truncatedplain0dsp::String="truncatedplain0dsp",
        name_truncatedshiftandadd::String="truncatedshiftandadd",
        luts_or_times::Int=1, #1 luts, 2 times, 3 both
    )
    lines = Vector{String}()

    open(filename, "r") do f
        lines = readlines(f)
    end

    filters_dict = Dict{String, Int}([
        "lp1x0" => 1,
        "lp1x1" => 2,
        "lp1x2" => 3,
        "lp1x3" => 4,
        "lp1x4" => 5,
        "lp1x5" => 6,
        "lp2x0" => 8,
        "lp2x1" => 9,
        "lp2x2" => 10,
        "lp2x3" => 11,
        "lp3x0" => 13,
        "lp3x1" => 14,
        "lp3x2" => 15,
        "lp3x3" => 16,
        "lp4" => 18,
        "hp0" => 20,
    ])

    dict_eps_kcm = Dict{String, Float64}()
    dict_kcmlut = Dict{String, Int}()
    dict_kcmtime = Dict{String, Float64}()
    list_filtername = Vector{String}()
    list_filters = Vector{String}()

    filtername = split(lines[1], ";")[1]
    for line in lines
        if occursin("fixiir", line) || filtername != split(line, ";")[1]
            filtername = split(line, ";")[1]
            push!(list_filtername, filtername)
            push!(list_filters, split(line, "_")[1])
            dict_kcmlut[filtername] = parse(Int, split(line, ";")[2])
            dict_kcmtime[filtername] = parse(Float64, split(line, ";")[4][1:(end-2)])
            dict_eps_kcm[filtername] = parse(Float64, split(line, ";")[5])
        end
    end

    str_filtername = ""
    unique!(list_filters)
    str_filtername *= "k="*split(list_filters[1], "x")[end]
    for filtername in list_filters[2:end]
        if !(filtername in ["lp4", "hp0"])
            str_filtername *= ","*"k="*split(filtername, "x")[end]
        else
            str_filtername *= ","*filtername
        end
    end

    figure_header_str = """
    \\begin{figure}$(luts_or_times == 3 ? "\n\\vspace{-3cm}" : "")
    \\centering
    """
    tikz_header_luts_str = """
    \\begin{tikzpicture}
        \\begin{axis}[
                ybar,
                ytick pos=left,
                xtick pos=bottom,
                legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},$(luts_or_times == 3 ? "\n            every axis legend/.code={\\let\\addlegendentry\\relax}," : "")
                ymin=0,
                width=\\textwidth,
                ylabel={\\#LUTs},
                xticklabels={$str_filtername},
                x tick label style={rotate=45, anchor=east},
                xtick=data,
                extra x ticks={3.5, 9.5, 14.5},
                extra x tick labels={lp1,lp2,lp3},
                extra x tick style={yshift=-25pt, major x tick style=transparent, tick label style={rotate=-45}},
                every node near coord/.append style={font=\\footnotesize, rotate=90, anchor=west},
                nodes near coords,
                nodes near coords align={vertical},
                point meta=explicit symbolic
            ]
    """
    tikz_header_times_str = """
    \\begin{tikzpicture}
        \\begin{axis}[
                ybar,
                ytick pos=left,
                xtick pos=bottom,
                legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},
                ymin=0,
                width=\\textwidth,
                ylabel={Time in ns},
                xticklabels={$str_filtername},
                x tick label style={rotate=45, anchor=east},
                xtick=data,
                extra x ticks={3.5, 9.5, 14.5},
                extra x tick labels={lp1,lp2,lp3},
                extra x tick style={yshift=-25pt, major x tick style=transparent, tick label style={rotate=-45}},
                every node near coord/.append style={rotate=90, anchor=west},
                nodes near coords,
                nodes near coords align={vertical},
                point meta=explicit symbolic
            ]
    """

    eps_axis_str = """
    \\end{axis}
    \\begin{axis}[
            %ymode=log,
            axis lines*=right,
            ymin=0.0,
            ytick pos=right,
            xmajorticks=false,
            width=\\textwidth,
            ylabel={error},
            xtick=data,
            every node near coord/.append style={rotate=90, anchor=west},
            nodes near coords,
            nodes near coords align={vertical},
            point meta=explicit symbolic,
        ]
    """

    tikz_footer_luts_str = tikz_footer_times_str = """
        \\end{axis}
    \\end{tikzpicture}
    """

    function figure_footer_str(caption_str)
        return "\\caption{$caption_str}\n\\end{figure}\n\n\\clearpage\n\n"
    end

    plots_str = ""

    initialize_xlabels_str = ""
    initialize_xlabels_eps_str = ""
    initialize_xlabels_str *= "("*string(filters_dict[list_filters[1]])*", 1)"
    initialize_xlabels_eps_str *= "("*string(filters_dict[list_filters[1]])*", 0)"
    for filtername in list_filters[2:end]
        initialize_xlabels_str *= " ("*string(filters_dict[filtername])*", 1)"
        initialize_xlabels_eps_str *= " ("*string(filters_dict[filtername])*", 0)"
    end

    for dw in ["08", "12", "16"]
        dict_kcmlut_str = ""
        dict_kcmtime_str = ""
        dict_eps_kcm_str = ""
        for filtername in list_filtername
            if occursin("_fixiir", filtername) && dw == split(filtername, "_")[2][3:4]
                dict_eps_kcm_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_eps_kcm[filtername])) "
                if acceptable_error < 0 || dict_eps_kcm[filtername] <= acceptable_error
                    dict_kcmlut_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_kcmlut[filtername])) "
                    dict_kcmtime_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_kcmtime[filtername])) "
                end
            end
        end

        plots_str *= figure_header_str
        if luts_or_times in [1,3]
            plots_str *= tikz_header_luts_str
            plots_str *= """
                    \\addplot[white, only marks, draw opacity=0, fill opacity=0] coordinates {$initialize_xlabels_str};
                    \\addplot[legend entry=$(name_fixiir), $(color_fixiir), fill=$(fill_fixiir)] coordinates {$dict_kcmlut_str};
            """

            if with_eps
                plots_str *= eps_axis_str
                plots_str *= """
                        \\addplot[draw=white,ultra thick,draw opacity=0] coordinates {$(initialize_xlabels_eps_str)};
                        \\addplot[draw=$(color_fixiir),ultra thick,mark=*,mark options={solid,draw=$(color_fixiir),fill=$(fill_fixiir)}] coordinates {$dict_eps_kcm_str};
                """
            end
            plots_str *= tikz_footer_luts_str
        end

        if luts_or_times in [2,3]
            plots_str *= tikz_header_times_str
            plots_str *= """
                    \\addplot[white, only marks, draw opacity=0, fill opacity=0] coordinates {$initialize_xlabels_str};
                    \\addplot[legend entry=$(name_fixiir), $(color_fixiir), fill=$(fill_fixiir)] coordinates {$dict_kcmtime_str};
            """

            if with_eps
                plots_str *= eps_axis_str
                plots_str *= """
                        \\addplot[draw=white,ultra thick,draw opacity=0] coordinates {$(initialize_xlabels_eps_str)};
                        \\addplot[draw=$(color_fixiir),ultra thick,mark=*,mark options={solid,draw=$(color_fixiir),fill=$(fill_fixiir)}] coordinates {$dict_eps_kcm_str};
                """
            end
            plots_str *= tikz_footer_times_str
        end

        plots_str *= figure_footer_str("Data wordlength: $(dw == "08" ? "8" : dw)")
    end

    return plots_str
end





function plots_dwcw(
        filename::String;
        acceptable_error::Union{Float64, Int}=-1,
        with_eps::Bool=false,
        wordlengths_range::Tuple{Int, Int}=(4,10),
        color_fixiir::String="blue",
        color_shiftandadd::String="green",
        color_truncatedshiftandadd::String="yellow",
        color_truncatedplain::String="red",
        color_truncatedplain0dsp::String="orange",
        fill_fixiir::String="blue!30!white",
        fill_shiftandadd::String="green!30!white",
        fill_truncatedshiftandadd::String="yellow!30!white",
        fill_truncatedplain::String="red!30!white",
        fill_truncatedplain0dsp::String="orange!30!white",
        name_fixiir::String="fixiir",
        name_shiftandadd::String="shiftandadd",
        name_truncatedplain::String="truncatedplain",
        name_truncatedplain0dsp::String="truncatedplain0dsp",
        name_truncatedshiftandadd::String="truncatedshiftandadd",
        with_fixiir::Bool=false,
        luts_or_times::Int=1, #1 luts, 2 times, 3 both
        printnbadders::Bool=true,
    )
    lines = Vector{String}()

    open(filename, "r") do f
        lines = readlines(f)
    end

    filters_dict = Dict{String, Int}([
        "lp1x0" => 1,
        "lp1x1" => 2,
        "lp1x2" => 3,
        "lp1x3" => 4,
        "lp1x4" => 5,
        "lp1x5" => 6,
        "lp2x0" => 8,
        "lp2x1" => 9,
        "lp2x2" => 10,
        "lp2x3" => 11,
        "lp3x0" => 13,
        "lp3x1" => 14,
        "lp3x2" => 15,
        "lp3x3" => 16,
        "lp4" => 18,
        "hp0" => 20,
    ])

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
    dict_adders_shiftandadd = Dict{String, Int}()
    dict_adders_truncatedshiftandadd = Dict{String, Int}()

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
                dict_adders_shiftandadd[filtername] = parse(Int, split(line, ";")[6])
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
                dict_adders_truncatedshiftandadd[filtername] = parse(Int, split(line, ";")[6])
            else
                error(line)
            end
        end
    end

    str_filtername = ""
    unique!(list_filters)
    str_filtername *= "k="*split(list_filters[1], "x")[end]
    for filtername in list_filters[2:end]
        if !(filtername in ["lp4", "hp0"])
            str_filtername *= ","*"k="*split(filtername, "x")[end]
        else
            str_filtername *= ","*filtername
        end
    end

    figure_header_str = """
    \\begin{figure}$(luts_or_times == 3 ? "\n\\vspace{-3cm}" : "")
    \\centering
    """

    tikz_header_luts_str = """
    \\begin{tikzpicture}
        \\begin{axis}[
                ybar,
                ytick pos=left,
                xtick pos=bottom,
                bar width=$(with_fixiir ? "0.12" : "0.15"),
                legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},$(luts_or_times == 3 ? "\n            every axis legend/.code={\\let\\addlegendentry\\relax}," : "")
                ymin=0,
                width=25cm,
                height=$(luts_or_times == 3 ? "7" : "15")cm,
                ylabel={\\#LUTs},
                xticklabels={$str_filtername},
                xtick=data,
                extra x ticks={3.5, 9.5, 14.5},
                extra x tick labels={lp1,lp2,lp3},
                extra x tick style={yshift=-25pt, major x tick style=transparent},
                point meta=explicit symbolic
            ]
    """

    tikz_header_times_str = """
    \\begin{tikzpicture}
        \\begin{axis}[
                ybar,
                ytick pos=left,
                xtick pos=bottom,
                bar width=$(with_fixiir ? "0.12" : "0.15"),
                legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},
                ymin=0,
                width=25cm,
                height=$(luts_or_times == 3 ? "7" : "15")cm,
                ylabel={Time in ns},
                xticklabels={$str_filtername},
                xtick=data,
                extra x ticks={3.5, 9.5, 14.5},
                extra x tick labels={lp1,lp2,lp3},
                extra x tick style={yshift=-25pt, major x tick style=transparent},
                point meta=explicit symbolic
            ]
    """

    eps_axis_str = """
    \\end{axis}
    \\begin{axis}[
            %ymode=log,
            axis lines*=right,
            ymin=0.0,
            ytick pos=right,
            xmajorticks=false,
            width=25cm,
            height=$(luts_or_times == 3 ? "7" : "15")cm,
            ylabel={error},
            xtick=data,
            every node near coord/.append style={rotate=90, anchor=west},
            nodes near coords,
            nodes near coords align={vertical},
            point meta=explicit symbolic,
        ]
    """

    tikz_footer_luts_str = tikz_footer_times_str = """
        \\end{axis}
    \\end{tikzpicture}
    """

    function figure_footer_str(caption_str)
        return "\\caption{$caption_str}\n\\end{figure}\n\n\\clearpage\n\n"
    end

    plots_str = ""

    initialize_xlabels_str = ""
    initialize_xlabels_eps_str = ""
    initialize_xlabels_str *= "("*string(filters_dict[list_filters[1]])*", 1)"
    initialize_xlabels_eps_str *= "("*string(filters_dict[list_filters[1]])*", 0)"
    for filtername in list_filters[2:end]
        initialize_xlabels_str *= " ("*string(filters_dict[filtername])*", 1)"
        initialize_xlabels_eps_str *= " ("*string(filters_dict[filtername])*", 0)"
    end

    for dw in ["08", "12", "16"]
        dict_kcmlut_str = ""
        dict_kcmtime_str = ""
        dict_eps_kcm_str = ""
        for filtername in list_filtername
            if occursin("_fixiir", filtername) && dw == split(filtername, "_")[2][3:4]
                dict_eps_kcm_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_eps_kcm[filtername])) "
                if acceptable_error < 0 || dict_eps_kcm[filtername] <= acceptable_error
                    dict_kcmlut_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_kcmlut[filtername])) "
                    dict_kcmtime_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_kcmtime[filtername])) "
                end
            end
        end

        for wordlength in wordlengths_range[1]:wordlengths_range[2]
            dict_luts_shiftandadd_str = ""
            dict_lutsdsp_truncatedplain_str = ""
            dict_luts_truncatedplain0dsp_str = ""
            dict_luts_truncatedshiftandadd_str = ""
            dict_times_shiftandadd_str = ""
            dict_times_truncatedplain_str = ""
            dict_times_truncatedplain0dsp_str = ""
            dict_times_truncatedshiftandadd_str = ""
            eps_shiftandadd_str = ""
            eps_truncatedplain_str = ""
            eps_truncatedplain0dsp_str = ""
            eps_truncatedshiftandadd_str = ""
            for filtername in list_filtername
                if occursin("_shiftandadd", filtername) && wordlength == parse(Int, split(filtername, "_")[3][3:4]) && dw == split(filtername, "_")[2][3:4] && get(dict_luts_shiftandadd, filtername, nothing) != nothing
                    eps_shiftandadd_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_eps_shiftandadd[filtername])) "
                    if acceptable_error < 0 || dict_eps_shiftandadd[filtername] <= acceptable_error
                        dict_luts_shiftandadd_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_luts_shiftandadd[filtername]))"
                        dict_times_shiftandadd_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_times_shiftandadd[filtername]))"
                        if printnbadders
                            dict_luts_shiftandadd_str *= "[\$$(dict_adders_shiftandadd[filtername])\$]"
                            dict_times_shiftandadd_str *= "[\$$(dict_adders_shiftandadd[filtername])\$]"
                        end
                        dict_luts_shiftandadd_str *= " "
                        dict_times_shiftandadd_str *= " "
                    end
                elseif occursin("_truncatedplain", filtername) && !occursin("_truncatedplain0dsp", filtername) && wordlength == parse(Int, split(filtername, "_")[3][3:4]) && dw == split(filtername, "_")[2][3:4]
                    eps_truncatedplain_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_eps_truncatedplain[filtername])) "
                    if acceptable_error < 0 || dict_eps_truncatedplain[filtername] <= acceptable_error
                        dict_lutsdsp_truncatedplain_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_luts_truncatedplain[filtername]))"
                        dict_times_truncatedplain_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_times_truncatedplain[filtername]))"
                        if dict_dsps_truncatedplain[filtername] != 0
                            dict_lutsdsp_truncatedplain_str *= "[\$+$(dict_dsps_truncatedplain[filtername])\$ DSPs]"
                        end
                        dict_lutsdsp_truncatedplain_str *= "\n            "
                    end
                elseif occursin("_truncatedplain0dsp", filtername) && wordlength == parse(Int, split(filtername, "_")[3][3:4]) && dw == split(filtername, "_")[2][3:4]
                    eps_truncatedplain_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_eps_truncatedplain0dsp[filtername])) "
                    if acceptable_error < 0 || dict_eps_truncatedplain0dsp[filtername] <= acceptable_error
                        dict_luts_truncatedplain0dsp_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_luts_truncatedplain0dsp[filtername]))"
                        dict_times_truncatedplain0dsp_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_times_truncatedplain0dsp[filtername]))"
                    end
                elseif occursin("_truncatedshiftandadd", filtername) && wordlength == parse(Int, split(filtername, "_")[3][3:4]) && dw == split(filtername, "_")[2][3:4]
                    eps_truncatedshiftandadd_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_eps_truncatedshiftandadd[filtername])) "
                    if acceptable_error < 0 || dict_eps_truncatedshiftandadd[filtername] <= acceptable_error
                        dict_luts_truncatedshiftandadd_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_luts_truncatedshiftandadd[filtername]))"
                        dict_times_truncatedshiftandadd_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_times_truncatedshiftandadd[filtername]))"
                        if printnbadders
                            dict_luts_truncatedshiftandadd_str *= "[\$$(dict_adders_truncatedshiftandadd[filtername])\$]"
                            dict_times_truncatedshiftandadd_str *= "[\$$(dict_adders_truncatedshiftandadd[filtername])\$]"
                        end
                        dict_luts_truncatedshiftandadd_str *= " "
                        dict_times_truncatedshiftandadd_str *= " "
                    end
                end
            end

            plots_str *= figure_header_str
            if luts_or_times in [1, 3]
                plots_str *= tikz_header_luts_str
                plots_str *= "        \\addplot[white, only marks, draw opacity=0, fill opacity=0] coordinates {$initialize_xlabels_str};\n"

                if with_fixiir
                    plots_str *= "        \\addplot[legend entry=$(name_fixiir), $(color_fixiir), fill=$(fill_fixiir)] coordinates {$dict_kcmlut_str};\n"
                end

                plots_str *= """
                        \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_shiftandadd), $(color_shiftandadd), fill=$(fill_shiftandadd)] coordinates {$dict_luts_shiftandadd_str};
                        \\addplot[every node near coord/.append style={font=\\footnotesize, rotate=90, anchor=west, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_truncatedplain), $(color_truncatedplain), fill=$(fill_truncatedplain)] coordinates {$dict_lutsdsp_truncatedplain_str};
                        \\addplot[legend entry=$(name_truncatedplain0dsp), $(color_truncatedplain0dsp), fill=$(fill_truncatedplain0dsp)] coordinates {$dict_luts_truncatedplain0dsp_str};
                        \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_truncatedshiftandadd), $(color_truncatedshiftandadd), fill=$(fill_truncatedshiftandadd)] coordinates {$dict_luts_truncatedshiftandadd_str};
                """

                if with_eps
                    plots_str *= eps_axis_str
                    plots_str *= "        \\addplot[draw=white,ultra thick,draw opacity=0] coordinates {$(initialize_xlabels_eps_str)};\n"
                    if with_fixiir
                        plots_str *= "        \\addplot[draw=$(color_fixiir),ultra thick,mark=*,mark options={solid,draw=$(color_fixiir),fill=$(fill_fixiir)}] coordinates {$dict_eps_kcm_str};\n"
                    end
                    plots_str *= """
                            \\addplot[draw=$(color_shiftandadd),ultra thick,mark=*,mark options={solid,draw=$(color_shiftandadd),fill=$(fill_shiftandadd)}] coordinates {$eps_shiftandadd_str};
                            \\addplot[draw=$(color_truncatedshiftandadd)!70!$(color_truncatedplain),ultra thick,mark=*,mark options={solid,draw=$(color_truncatedshiftandadd)!70!$(color_truncatedplain),fill=$(fill_truncatedshiftandadd)!70!$(fill_truncatedplain)}] coordinates {$eps_truncatedshiftandadd_str};
                    """
                end
                plots_str *= tikz_footer_luts_str
            end

            if luts_or_times in [2, 3]
                plots_str *= tikz_header_times_str
                plots_str *= "        \\addplot[white, only marks, draw opacity=0, fill opacity=0] coordinates {$initialize_xlabels_str};\n"

                if with_fixiir
                    plots_str *= "        \\addplot[legend entry=$(name_fixiir), $(color_fixiir), fill=$(fill_fixiir)] coordinates {$dict_kcmtime_str};\n"
                end

                plots_str *= """
                        \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_shiftandadd), $(color_shiftandadd), fill=$(fill_shiftandadd)] coordinates {$dict_times_shiftandadd_str};
                        \\addplot[every node near coord/.append style={font=\\footnotesize, rotate=90, anchor=west, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_truncatedplain), $(color_truncatedplain), fill=$(fill_truncatedplain)] coordinates {$dict_times_truncatedplain_str};
                        \\addplot[legend entry=$(name_truncatedplain0dsp), $(color_truncatedplain0dsp), fill=$(fill_truncatedplain0dsp)] coordinates {$dict_times_truncatedplain0dsp_str};
                        \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_truncatedshiftandadd), $(color_truncatedshiftandadd), fill=$(fill_truncatedshiftandadd)] coordinates {$dict_times_truncatedshiftandadd_str};
                """

                if with_eps
                    plots_str *= eps_axis_str
                    plots_str *= "        \\addplot[draw=white,ultra thick,draw opacity=0] coordinates {$(initialize_xlabels_eps_str)};\n"
                    if with_fixiir
                        plots_str *= "        \\addplot[draw=$(color_fixiir),ultra thick,mark=*,mark options={solid,draw=$(color_fixiir),fill=$(fill_fixiir)}] coordinates {$dict_eps_kcm_str};\n"
                    end
                    plots_str *= """
                            \\addplot[draw=$(color_shiftandadd),ultra thick,mark=*,mark options={solid,draw=$(color_shiftandadd),fill=$(fill_shiftandadd)}] coordinates {$eps_shiftandadd_str};
                            \\addplot[draw=$(color_truncatedshiftandadd)!70!$(color_truncatedplain),ultra thick,mark=*,mark options={solid,draw=$(color_truncatedshiftandadd)!70!$(color_truncatedplain),fill=$(fill_truncatedshiftandadd)!70!$(fill_truncatedplain)}] coordinates {$eps_truncatedshiftandadd_str};
                    """
                end
                plots_str *= tikz_footer_times_str
            end

            plots_str *= figure_footer_str("Data wordlength: $(dw == "08" ? "8" : dw) -- Coefficients wordlength: $wordlength")
        end
    end

    return plots_str
end







function plots_first_small_eps(
        filename::String;
        acceptable_error::Union{Float64, Int}=-1,
        wordlengths_range::Tuple{Int, Int}=(4,10),
        color_fixiir::String="blue",
        color_shiftandadd::String="green",
        color_truncatedshiftandadd::String="yellow",
        color_truncatedplain::String="red",
        color_truncatedplain0dsp::String="orange",
        fill_fixiir::String="blue!30!white",
        fill_shiftandadd::String="green!30!white",
        fill_truncatedshiftandadd::String="yellow!30!white",
        fill_truncatedplain::String="red!30!white",
        fill_truncatedplain0dsp::String="orange!30!white",
        name_fixiir::String="fixiir",
        name_shiftandadd::String="shiftandadd",
        name_truncatedplain::String="truncatedplain",
        name_truncatedplain0dsp::String="truncatedplain0dsp",
        name_truncatedshiftandadd::String="truncatedshiftandadd",
        printnbadders::Bool=true,
    )
    lines = Vector{String}()

    open(filename, "r") do f
        lines = readlines(f)
    end

    filters_dict = Dict{String, String}([
        "lp1x0" => "1",
        "lp1x1" => "1+1*\\xtickspaces",
        "lp1x2" => "1+2*\\xtickspaces",
        "lp1x3" => "1+3*\\xtickspaces",
        "lp1x4" => "1+4*\\xtickspaces",
        "lp1x5" => "1+5*\\xtickspaces",
        "lp2x0" => "1+5*\\xtickspaces+1*\\xtickblockspaces",
        "lp2x1" => "1+6*\\xtickspaces+1*\\xtickblockspaces",
        "lp2x2" => "1+7*\\xtickspaces+1*\\xtickblockspaces",
        "lp2x3" => "1+8*\\xtickspaces+1*\\xtickblockspaces",
        "lp3x0" => "1+8*\\xtickspaces+2*\\xtickblockspaces",
        "lp3x1" => "1+9*\\xtickspaces+2*\\xtickblockspaces",
        "lp3x2" => "1+10*\\xtickspaces+2*\\xtickblockspaces",
        "lp3x3" => "1+11*\\xtickspaces+2*\\xtickblockspaces",
        "lp4" => "1+11*\\xtickspaces+3*\\xtickblockspaces",
        "hp0" => "1+11*\\xtickspaces+4*\\xtickblockspaces",
    ])

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
    dict_adders_shiftandadd = Dict{String, Int}()
    dict_adders_truncatedshiftandadd = Dict{String, Int}()

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
                dict_adders_shiftandadd[filtername] = parse(Int, split(line, ";")[6])
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
                dict_adders_truncatedshiftandadd[filtername] = parse(Int, split(line, ";")[6])
            else
                error(line)
            end
        end
    end

    str_filtername = ""
    unique!(list_filters)
    if !(list_filters[1] in ["lp4", "hp0"])
        str_filtername *= split(list_filters[1], "x")[1]*"\$_"*split(list_filters[1], "x")[2]*"\$"
    else
        str_filtername *= list_filters[1]
    end
    for filtername in list_filters[2:end]
        if !(filtername in ["lp4", "hp0"])
            str_filtername *= ","*split(filtername, "x")[1]*"\$_"*split(filtername, "x")[2]*"\$"
        else
            str_filtername *= ","*filtername
        end
    end


    function figure_header_str(str_filtername, luts_or_times_str, filters_dict, wordlength_shiftandadd, wordlength_truncated)
        return """
        \\begin{figure}
        \\newcommand\\xtickspaces{1}
        \\newcommand\\xtickblockspaces{1.5}
        \\def\\varwidth{25cm}
	    \\def\\varheight{15cm}
        \\def\\barwidth{0.11}
        \\def\\varyshift{17}
        \\def\\varynegshift{-15}
        \\centering
        \\begin{tikzpicture}
            \\begin{axis}[
                    ybar,
                    ytick pos=left,
                    xtick pos=bottom,
                    bar width=\\barwidth,
                    legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},
                    ymin=0,
                    width=\\varwidth,
                    height=\\varheight,
                    ylabel={$luts_or_times_str},
                    xticklabels={$str_filtername},
                    xtick=data,
                    x tick label style={yshift=\\varynegshift},
                    extra x ticks={$(filters_dict["lp1x0"]), $(filters_dict["lp1x1"]), $(filters_dict["lp1x2"]), $(filters_dict["lp1x3"]), $(filters_dict["lp1x4"]), $(filters_dict["lp1x5"]), $(filters_dict["lp2x0"]), $(filters_dict["lp2x1"]), $(filters_dict["lp2x2"]), $(filters_dict["lp2x3"]), $(filters_dict["lp3x0"]), $(filters_dict["lp3x1"]), $(filters_dict["lp3x2"]), $(filters_dict["lp3x3"]), $(filters_dict["lp4"]), $(filters_dict["hp0"])},
                    extra x tick labels={
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp1x0"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp1x0"])$(wordlength_shiftandadd["lp1x0"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp1x0"] < 10 ? "\\;" : "")$(wordlength_truncated["lp1x0"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp1x0"])$(wordlength_truncated["lp1x0"] == 0 ? "}" : "")$(wordlength_truncated["lp1x0"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp1x1"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp1x1"])$(wordlength_shiftandadd["lp1x1"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp1x1"] < 10 ? "\\;" : "")$(wordlength_truncated["lp1x1"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp1x1"])$(wordlength_truncated["lp1x1"] == 0 ? "}" : "")$(wordlength_truncated["lp1x1"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp1x2"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp1x2"])$(wordlength_shiftandadd["lp1x2"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp1x2"] < 10 ? "\\;" : "")$(wordlength_truncated["lp1x2"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp1x2"])$(wordlength_truncated["lp1x2"] == 0 ? "}" : "")$(wordlength_truncated["lp1x2"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp1x3"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp1x3"])$(wordlength_shiftandadd["lp1x3"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp1x3"] < 10 ? "\\;" : "")$(wordlength_truncated["lp1x3"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp1x3"])$(wordlength_truncated["lp1x3"] == 0 ? "}" : "")$(wordlength_truncated["lp1x3"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp1x4"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp1x4"])$(wordlength_shiftandadd["lp1x4"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp1x4"] < 10 ? "\\;" : "")$(wordlength_truncated["lp1x4"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp1x4"])$(wordlength_truncated["lp1x4"] == 0 ? "}" : "")$(wordlength_truncated["lp1x4"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp1x5"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp1x5"])$(wordlength_shiftandadd["lp1x5"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp1x5"] < 10 ? "\\;" : "")$(wordlength_truncated["lp1x5"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp1x5"])$(wordlength_truncated["lp1x5"] == 0 ? "}" : "")$(wordlength_truncated["lp1x5"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp2x0"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp2x0"])$(wordlength_shiftandadd["lp2x0"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp2x0"] < 10 ? "\\;" : "")$(wordlength_truncated["lp2x0"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp2x0"])$(wordlength_truncated["lp2x0"] == 0 ? "}" : "")$(wordlength_truncated["lp2x0"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp2x1"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp2x1"])$(wordlength_shiftandadd["lp2x1"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp2x1"] < 10 ? "\\;" : "")$(wordlength_truncated["lp2x1"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp2x1"])$(wordlength_truncated["lp2x1"] == 0 ? "}" : "")$(wordlength_truncated["lp2x1"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp2x2"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp2x2"])$(wordlength_shiftandadd["lp2x2"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp2x2"] < 10 ? "\\;" : "")$(wordlength_truncated["lp2x2"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp2x2"])$(wordlength_truncated["lp2x2"] == 0 ? "}" : "")$(wordlength_truncated["lp2x2"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp2x3"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp2x3"])$(wordlength_shiftandadd["lp2x3"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp2x3"] < 10 ? "\\;" : "")$(wordlength_truncated["lp2x3"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp2x3"])$(wordlength_truncated["lp2x3"] == 0 ? "}" : "")$(wordlength_truncated["lp2x3"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp3x0"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp3x0"])$(wordlength_shiftandadd["lp3x0"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp3x0"] < 10 ? "\\;" : "")$(wordlength_truncated["lp3x0"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp3x0"])$(wordlength_truncated["lp3x0"] == 0 ? "}" : "")$(wordlength_truncated["lp3x0"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp3x1"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp3x1"])$(wordlength_shiftandadd["lp3x1"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp3x1"] < 10 ? "\\;" : "")$(wordlength_truncated["lp3x1"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp3x1"])$(wordlength_truncated["lp3x1"] == 0 ? "}" : "")$(wordlength_truncated["lp3x1"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp3x2"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp3x2"])$(wordlength_shiftandadd["lp3x2"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp3x2"] < 10 ? "\\;" : "")$(wordlength_truncated["lp3x2"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp3x2"])$(wordlength_truncated["lp3x2"] == 0 ? "}" : "")$(wordlength_truncated["lp3x2"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp3x3"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp3x3"])$(wordlength_shiftandadd["lp3x3"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp3x3"] < 10 ? "\\;" : "")$(wordlength_truncated["lp3x3"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp3x3"])$(wordlength_truncated["lp3x3"] == 0 ? "}" : "")$(wordlength_truncated["lp3x3"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["lp4"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["lp4"])$(wordlength_shiftandadd["lp4"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["lp4"] < 10 ? "\\;" : "")$(wordlength_truncated["lp4"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["lp4"])$(wordlength_truncated["lp4"] == 0 ? "}" : "")$(wordlength_truncated["lp4"] < 10 ? "\\;" : "")\\;},
                        \\small{\\phantom{K}\\;$(wordlength_shiftandadd["hp0"] == 0 ? "\\phantom{" : "")$(wordlength_shiftandadd["hp0"])$(wordlength_shiftandadd["hp0"] == 0 ? "}" : "")\\;\\;$(wordlength_truncated["hp0"] < 10 ? "\\;" : "")$(wordlength_truncated["hp0"] == 0 ? "\\phantom{" : "")$(wordlength_truncated["hp0"])$(wordlength_truncated["hp0"] == 0 ? "}" : "")$(wordlength_truncated["hp0"] < 10 ? "\\;" : "")\\;}},
                    extra x tick style={yshift=\\varyshift, major x tick style=transparent},
                    point meta=explicit symbolic
                ]
        """
    end

    function figure_footer_str(caption_str)
        return """
            \\end{axis}
        \\end{tikzpicture}
        \\caption{$caption_str}
        \\end{figure}\n\n\\clearpage\n\n
        """
    end

    plots_str = ""

    initialize_xlabels_str = ""
    initialize_xlabels_str *= "("*string(filters_dict[list_filters[1]])*", 1)"
    for filtername in list_filters[2:end]
        initialize_xlabels_str *= " ("*string(filters_dict[filtername])*", 1)"
    end

    for dw in ["08", "12", "16"]
        dict_kcmlut_str = ""
        dict_kcmtime_str = ""
        for filtername in list_filtername
            if occursin("_fixiir", filtername) && dw == split(filtername, "_")[2][3:4]
                if acceptable_error < 0 || dict_eps_kcm[filtername] <= acceptable_error
                    dict_kcmlut_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_kcmlut[filtername])) "
                    dict_kcmtime_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_kcmtime[filtername])) "
                end
            end
        end

        dict_luts_shiftandadd_str = ""
        dict_lutsdsp_truncatedplain_str = ""
        dict_luts_truncatedplain0dsp_str = ""
        dict_luts_truncatedshiftandadd_str = ""
        dict_times_shiftandadd_str = ""
        dict_times_truncatedplain_str = ""
        dict_times_truncatedplain0dsp_str = ""
        dict_times_truncatedshiftandadd_str = ""
        minimumfound_shiftandadd = Dict{String, Int}()
        minimumfound_truncatedplain = Dict{String, Int}()
        minimumfound_truncatedplain0dsp = Dict{String, Int}()
        minimumfound_truncatedshiftandadd = Dict{String, Int}()
        for (key_dict, val_tmp) in filters_dict
            minimumfound_shiftandadd[key_dict] = 0
            minimumfound_truncatedplain[key_dict] = 0
            minimumfound_truncatedplain0dsp[key_dict] = 0
            minimumfound_truncatedshiftandadd[key_dict] = 0
        end
        for wordlength in wordlengths_range[1]:wordlengths_range[2]
            for filtername in list_filtername
                if minimumfound_shiftandadd[split(filtername, "_")[1]] == 0 && occursin("_shiftandadd", filtername) && wordlength == parse(Int, split(filtername, "_")[3][3:4]) && dw == split(filtername, "_")[2][3:4] && get(dict_luts_shiftandadd, filtername, nothing) != nothing
                    if acceptable_error < 0 || dict_eps_shiftandadd[filtername] <= acceptable_error
                        dict_luts_shiftandadd_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_luts_shiftandadd[filtername]))"
                        dict_times_shiftandadd_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_times_shiftandadd[filtername]))"
                        if printnbadders
                            dict_luts_shiftandadd_str *= "[\$$(dict_adders_shiftandadd[filtername])\$]"
                            dict_times_shiftandadd_str *= "[\$$(dict_adders_shiftandadd[filtername])\$]"
                        end
                        dict_luts_shiftandadd_str *= " "
                        dict_times_shiftandadd_str *= " "
                        minimumfound_shiftandadd[split(filtername, "_")[1]] = wordlength
                    end
                elseif minimumfound_truncatedplain[split(filtername, "_")[1]] == 0 && occursin("_truncatedplain", filtername) && !occursin("_truncatedplain0dsp", filtername) && wordlength == parse(Int, split(filtername, "_")[3][3:4]) && dw == split(filtername, "_")[2][3:4]
                    if acceptable_error < 0 || dict_eps_truncatedplain[filtername] <= acceptable_error
                        dict_lutsdsp_truncatedplain_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_luts_truncatedplain[filtername]))"
                        dict_times_truncatedplain_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_times_truncatedplain[filtername])) "
                        if dict_dsps_truncatedplain[filtername] != 0
                            dict_lutsdsp_truncatedplain_str *= "[\\textcolor{black}{\$+$(dict_dsps_truncatedplain[filtername])\$ DSPs}]"
                        end
                        minimumfound_truncatedplain[split(filtername, "_")[1]] = wordlength
                        dict_lutsdsp_truncatedplain_str *= "\n            "
                    end
                elseif minimumfound_truncatedplain0dsp[split(filtername, "_")[1]] == 0 && occursin("_truncatedplain0dsp", filtername) && wordlength == parse(Int, split(filtername, "_")[3][3:4]) && dw == split(filtername, "_")[2][3:4]
                    if acceptable_error < 0 || dict_eps_truncatedplain0dsp[filtername] <= acceptable_error
                        dict_luts_truncatedplain0dsp_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_luts_truncatedplain0dsp[filtername])) "
                        dict_times_truncatedplain0dsp_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_times_truncatedplain0dsp[filtername])) "
                        minimumfound_truncatedplain0dsp[split(filtername, "_")[1]] = wordlength
                    end
                elseif minimumfound_truncatedshiftandadd[split(filtername, "_")[1]] == 0 && occursin("_truncatedshiftandadd", filtername) && wordlength == parse(Int, split(filtername, "_")[3][3:4]) && dw == split(filtername, "_")[2][3:4]
                    if acceptable_error < 0 || dict_eps_truncatedshiftandadd[filtername] <= acceptable_error
                        dict_luts_truncatedshiftandadd_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_luts_truncatedshiftandadd[filtername]))"
                        dict_times_truncatedshiftandadd_str *= "($(string(filters_dict[split(filtername, "_")[1]])), $(dict_times_truncatedshiftandadd[filtername]))"
                        if printnbadders
                            dict_luts_truncatedshiftandadd_str *= "[\$$(dict_adders_truncatedshiftandadd[filtername])\$]"
                            dict_times_truncatedshiftandadd_str *= "[\$$(dict_adders_truncatedshiftandadd[filtername])\$]"
                        end
                        dict_luts_truncatedshiftandadd_str *= " "
                        dict_times_truncatedshiftandadd_str *= " "
                        minimumfound_truncatedshiftandadd[split(filtername, "_")[1]] = wordlength
                    end
                end
            end
        end

        plots_str *= figure_header_str(str_filtername, "\\#LUTs", filters_dict, minimumfound_shiftandadd, minimumfound_truncatedplain)
        plots_str *= "        \\addplot[white, only marks, draw opacity=0, fill opacity=0] coordinates {$initialize_xlabels_str};\n"

        plots_str *= """
                \\addplot[legend entry=$(name_fixiir), $(color_fixiir), fill=$(fill_fixiir)] coordinates {$dict_kcmlut_str};
                \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_shiftandadd), $(color_shiftandadd), fill=$(fill_shiftandadd)] coordinates {$dict_luts_shiftandadd_str};
                \\addplot[every node near coord/.append style={font=\\footnotesize, rotate=90, anchor=west, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_truncatedplain), $(color_truncatedplain), fill=$(fill_truncatedplain)] coordinates {$dict_lutsdsp_truncatedplain_str};
                \\addplot[legend entry=$(name_truncatedplain0dsp), $(color_truncatedplain0dsp), fill=$(fill_truncatedplain0dsp)] coordinates {$dict_luts_truncatedplain0dsp_str};
                \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_truncatedshiftandadd), $(color_truncatedshiftandadd), fill=$(fill_truncatedshiftandadd)] coordinates {$dict_luts_truncatedshiftandadd_str};
        """

        plots_str *= figure_footer_str("Data wordlength: $(dw == "08" ? "8" : dw)")

        plots_str *= figure_header_str(str_filtername, "Times in ns", filters_dict, minimumfound_shiftandadd, minimumfound_truncatedplain)
        plots_str *= "        \\addplot[white, only marks, draw opacity=0, fill opacity=0] coordinates {$initialize_xlabels_str};\n"

        plots_str *= """
                \\addplot[legend entry=$(name_fixiir), $(color_fixiir), fill=$(fill_fixiir)] coordinates {$dict_kcmtime_str};
                \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_shiftandadd), $(color_shiftandadd), fill=$(fill_shiftandadd)] coordinates {$dict_times_shiftandadd_str};
                \\addplot[legend entry=$(name_truncatedplain), $(color_truncatedplain), fill=$(fill_truncatedplain)] coordinates {$dict_times_truncatedplain_str};
                \\addplot[legend entry=$(name_truncatedplain0dsp), $(color_truncatedplain0dsp), fill=$(fill_truncatedplain0dsp)] coordinates {$dict_times_truncatedplain0dsp_str};
                \\addplot[every node near coord/.append style={font=\\footnotesize, anchor=north, color=black}, nodes near coords, nodes near coords align={vertical}, legend entry=$(name_truncatedshiftandadd), $(color_truncatedshiftandadd), fill=$(fill_truncatedshiftandadd)] coordinates {$dict_times_truncatedshiftandadd_str};
        """

        plots_str *= figure_footer_str("Data wordlength: $(dw == "08" ? "8" : dw)")
    end

    return plots_str
end














function plots_transfer_functions(
        filename::String;
        wordlengths_range::Tuple{Int, Int}=(4,10),
        color_fixiir::String="blue",
        color_shiftandadd::String="green",
        color_truncated::String="orange",
        name_fixiir::String="fixiir",
        name_shiftandadd::String="shiftandadd",
        name_truncated::String="truncated",
        filterstype::Vector{String}=["lp1", "lp2", "lp3", "lp4"],
    )
    lines = Vector{String}()

    open(filename, "r") do f
        lines = readlines(f)
    end

    dict_filter_specifications = Dict{String, Tuple{Float64, Float64, Float64}}([
        "lp1x0" => (0.3, 0.7, 0.1),
        "lp1x1" => (0.3, 0.7, 0.09),
        "lp1x2" => (0.3, 0.7, 0.08),
        "lp1x3" => (0.3, 0.7, 0.07),
        "lp1x4" => (0.3, 0.7, 0.06),
        "lp1x5" => (0.3, 0.7, 0.05),
        "lp2x0" => (0.3, 0.7, 0.1),
        "lp2x1" => (0.3, 0.65, 0.1),
        "lp2x2" => (0.3, 0.6, 0.1),
        "lp2x3" => (0.3, 0.55, 0.1),
        "lp3x0" => (0.3, 0.7, 0.1),
        "lp3x1" => (0.35, 0.7, 0.1),
        "lp3x2" => (0.4, 0.7, 0.1),
        "lp3x3" => (0.45, 0.7, 0.1),
        "lp4" => (0.5, 0.9, 0.1),
    ])

    dict_coeffs_shiftandadd = Dict{String, Tuple{Int, Int, Int, Int, Int}}()
    dict_coeffs_truncatedplain = Dict{String, Tuple{Int, Int, Int, Int, Int}}()
    #dict_coeffs_truncatedplain0dsp = Dict{String, Tuple{Int, Int, Int, Int, Int}}()
    #dict_coeffs_truncatedshiftandadd = Dict{String, Tuple{Int, Int, Int, Int, Int}}()
    dict_shifts_shiftandadd = Dict{String, Tuple{Int, Int}}()
    dict_shifts_truncatedplain = Dict{String, Tuple{Int, Int}}()
    #dict_shifts_truncatedplain0dsp = Dict{String, Tuple{Int, Int}}()
    #dict_shifts_truncatedshiftandadd = Dict{String, Tuple{Int, Int}}()
    dict_coeffs_kcm = Dict{String, Tuple{Float64, Float64, Float64, Float64, Float64}}()
    list_filtername = Vector{String}()
    list_filters = Vector{String}()

    for line in lines
        line_vec = split(line)
        filtername = ""
        for val in line_vec
            if startswith(val, "outputFile=")
                filtername = val[(length("outputFile=")+1):end]
            end
        end
        push!(list_filtername, filtername)
        if occursin("fixiir", line)
            push!(list_filters, split(filtername, "_")[1])
            dict_coeffs_kcm[filtername] = Tuple(append!(getindex.(reinterpret.(Float64, reverse.(hex2bytes.(lstrip.(lstrip.(split(split(split(line, "coeffb=\"")[2], "\"")[1], ":"), '0'), 'x')))), 1), getindex.(reinterpret.(Float64, reverse.(hex2bytes.(lstrip.(lstrip.(split(split(split(line, "coeffa=\"")[2], "\"")[1], ":"), '0'), 'x')))), 1)))
        else
            if filtername != split(line, ";")[1]
                push!(list_filters, split(filtername, "_")[1])
            end
            if occursin("_shiftandadd", line)
                dict_coeffs_shiftandadd[filtername] = Tuple(append!(parse.(Int, split(split(split(line, "coeffb=\"")[2], "\"")[1], ":")), parse.(Int, split(split(split(line, "coeffa=\"")[2], "\"")[1], ":"))))
                dict_shifts_shiftandadd[filtername] = (parse(Int, split(split(line, "shiftb=")[2])[1]), parse(Int, split(split(line, "shifta=")[2])[1]))
            elseif occursin("_truncatedplain", line) && !occursin("_truncatedplain0dsp", line)
                dict_coeffs_truncatedplain[filtername] = Tuple(append!(parse.(Int, split(split(split(line, "coeffb=\"")[2], "\"")[1], ":")), parse.(Int, split(split(split(line, "coeffa=\"")[2], "\"")[1], ":"))))
                dict_shifts_truncatedplain[filtername] = (parse(Int, split(split(line, "shiftb=")[2])[1]), parse(Int, split(split(line, "shifta=")[2])[1]))
            elseif occursin("_truncatedplain0dsp", line)
                #nothing
            elseif occursin("_truncatedshiftandadd", line)
                #nothing
            else
                error(line)
            end
        end
    end
    unique!(list_filters)
    # Add special cases
    push!(list_filters, "lp1x5")
    push!(list_filters, "lp2x3")
    push!(list_filters, "lp3x3")

    figure_header_str = """
    \\begin{figure}
    \\centering
    \\begin{tikzpicture}
        \\begin{axis}[
                ytick pos=left,
                xtick pos=bottom,
                trig format=rad,
                legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},
                ymin=0,
                xmin=0,
                xmax=1,
                width=\\textwidth,
                ylabel={\$\\left|H\\left(e^{i\\omega}\\right)\\right|\$},
                ylabel style={rotate=3*pi/2},
                xlabel={\$\\omega\$},
            ]
    """

    function figure_footer_str(caption_str)
        return """
            \\end{axis}
        \\end{tikzpicture}
        \\caption{$caption_str}
        \\end{figure}\n\n\\clearpage\n\n
        """
    end

    plots_str = ""

    for current_filter in list_filters
        if !(split(current_filter, "x")[1] in filterstype)
            continue
        end
        for wordlength in wordlengths_range[1]:wordlengths_range[2]
            plots_str *= figure_header_str
            plots_str *= """
                    \\addplot [domain=0:$(dict_filter_specifications[current_filter][1]),color=black] {1+$(dict_filter_specifications[current_filter][3])};
                    \\addplot [domain=0:$(dict_filter_specifications[current_filter][1]),color=black] {1-$(dict_filter_specifications[current_filter][3])};
                    \\addplot [domain=$(dict_filter_specifications[current_filter][2]):1,color=black] {$(dict_filter_specifications[current_filter][3])};
            """

            for filtername in list_filtername
                if 8 != parse(Int, split(filtername, "_")[2][3:end])
                    continue
                end
                if wordlength != parse(Int, split(filtername, "_")[3][3:end]) && parse(Int, split(filtername, "_")[3][3:end]) != 0
                    continue
                end
                if current_filter != split(filtername, "_")[1]
                    continue
                end

                if occursin("_fixiir", filtername)
                    b0 = dict_coeffs_kcm[filtername][1]
                    b1 = dict_coeffs_kcm[filtername][2]
                    b2 = dict_coeffs_kcm[filtername][3]
                    a0 = 1.0
                    a1 = dict_coeffs_kcm[filtername][4]
                    a2 = dict_coeffs_kcm[filtername][5]

                    plots_str *= """
                            \\addplot [
                                    legend entry=$(name_fixiir),
                                    domain=0:1,
                                    samples=200,
                                    color=$color_fixiir,
                                ]
                                {sqrt(($(b0^2+b1^2+b2^2)+$(2*b0*b2)*cos(2*pi*x)+$(2*b1*b2)*cos(pi*x)+$(2*b0*b1)*cos(pi*x))/($(a0^2+a1^2+a2^2)+$(2*a0*a2)*cos(2*pi*x)+$(2*a1*a2)*cos(pi*x)+$(2*a0*a1)*cos(pi*x)))};
                    """
                elseif occursin("_shiftandadd", filtername)
                    b0 = dict_coeffs_shiftandadd[filtername][1]/(1 << dict_shifts_shiftandadd[filtername][1])
                    b1 = dict_coeffs_shiftandadd[filtername][2]/(1 << dict_shifts_shiftandadd[filtername][1])
                    b2 = dict_coeffs_shiftandadd[filtername][3]/(1 << dict_shifts_shiftandadd[filtername][1])
                    a0 = 1.0
                    a1 = dict_coeffs_shiftandadd[filtername][4]/(1 << dict_shifts_shiftandadd[filtername][2])
                    a2 = dict_coeffs_shiftandadd[filtername][5]/(1 << dict_shifts_shiftandadd[filtername][2])
                    plots_str *= """
                            \\addplot [
                                    legend entry=$(name_shiftandadd),
                                    domain=0:1,
                                    samples=200,
                                    color=$color_shiftandadd,
                                ]
                                {sqrt(($(b0^2+b1^2+b2^2)+$(2*b0*b2)*cos(2*pi*x)+$(2*b1*b2)*cos(pi*x)+$(2*b0*b1)*cos(pi*x))/($(a0^2+a1^2+a2^2)+$(2*a0*a2)*cos(2*pi*x)+$(2*a1*a2)*cos(pi*x)+$(2*a0*a1)*cos(pi*x)))};
                    """
                elseif occursin("_truncatedplain", filtername) && !occursin("_truncatedplain0dsp", filtername)
                    b0 = dict_coeffs_truncatedplain[filtername][1]/(1 << dict_shifts_truncatedplain[filtername][1])
                    b1 = dict_coeffs_truncatedplain[filtername][2]/(1 << dict_shifts_truncatedplain[filtername][1])
                    b2 = dict_coeffs_truncatedplain[filtername][3]/(1 << dict_shifts_truncatedplain[filtername][1])
                    a0 = 1.0
                    a1 = dict_coeffs_truncatedplain[filtername][4]/(1 << dict_shifts_truncatedplain[filtername][2])
                    a2 = dict_coeffs_truncatedplain[filtername][5]/(1 << dict_shifts_truncatedplain[filtername][2])
                    plots_str *= """
                            \\addplot [
                                    legend entry=$(name_truncated),
                                    domain=0:1,
                                    samples=200,
                                    color=$color_truncated,
                                ]
                                {sqrt(($(b0^2+b1^2+b2^2)+$(2*b0*b2)*cos(2*pi*x)+$(2*b1*b2)*cos(pi*x)+$(2*b0*b1)*cos(pi*x))/($(a0^2+a1^2+a2^2)+$(2*a0*a2)*cos(2*pi*x)+$(2*a1*a2)*cos(pi*x)+$(2*a0*a1)*cos(pi*x)))};
                    """
                elseif occursin("_truncatedplain0dsp", filtername)
                    #nothing
                elseif occursin("_truncatedshiftandadd", filtername)
                    #nothing
                else
                    error(filtername)
                end


            end
            plots_str *= figure_footer_str("$current_filter -- Coefficient wordlength: $wordlength")
        end
    end

    return plots_str
end
























function plots_eps(
        filename::String;
        acceptable_errors::Union{Vector{Union{Float64, Int}}, Vector{Float64}, Vector{Int}}=[-1],
        color_fixiir::String="blue",
        color_shiftandadd::String="green",
        color_truncatedshiftandadd::String="yellow",
        color_truncatedplain::String="red",
        color_truncatedplain0dsp::String="orange",
        fill_fixiir::String="blue!30!white",
        fill_shiftandadd::String="green!30!white",
        fill_truncatedshiftandadd::String="yellow!30!white",
        fill_truncatedplain::String="red!30!white",
        fill_truncatedplain0dsp::String="orange!30!white",
        name_fixiir::String="fixiir",
        name_shiftandadd::String="shiftandadd",
        name_truncatedplain::String="truncatedplain",
        name_truncatedplain0dsp::String="truncatedplain0dsp",
        name_truncatedshiftandadd::String="truncatedshiftandadd",
        lp_shapes::Dict{String, String}=Dict{String, String}(["lp1"=>"*", "lp2"=>"triangle*", "lp3"=>"square*", "lp4"=>"diamond*"]),
        benchmarks_print::Vector{String}=["lp1", "lp2", "lp3", "lp4"],
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

    filtername = split(lines[1], ";")[1]
    for line in lines
        if occursin("fixiir", line)
            filtername = split(line, ";")[1]
            push!(list_filtername, filtername)
            push!(list_filters, split(line, "_")[1])
            dict_kcmlut[filtername] = parse(Int, split(line, ";")[2])
            dict_kcmtime[filtername] = parse(Float64, split(line, ";")[4][1:(end-2)])
            dict_eps_kcm[filtername] = parse(Float64, split(line, ";")[5])
        else
            if filtername != split(line, ";")[1]
                filtername = split(line, ";")[1]
                push!(list_filtername, filtername)
                push!(list_filters, split(line, "_")[1])
            end
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

    str_filtername = ""
    unique!(list_filters)
    if !(list_filters[1] in ["lp4", "hp0"])
        str_filtername *= split(list_filters[1], "x")[1]*"\$_"*split(list_filters[1], "x")[2]*"\$"
    else
        str_filtername *= list_filters[1]
    end
    for filtername in list_filters[2:end]
        if !(filtername in ["lp4", "hp0"])
            str_filtername *= ","*split(filtername, "x")[1]*"\$_"*split(filtername, "x")[2]*"\$"
        else
            str_filtername *= ","*filtername
        end
    end

    figure_header_luts_str = """
    \\begin{figure}
    \\centering
    \\begin{tikzpicture}
        \\begin{axis}[
                ytick pos=left,
                xtick pos=bottom,
                legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},
                width=\\textwidth,
                ylabel={error},
                xlabel={\\#LUTs},
            ]
    """

    figure_header_times_str = """
    \\begin{figure}
    \\centering
    \\begin{tikzpicture}
        \\begin{axis}[
                ytick pos=left,
                xtick pos=bottom,
                legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},
                width=\\textwidth,
                ylabel={error},
                xlabel={Time in ns},
            ]
    """

    function figure_footer_str(caption_str)
        return """
            \\end{axis}
        \\end{tikzpicture}
        \\caption{$caption_str}
        \\end{figure}\n\n\\clearpage\n\n
        """
    end

    plots_str = ""

    for acceptable_error in acceptable_errors
        for dwordlength in 8:4:16
            dict_kcmlut_str = Dict{String, String}()
            for filtertype in benchmarks_print
                dict_kcmlut_str[filtertype] = ""
            end
            dict_luts_shiftandadd_str = Dict{String, String}()
            for filtertype in benchmarks_print
                dict_luts_shiftandadd_str[filtertype] = ""
            end
            dict_luts_truncatedplain_str = Dict{String, String}()
            for filtertype in benchmarks_print
                dict_luts_truncatedplain_str[filtertype] = ""
            end
            dict_luts_truncatedplain0dsp_str = Dict{String, String}()
            for filtertype in benchmarks_print
                dict_luts_truncatedplain0dsp_str[filtertype] = ""
            end
            dict_luts_truncatedshiftandadd_str = Dict{String, String}()
            for filtertype in benchmarks_print
                dict_luts_truncatedshiftandadd_str[filtertype] = ""
            end

            dict_kcmtime_str = Dict{String, String}()
            for filtertype in benchmarks_print
                dict_kcmtime_str[filtertype] = ""
            end
            dict_times_shiftandadd_str = Dict{String, String}()
            for filtertype in benchmarks_print
                dict_times_shiftandadd_str[filtertype] = ""
            end
            dict_times_truncatedplain_str = Dict{String, String}()
            for filtertype in benchmarks_print
                dict_times_truncatedplain_str[filtertype] = ""
            end
            dict_times_truncatedplain0dsp_str = Dict{String, String}()
            for filtertype in benchmarks_print
                dict_times_truncatedplain0dsp_str[filtertype] = ""
            end
            dict_times_truncatedshiftandadd_str = Dict{String, String}()
            for filtertype in benchmarks_print
                dict_times_truncatedshiftandadd_str[filtertype] = ""
            end

            for filtername in list_filtername
                if !occursin(string(dwordlength), split(filtername, "_dw")[2][1:2])
                    continue
                end
                if occursin("_fixiir", filtername) && (acceptable_error < 0 || dict_eps_kcm[filtername] <= acceptable_error)
                    for filtertype in benchmarks_print
                        if occursin(filtertype, filtername)
                            dict_kcmlut_str[filtertype] *= "($(dict_kcmlut[filtername]), $(dict_eps_kcm[filtername])) "
                            dict_kcmtime_str[filtertype] *= "($(dict_kcmtime[filtername]), $(dict_eps_kcm[filtername])) "
                        end
                    end
                elseif occursin("_shiftandadd", filtername) && (acceptable_error < 0 || dict_eps_shiftandadd[filtername] <= acceptable_error)
                    for filtertype in benchmarks_print
                        if occursin(filtertype, filtername)
                            dict_luts_shiftandadd_str[filtertype] *= "($(dict_luts_shiftandadd[filtername]), $(dict_eps_shiftandadd[filtername])) "
                            dict_times_shiftandadd_str[filtertype] *= "($(dict_times_shiftandadd[filtername]), $(dict_eps_shiftandadd[filtername])) "
                        end
                    end
                elseif occursin("_truncatedplain", filtername) && !occursin("_truncatedplain0dsp", filtername) && (acceptable_error < 0 || dict_eps_truncatedplain[filtername] <= acceptable_error)
                    for filtertype in benchmarks_print
                        if occursin(filtertype, filtername)
                            dict_luts_truncatedplain_str[filtertype] *= "($(dict_luts_truncatedplain[filtername]), $(dict_eps_truncatedplain[filtername])) "
                            dict_times_truncatedplain_str[filtertype] *= "($(dict_times_truncatedplain[filtername]), $(dict_eps_truncatedplain[filtername])) "
                        end
                    end
                elseif occursin("_truncatedplain0dsp", filtername) && (acceptable_error < 0 || dict_eps_truncatedplain0dsp[filtername] <= acceptable_error)
                    for filtertype in benchmarks_print
                        if occursin(filtertype, filtername)
                            dict_luts_truncatedplain0dsp_str[filtertype] *= "($(dict_luts_truncatedplain0dsp[filtername]), $(dict_eps_truncatedplain0dsp[filtername])) "
                            dict_times_truncatedplain0dsp_str[filtertype] *= "($(dict_times_truncatedplain0dsp[filtername]), $(dict_eps_truncatedplain0dsp[filtername])) "
                        end
                    end
                elseif occursin("_truncatedshiftandadd", filtername) && (acceptable_error < 0 || dict_eps_truncatedshiftandadd[filtername] <= acceptable_error)
                    for filtertype in benchmarks_print
                        if occursin(filtertype, filtername)
                            dict_luts_truncatedshiftandadd_str[filtertype] *= "($(dict_luts_truncatedshiftandadd[filtername]), $(dict_eps_truncatedshiftandadd[filtername])) "
                            dict_times_truncatedshiftandadd_str[filtertype] *= "($(dict_times_truncatedshiftandadd[filtername]), $(dict_eps_truncatedshiftandadd[filtername])) "
                        end
                    end
                end
            end

            plots_str *= figure_header_luts_str

            for filtertype in benchmarks_print[1:(end-1)]
                plots_str *= "        \\addplot[mark=$(lp_shapes[filtertype]), legend entry=\\phantom{0}, only marks, color=$(color_fixiir), fill=$(fill_fixiir)] coordinates {$(dict_kcmlut_str[filtertype])};\n"
            end
            plots_str *= "        \\addplot[mark=$(lp_shapes[benchmarks_print[end]]), legend entry=$(name_fixiir), only marks, color=$(color_fixiir), fill=$(fill_fixiir)] coordinates {$(dict_kcmlut_str[benchmarks_print[end]])};\n"
            for filtertype in benchmarks_print[1:(end-1)]
                plots_str *= "        \\addplot[mark=$(lp_shapes[filtertype]), legend entry=\\phantom{0}, only marks, color=$(color_shiftandadd), fill=$(fill_shiftandadd)] coordinates {$(dict_luts_shiftandadd_str[filtertype])};\n"
            end
            plots_str *= "        \\addplot[mark=$(lp_shapes[benchmarks_print[end]]), legend entry=$(name_shiftandadd), only marks, color=$(color_shiftandadd), fill=$(fill_shiftandadd)] coordinates {$(dict_luts_shiftandadd_str[benchmarks_print[end]])};\n"
            for filtertype in benchmarks_print[1:(end-1)]
                plots_str *= "        \\addplot[mark=$(lp_shapes[filtertype]), legend entry=\\phantom{0}, only marks, color=$(color_truncatedplain), fill=$(fill_truncatedplain)] coordinates {$(dict_luts_truncatedplain_str[filtertype])};\n"
            end
            plots_str *= "        \\addplot[mark=$(lp_shapes[benchmarks_print[end]]), legend entry=$(name_truncatedplain), only marks, color=$(color_truncatedplain), fill=$(fill_truncatedplain)] coordinates {$(dict_luts_truncatedplain_str[benchmarks_print[end]])};\n"
            for filtertype in benchmarks_print[1:(end-1)]
                plots_str *= "        \\addplot[mark=$(lp_shapes[filtertype]), legend entry=\\phantom{0}, only marks, color=$(color_truncatedplain0dsp), fill=$(fill_truncatedplain0dsp)] coordinates {$(dict_luts_truncatedplain0dsp_str[filtertype])};\n"
            end
            plots_str *= "        \\addplot[mark=$(lp_shapes[benchmarks_print[end]]), legend entry=$(name_truncatedplain0dsp), only marks, color=$(color_truncatedplain0dsp), fill=$(fill_truncatedplain0dsp)] coordinates {$(dict_luts_truncatedplain0dsp_str[benchmarks_print[end]])};\n"
            for filtertype in benchmarks_print[1:(end-1)]
                plots_str *= "        \\addplot[mark=$(lp_shapes[filtertype]), legend entry=\\phantom{0}, only marks, color=$(color_truncatedshiftandadd), fill=$(fill_truncatedshiftandadd)] coordinates {$(dict_luts_truncatedshiftandadd_str[filtertype])};\n"
            end
            plots_str *= "        \\addplot[mark=$(lp_shapes[benchmarks_print[end]]), legend entry=$(name_truncatedshiftandadd), only marks, color=$(color_truncatedshiftandadd), fill=$(fill_truncatedshiftandadd)] coordinates {$(dict_luts_truncatedshiftandadd_str[benchmarks_print[end]])};\n"

            plots_str *= figure_footer_str("Error/LUTs -- $(acceptable_error) -- $benchmarks_print")


            plots_str *= figure_header_times_str

            for filtertype in benchmarks_print[1:(end-1)]
                plots_str *= "        \\addplot[mark=$(lp_shapes[filtertype]), legend entry=\\phantom{0}, only marks, color=$(color_fixiir), fill=$(fill_fixiir)] coordinates {$(dict_kcmtime_str[filtertype])};\n"
            end
            plots_str *= "        \\addplot[mark=$(lp_shapes[benchmarks_print[end]]), legend entry=$(name_fixiir), only marks, color=$(color_fixiir), fill=$(fill_fixiir)] coordinates {$(dict_kcmtime_str[benchmarks_print[end]])};\n"
            for filtertype in benchmarks_print[1:(end-1)]
                plots_str *= "        \\addplot[mark=$(lp_shapes[filtertype]), legend entry=\\phantom{0}, only marks, color=$(color_shiftandadd), fill=$(fill_shiftandadd)] coordinates {$(dict_times_shiftandadd_str[filtertype])};\n"
            end
            plots_str *= "        \\addplot[mark=$(lp_shapes[benchmarks_print[end]]), legend entry=$(name_shiftandadd), only marks, color=$(color_shiftandadd), fill=$(fill_shiftandadd)] coordinates {$(dict_times_shiftandadd_str[benchmarks_print[end]])};\n"
            for filtertype in benchmarks_print[1:(end-1)]
                plots_str *= "        \\addplot[mark=$(lp_shapes[filtertype]), legend entry=\\phantom{0}, only marks, color=$(color_truncatedplain), fill=$(fill_truncatedplain)] coordinates {$(dict_times_truncatedplain_str[filtertype])};\n"
            end
            plots_str *= "        \\addplot[mark=$(lp_shapes[benchmarks_print[end]]), legend entry=$(name_truncatedplain), only marks, color=$(color_truncatedplain), fill=$(fill_truncatedplain)] coordinates {$(dict_times_truncatedplain_str[benchmarks_print[end]])};\n"
            for filtertype in benchmarks_print[1:(end-1)]
                plots_str *= "        \\addplot[mark=$(lp_shapes[filtertype]), legend entry=\\phantom{0}, only marks, color=$(color_truncatedplain0dsp), fill=$(fill_truncatedplain0dsp)] coordinates {$(dict_times_truncatedplain0dsp_str[filtertype])};\n"
            end
            plots_str *= "        \\addplot[mark=$(lp_shapes[benchmarks_print[end]]), legend entry=$(name_truncatedplain0dsp), only marks, color=$(color_truncatedplain0dsp), fill=$(fill_truncatedplain0dsp)] coordinates {$(dict_times_truncatedplain0dsp_str[benchmarks_print[end]])};\n"
            for filtertype in benchmarks_print[1:(end-1)]
                plots_str *= "        \\addplot[mark=$(lp_shapes[filtertype]), legend entry=\\phantom{0}, only marks, color=$(color_truncatedshiftandadd), fill=$(fill_truncatedshiftandadd)] coordinates {$(dict_times_truncatedshiftandadd_str[filtertype])};\n"
            end
            plots_str *= "        \\addplot[mark=$(lp_shapes[benchmarks_print[end]]), legend entry=$(name_truncatedshiftandadd), only marks, color=$(color_truncatedshiftandadd), fill=$(fill_truncatedshiftandadd)] coordinates {$(dict_times_truncatedshiftandadd_str[benchmarks_print[end]])};\n"

            plots_str *= figure_footer_str("Error/Time -- Data wordlength: $dwordlength -- $(acceptable_error) -- $benchmarks_print")
        end
    end

    return plots_str
end



















function main_plots(filename::String)
    latex_header_str = """
    \\documentclass{article}
    \\usepackage{a4wide}
    \\usepackage[utf8]{inputenc}
    \\usepackage{graphicx} % Allows including images
    \\usepackage{amssymb}
    \\usepackage[usenames, dvipsnames]{xcolor}
    \\usepackage{pgfplots}
    \\pgfplotsset{compat=1.12}
    \\usepackage{amsmath}
    \\usepackage{tikz}

    \\pgfplotsset{
        legend entry/.initial=,
        every axis plot post/.code={%
            \\pgfkeysgetvalue{/pgfplots/legend entry}\\tempValue
            \\ifx\\tempValue\\empty
            \\pgfkeysalso{/pgfplots/forget plot}%
            \\else
            \\expandafter\\addlegendentry\\expandafter{\\tempValue}%
            \\fi
        },
    }

    \\definecolor{twostepkcm}{RGB}{0, 0, 255}
    \\definecolor{twostepkcmfill}{RGB}{100, 100, 255}
    \\definecolor{twostepmcm}{RGB}{230, 230, 0}
    \\definecolor{twostepmcmfill}{RGB}{230, 230, 100}
    \\definecolor{twostepgeneric}{RGB}{255, 0, 0}
    \\definecolor{twostepgenericfill}{RGB}{255, 100, 100}
    \\definecolor{twostepgenericnodsp}{RGB}{255, 127, 0}
    \\definecolor{twostepgenericnodspfill}{RGB}{255, 227, 100}
    \\definecolor{ours}{RGB}{0, 255, 0}
    \\definecolor{oursfill}{RGB}{100, 255, 100}
    \\definecolor{errorplotlines}{RGB}{255, 0, 0}
    \\definecolor{errorplotcenterdots}{RGB}{0, 0, 0}
    \\definecolor{errorplotdots}{RGB}{255, 0, 0}

    \\begin{document}\n\n
    """

    latex_footer_str = "\n\n\\end{document}"

    latex_header_ls_str = """
    \\documentclass{article}
    \\usepackage{a4wide}
    \\usepackage[utf8]{inputenc}
    \\usepackage{graphicx} % Allows including images
    \\usepackage{amssymb}
    \\usepackage[usenames, dvipsnames]{xcolor}
    \\usepackage{pgfplots}
    \\pgfplotsset{compat=1.12}
    \\usepackage{amsmath}
    \\usepackage{tikz}
    \\usepackage{lscape}

    \\pgfplotsset{
        legend entry/.initial=,
        every axis plot post/.code={%
            \\pgfkeysgetvalue{/pgfplots/legend entry}\\tempValue
            \\ifx\\tempValue\\empty
            \\pgfkeysalso{/pgfplots/forget plot}%
            \\else
            \\expandafter\\addlegendentry\\expandafter{\\tempValue}%
            \\fi
        },
    }

    \\definecolor{twostepkcm}{RGB}{0, 0, 255}
    \\definecolor{twostepkcmfill}{RGB}{100, 100, 255}
    \\definecolor{twostepmcm}{RGB}{230, 230, 0}
    \\definecolor{twostepmcmfill}{RGB}{230, 230, 100}
    \\definecolor{twostepgeneric}{RGB}{255, 0, 0}
    \\definecolor{twostepgenericfill}{RGB}{255, 100, 100}
    \\definecolor{twostepgenericnodsp}{RGB}{255, 127, 0}
    \\definecolor{twostepgenericnodspfill}{RGB}{255, 227, 100}
    \\definecolor{ours}{RGB}{0, 255, 0}
    \\definecolor{oursfill}{RGB}{100, 255, 100}
    \\definecolor{errorplotlines}{RGB}{255, 0, 0}
    \\definecolor{errorplotcenterdots}{RGB}{0, 0, 0}
    \\definecolor{errorplotdots}{RGB}{255, 0, 0}

    \\begin{document}
    \\begin{landscape}\n\n
    """

    latex_footer_ls_str = "\n\n\\end{landscape}\n\\end{document}"



    plots_folder = (@__DIR__)*"/../plots/"


    options = [
        ("_luts", "", -1, false, 1),
        ("_luts", "_eps", 1e-7, false, 1),
        ("_luts", "_with_eps", -1, true, 1),
        ("_times", "", -1, false, 2),
        ("_times", "_eps", 1e-7, false, 2),
        ("_times", "_with_eps", -1, true, 2),
        ("_lutstimes", "", -1, false, 3),
        ("_lutstimes", "_eps", 1e-7, false, 3),
        ("_lutstimes", "_with_eps", -1, true, 3),
    ]

    filenamewrite = plots_folder*"all_plots_filters"
    for (name_type, name_type_error, acceptable_error, with_eps, luts_or_times) in options
        plots_str = plots_for_each_filter(filename, acceptable_error=acceptable_error, with_eps=with_eps, luts_or_times=luts_or_times,
            color_fixiir="twostepkcm", fill_fixiir="twostepkcmfill",
            color_shiftandadd="ours", fill_shiftandadd="oursfill",
            color_truncatedshiftandadd="twostepmcm", fill_truncatedshiftandadd="twostepmcmfill",
            color_truncatedplain="twostepgeneric", fill_truncatedplain="twostepgenericfill",
            color_truncatedplain0dsp="twostepgenericnodsp", fill_truncatedplain0dsp="twostepgenericnodspfill",
            name_fixiir="2-step KCM",
            name_shiftandadd="Ours",
            name_truncatedplain="2-step Generic",
            name_truncatedplain0dsp="2-step Generic NoDSP",
            name_truncatedshiftandadd="2-step MCM"
        )
        open(filenamewrite*name_type*name_type_error*".tex", "w") do f
            write(f, latex_header_str)
            write(f, plots_str)
            write(f, latex_footer_str)
        end
    end


    filenamewrite = plots_folder*"all_plots_dwcw"
    for (name_type, name_type_error, acceptable_error, with_eps, luts_or_times) in options
        plots_str = plots_dwcw(filename, acceptable_error=acceptable_error, with_eps=with_eps, luts_or_times=luts_or_times,
            color_fixiir="twostepkcm", fill_fixiir="twostepkcmfill",
            color_shiftandadd="ours", fill_shiftandadd="oursfill",
            color_truncatedshiftandadd="twostepmcm", fill_truncatedshiftandadd="twostepmcmfill",
            color_truncatedplain="twostepgeneric", fill_truncatedplain="twostepgenericfill",
            color_truncatedplain0dsp="twostepgenericnodsp", fill_truncatedplain0dsp="twostepgenericnodspfill",
            name_fixiir="2-step KCM",
            name_shiftandadd="Ours",
            name_truncatedplain="2-step Generic",
            name_truncatedplain0dsp="2-step Generic NoDSP",
            name_truncatedshiftandadd="2-step MCM"
        )
        open(filenamewrite*name_type*name_type_error*".tex", "w") do f
            write(f, latex_header_ls_str)
            write(f, plots_str)
            write(f, latex_footer_ls_str)
        end
    end


    filenamewrite = plots_folder*"all_plots_all"
    for (name_type, name_type_error, acceptable_error, with_eps, luts_or_times) in options
        plots_str = plots_dwcw(filename, acceptable_error=acceptable_error, with_eps=with_eps, luts_or_times=luts_or_times, with_fixiir=true,
            color_fixiir="twostepkcm", fill_fixiir="twostepkcmfill",
            color_shiftandadd="ours", fill_shiftandadd="oursfill",
            color_truncatedshiftandadd="twostepmcm", fill_truncatedshiftandadd="twostepmcmfill",
            color_truncatedplain="twostepgeneric", fill_truncatedplain="twostepgenericfill",
            color_truncatedplain0dsp="twostepgenericnodsp", fill_truncatedplain0dsp="twostepgenericnodspfill",
            name_fixiir="2-step KCM",
            name_shiftandadd="Ours",
            name_truncatedplain="2-step Generic",
            name_truncatedplain0dsp="2-step Generic NoDSP",
            name_truncatedshiftandadd="2-step MCM"
        )
        open(filenamewrite*name_type*name_type_error*".tex", "w") do f
            write(f, latex_header_ls_str)
            write(f, plots_str)
            write(f, latex_footer_ls_str)
        end
    end


    filenamewrite = plots_folder*"all_plots_fixiir"
    for (name_type, name_type_error, acceptable_error, with_eps, luts_or_times) in options
        plots_str = plots_fixiir(filename, acceptable_error=acceptable_error, with_eps=with_eps, luts_or_times=luts_or_times,
            color_fixiir="twostepkcm", fill_fixiir="twostepkcmfill",
            color_shiftandadd="ours", fill_shiftandadd="oursfill",
            color_truncatedshiftandadd="twostepmcm", fill_truncatedshiftandadd="twostepmcmfill",
            color_truncatedplain="twostepgeneric", fill_truncatedplain="twostepgenericfill",
            color_truncatedplain0dsp="twostepgenericnodsp", fill_truncatedplain0dsp="twostepgenericnodspfill",
            name_fixiir="2-step KCM",
            name_shiftandadd="Ours",
            name_truncatedplain="2-step Generic",
            name_truncatedplain0dsp="2-step Generic NoDSP",
            name_truncatedshiftandadd="2-step MCM"
        )
        open(filenamewrite*name_type*name_type_error*".tex", "w") do f
            write(f, latex_header_str)
            write(f, plots_str)
            write(f, latex_footer_str)
        end
    end


    plots_str = plots_first_small_eps(filename, acceptable_error=1e-7,wordlengths_range=(4,32), printnbadders=false,
        color_fixiir="twostepkcm", fill_fixiir="twostepkcmfill",
        color_shiftandadd="ours", fill_shiftandadd="oursfill",
        color_truncatedshiftandadd="twostepmcm", fill_truncatedshiftandadd="twostepmcmfill",
        color_truncatedplain="twostepgeneric", fill_truncatedplain="twostepgenericfill",
        color_truncatedplain0dsp="twostepgenericnodsp", fill_truncatedplain0dsp="twostepgenericnodspfill",
        name_fixiir="2-step KCM",
        name_shiftandadd="Ours",
        name_truncatedplain="2-step Generic",
        name_truncatedplain0dsp="2-step Generic NoDSP",
        name_truncatedshiftandadd="2-step MCM"
    )
    filenamewrite = plots_folder*"all_plots_recap"
    open(filenamewrite*".tex", "w") do f
        write(f, latex_header_ls_str)
        write(f, plots_str)
        write(f, latex_footer_ls_str)
    end

    filenamewrite = plots_folder*"plots_eps"
    open(filenamewrite*".tex", "w") do f
        write(f, latex_header_str)
        plots_str = plots_eps(filename, acceptable_errors=[-1, 1e-7, 1e-2],
            color_fixiir="twostepkcm", fill_fixiir="twostepkcmfill",
            color_shiftandadd="ours", fill_shiftandadd="oursfill",
            color_truncatedshiftandadd="twostepmcm", fill_truncatedshiftandadd="twostepmcmfill",
            color_truncatedplain="twostepgeneric", fill_truncatedplain="twostepgenericfill",
            color_truncatedplain0dsp="twostepgenericnodsp", fill_truncatedplain0dsp="twostepgenericnodspfill",
            name_fixiir="2-step KCM",
            name_shiftandadd="Ours",
            name_truncatedplain="2-step Generic",
            name_truncatedplain0dsp="2-step Generic NoDSP",
            name_truncatedshiftandadd="2-step MCM"
        )
        write(f, plots_str)
        lp_shapes = Dict{String, String}(["lp1"=>"*", "lp2"=>"triangle*", "lp3"=>"square*", "lp4"=>"diamond*"])
        for filtertype in ["lp1", "lp2", "lp3", "lp4", "lp1x0", "lp1x1", "lp1x2", "lp1x3", "lp1x4", "lp1x5", "lp2x0", "lp2x1", "lp2x2", "lp2x3", "lp3x0", "lp3x1", "lp3x2", "lp3x3"]
            lp_shapes[filtertype] = lp_shapes[string(split(filtertype, "x")[1])]
            plots_str = plots_eps(filename, acceptable_errors=[-1, 1e-7, 1e-2], lp_shapes=lp_shapes, benchmarks_print=[filtertype],
                color_fixiir="twostepkcm", fill_fixiir="twostepkcmfill",
                color_shiftandadd="ours", fill_shiftandadd="oursfill",
                color_truncatedshiftandadd="twostepmcm", fill_truncatedshiftandadd="twostepmcmfill",
                color_truncatedplain="twostepgeneric", fill_truncatedplain="twostepgenericfill",
                color_truncatedplain0dsp="twostepgenericnodsp", fill_truncatedplain0dsp="twostepgenericnodspfill",
                name_fixiir="2-step KCM",
                name_shiftandadd="Ours",
                name_truncatedplain="2-step Generic",
                name_truncatedplain0dsp="2-step Generic NoDSP",
                name_truncatedshiftandadd="2-step MCM"
            )
            write(f, plots_str)
        end
        write(f, latex_footer_str)
    end



    filenamewrite = plots_folder*"plots_transfer_functions"
    for filtertype in ["lp1", "lp2", "lp3", "lp4"]
        plots_str = plots_transfer_functions((@__DIR__)*"/../../filtergeneration/flopococalls/FloPoCoCalls.txt", filterstype=[filtertype])
        open(filenamewrite*"_"*filtertype*".tex", "w") do f
            write(f, latex_header_str)
            write(f, plots_str)
            write(f, latex_footer_str)
        end
    end

    return nothing
end


function only_transfer_function_plots()
    latex_header_str = """
    \\documentclass{article}
    \\usepackage{a4wide}
    \\usepackage[utf8]{inputenc}
    \\usepackage{graphicx} % Allows including images
    \\usepackage{amssymb}
    \\usepackage[usenames, dvipsnames]{xcolor}
    \\usepackage{pgfplots}
    \\pgfplotsset{compat=1.12}
    \\usepackage{amsmath}
    \\usepackage{tikz}

    \\pgfplotsset{
        legend entry/.initial=,
        every axis plot post/.code={%
            \\pgfkeysgetvalue{/pgfplots/legend entry}\\tempValue
            \\ifx\\tempValue\\empty
            \\pgfkeysalso{/pgfplots/forget plot}%
            \\else
            \\expandafter\\addlegendentry\\expandafter{\\tempValue}%
            \\fi
        },
    }

    \\begin{document}\n\n
    """

    latex_footer_str = "\n\n\\end{document}"

    plots_folder = (@__DIR__)*"/../plots/"

    filenamewrite = plots_folder*"plots_transfer_functions"
    for filtertype in ["lp1", "lp2", "lp3", "lp4"]
        plots_str = plots_transfer_functions((@__DIR__)*"/../../filtergeneration/flopococalls/FloPoCoCalls.txt", filterstype=[filtertype])
        open(filenamewrite*"_"*filtertype*".tex", "w") do f
            write(f, latex_header_str)
            write(f, plots_str)
            write(f, latex_footer_str)
        end
    end

    return nothing
end




function solutions_to_plots()
    latex_header_str = """
    \\documentclass{article}
    \\usepackage{a4wide}
    \\usepackage[utf8]{inputenc}
    \\usepackage{graphicx} % Allows including images
    \\usepackage{amssymb}
    \\usepackage[usenames, dvipsnames]{xcolor}
    \\usepackage{pgfplots}
    \\pgfplotsset{compat=1.12}
    \\usepackage{amsmath}
    \\usepackage{tikz}

    \\pgfplotsset{
        legend entry/.initial=,
        every axis plot post/.code={%
            \\pgfkeysgetvalue{/pgfplots/legend entry}\\tempValue
            \\ifx\\tempValue\\empty
            \\pgfkeysalso{/pgfplots/forget plot}%
            \\else
            \\expandafter\\addlegendentry\\expandafter{\\tempValue}%
            \\fi
        },
    }

    \\begin{document}\n\n
    """

    latex_footer_str = "\n\n\\end{document}"

    tikz_header_str = """
    \\begin{figure}
    \\begin{tikzpicture}
        \\begin{axis}[
                ybar,
                bar width=0.15,
                legend style={at={(0.5,-0.15)}, anchor=north, legend columns=-1},
                ymin=0,
                ytick pos=left,
                xtick pos=bottom,
                xtick=data,
                width=\\textwidth,
                ylabel={\\#solutions},
            ]
    """
    function tikz_footer_str(caption::String)
        return "    \\end{axis}\n\\end{tikzpicture}\n\\caption{$caption}\n\\end{figure}"
    end

    filenames = [
        "lp1x4",
    ]

    open((@__DIR__)*"/../plots/mcm_solutions_plots.tex", "w") do fwrite
        write(fwrite, latex_header_str)
        nbAddersDict = Dict{Set{Int}, Int}()
        for filename in filenames
            for wordlength in 4:10
                whole_filename = (@__DIR__)*"/../../filtergeneration/solutions/$(filename)_$(wordlength).txt"
                if isfile(whole_filename)
                    println("-- $(filename)_$(wordlength) --")
                    nbAddersSum = Dict{Int, Int}()
                    open(whole_filename) do fread
                        lines = readlines(fread)
                        for line in lines[2:end]
                            a1, a2, b0, b1, b2, _, _ = odd.(abs.(parse.(Int, split(line, ";"))))
                            modela = Model(CPLEX.Optimizer)
                            modelb = Model(CPLEX.Optimizer)
                            set_silent(modela)
                            set_silent(modelb)
                            nbAdders = get!(nbAddersDict, Set{Int}([a1, a2]),
                                length(mcm(modela, [a1, a2], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true))) +
                                       get!(nbAddersDict, Set{Int}([b0, b1, b2]),
                                length(mcm(modelb, [b0, b1, b2], wordlength=wordlength, use_nlp=false, ilp=1, use_big_m=true))) +
                                sum(coef != 0 ? 1 : 0 for coef in [a1, a2, b0, b1, b2]) - 1
                            if !haskey(nbAddersSum, nbAdders)
                                nbAddersSum[nbAdders] = 0
                            end
                            nbAddersSum[nbAdders] += 1
                        end
                    end
                    if !isempty(nbAddersSum)
                        write(fwrite, "\n\n")
                        write(fwrite, tikz_header_str)
                        str_nb_adders = ""
                        for (key, val) in nbAddersSum
                            str_nb_adders *= "($(key), $(val)) "
                        end
                        write(fwrite, "        \\addplot[blue, fill=blue!30!white] coordinates {$(str_nb_adders)};\n")
                        write(fwrite, tikz_footer_str("$(filename)\\_$(wordlength)"))
                        write(fwrite, "\n\n")
                    end
                end
            end
        end
        write(fwrite, latex_footer_str)
    end

    return nothing
end
