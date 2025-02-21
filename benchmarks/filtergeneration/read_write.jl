const flopococalls_base_name = (@__DIR__)*"/flopococalls/"*"FloPoCoCalls"
const rpagcalls_base_name = (@__DIR__)*"/rpagcalls/"
const error_epsilon_base_name = (@__DIR__)*"/epsilon/table_epsilon.csv"

function initialize_files()
    close(open(flopococalls_base_name*"_plain.txt", "w"))
    close(open(flopococalls_base_name*"_shiftandadd.txt", "w"))
    close(open(flopococalls_base_name*"_truncated_plain.txt", "w"))
    close(open(flopococalls_base_name*"_truncated_shiftandadd.txt", "w"))
    close(open(flopococalls_base_name*"_fixiir.txt", "w"))
    close(open(flopococalls_base_name*".txt", "w"))
    close(open(flopococalls_base_name*"_with_python.txt", "w"))
    close(open(error_epsilon_base_name, "w"))
    close(open(rpagcalls_base_name*"truncated.csv", "w"))
    close(open(rpagcalls_base_name*"truncated_calls.csv", "w"))
    return nothing
end

################################################################################
function write_flopoco_solution_plain(solution::Iir2AdderGraphs, instance::String, wordlength::Int, wcpg_val::Union{Float64, Int}, wcpg_eps_val::Union{Float64, Int})
    open(flopococalls_base_name*"_plain.txt", "a") do io
        for lsbin in -16:4:-8
            write(io, "flopoco FixIIRShiftAdd msbIn=-1 lsbIn=$(lsbin) lsbOut=$(lsbin) ")
            write(io, "guardBits=-1 ")
            write(io, "coeffa=\"$(solution.coefficients_a[1]):$(solution.coefficients_a[2])\" ")
            write(io, "coeffb=\"$(solution.coefficients_b[1]):$(solution.coefficients_b[2]):$(solution.coefficients_b[3])\" ")
            write(io, "shifta=$(solution.shifts[1]) shiftb=$(solution.shifts[2]) ")
            write(io, "method=\"plain\" ")
            write(io, "outputFile=$(instance)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_plain.vhd")
            if wcpg_val != 0 || wcpg_eps_val != 0
                write(io, " H=$(wcpg_val) Heps=$(wcpg_eps_val)")
            end
            write(io, "\n")
        end
    end
    return nothing
end

function write_flopoco_solution_shiftandadd(solution::Iir2AdderGraphs, instance::String, wordlength::Int, wcpg_val::Union{Float64, Int}, wcpg_eps_val::Union{Float64, Int})
    open(flopococalls_base_name*"_shiftandadd.txt", "a") do io
        for lsbin in -16:4:-8
            write(io, "flopoco FixIIRShiftAdd msbIn=-1 lsbIn=$(lsbin) lsbOut=$(lsbin) ")
            write(io, "guardBits=-1 ")
            write(io, "coeffa=\"$(solution.coefficients_a[1]):$(solution.coefficients_a[2])\" ")
            write(io, "coeffb=\"$(solution.coefficients_b[1]):$(solution.coefficients_b[2]):$(solution.coefficients_b[3])\" ")
            write(io, "shifta=$(solution.shifts[1]) shiftb=$(solution.shifts[2]) ")
            write(io, "graphb=\"$(write_addergraph(solution.addergraph_b))\" ")
            write(io, "grapha=\"$(write_addergraph(solution.addergraph_a))\" ")
            write(io, "method=\"multiplierless\" ")
            write(io, "outputFile=$(instance)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_shiftandadd.vhd")
            if wcpg_val != 0 || wcpg_eps_val != 0
                write(io, " H=$(wcpg_val) Heps=$(wcpg_eps_val)")
            end
            write(io, "\n")
        end
    end
    return nothing
end

function write_flopoco_solution(solution::Iir2AdderGraphs, instance::String, wordlength::Int, wcpg_val::Union{Float64, Int}, wcpg_eps_val::Union{Float64, Int})
    write_flopoco_solution_plain(solution, instance, wordlength, wcpg_val, wcpg_eps_val)
    write_flopoco_solution_shiftandadd(solution, instance, wordlength, wcpg_val, wcpg_eps_val)
    return nothing
end




################################################################################
function write_flopoco_solution_plain(solution::Vector{Any}, instance::String, wordlength::Int, shifta::Int, shiftb::Int)
    open(flopococalls_base_name*"_truncated_plain.txt", "a") do io
        for lsbin in -16:4:-8
            write(io, "flopoco FixIIRShiftAdd msbIn=-1 lsbIn=$(lsbin) lsbOut=$(lsbin) ")
            write(io, "guardBits=-1 ")
            write(io, "coeffa=\"$(solution[1][1]):$(solution[1][2])\" ")
            write(io, "coeffb=\"$(solution[2][1]):$(solution[2][2]):$(solution[2][3])\" ")
            write(io, "shifta=$(shifta) shiftb=$(shiftb) ")
            write(io, "method=\"plain\" ")
            write(io, "outputFile=$(instance)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain.vhd\n")
        end
    end
    return nothing
end

function write_flopoco_solution_shiftandadd(solution::Vector{Any}, instance::String, wordlength::Int, shifta::Int, shiftb::Int)
    open(flopococalls_base_name*"_truncated_shiftandadd.txt", "a") do io
        for lsbin in -16:4:-8
            write(io, "flopoco FixIIRShiftAdd msbIn=-1 lsbIn=$(lsbin) lsbOut=$(lsbin) ")
            write(io, "guardBits=-1 ")
            write(io, "coeffa=\"$(solution[1][1]):$(solution[1][2])\" ")
            write(io, "coeffb=\"$(solution[2][1]):$(solution[2][2]):$(solution[2][3])\" ")
            write(io, "shifta=$(shifta) shiftb=$(shiftb) ")
            write(io, "graphb=\"$(write_addergraph(solution[4]))\" ")
            write(io, "grapha=\"$(write_addergraph(solution[3]))\" ")
            write(io, "method=\"multiplierless\" ")
            write(io, "outputFile=$(instance)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd.vhd\n")
        end
    end
    return nothing
end

function write_flopoco_solution(solution::Vector{Any}, instance::String, wordlength::Int, shifta::Int, shiftb::Int)
    write_flopoco_solution_plain(solution, instance, wordlength, shifta, shiftb)
    write_flopoco_solution_shiftandadd(solution, instance, wordlength, shifta, shiftb)
    return nothing
end


function write_flopoco_solution(b0::Float64, b1::Float64, b2::Float64, a1::Float64, a2::Float64, instancename::String)
    filename = flopococalls_base_name*"_fixiir.txt"
    open(filename, "a") do f
        for lsbin in -16:4:-8
            write(f, "flopoco FixIIR lsbIn=$lsbin lsbOut=$lsbin ")
            write(f, "guardBits=-1 ")
            write(f, "coeffb=\"")
            write(f, "$b0:")
            write(f, "$b1:")
            write(f, "$b2\" ")
            write(f, "coeffa=\"")
            write(f, "$a1:")
            write(f, "$a2\" ")
            write(f, "outputFile=$(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw00_fixiir.vhd\n")
        end
    end
end


function write_flopoco_solution(b0_hex::String, b1_hex::String, b2_hex::String, a1_hex::String, a2_hex::String, instancename::String)
    filename = flopococalls_base_name*"_fixiir.txt"
    open(filename, "a") do f
        for lsbin in -16:4:-8
            write(f, "flopoco FixIIR lsbIn=$lsbin lsbOut=$lsbin ")
            write(f, "guardBits=-1 ")
            write(f, "coeffb=\"")
            write(f, "0x$b0_hex:")
            write(f, "0x$b1_hex:")
            write(f, "0x$b2_hex\" ")
            write(f, "coeffa=\"")
            write(f, "0x$a1_hex:")
            write(f, "0x$a2_hex\" ")
            write(f, "outputFile=$(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw00_fixiir.vhd\n")
        end
    end
end



################################################################################
function write_table_epsilon(instance::String, NA::Int, varepsilon::Float64, wordlength::Int)
    open(error_epsilon_base_name, "a") do f
        for lsbin in 16:-4:8
            write(f, "$(replace(instance, "dw00" => "dw$(lsbin < 10 ? "0" : "")$lsbin"));$NA;$varepsilon;$wordlength;$lsbin")
            write(f, "\n")
        end
    end
    return nothing
end





################################################################################
function write_rpag_calls(instancename::String, b0, b1, b2, a1, a2, wordlength, varepsilon, shifta, shiftb)
    filename = rpagcalls_base_name*"truncated_calls.csv"
    open(filename, "a") do f
        write(f, "\"$instancename\";")
        write(f, "\"coeffa\";")
        write(f, "\"$a1, $a2\"")
        write(f, "\n")
        write(f, "\"$instancename\";")
        write(f, "\"coeffb\";")
        write(f, "\"$b0, $b1, $b2\";")
        write(f, "\n")
    end

    filename = rpagcalls_base_name*"truncated.csv"
    open(filename, "a") do f
        write(f, "\"$instancename\";")
        write(f, "\"coeffa\";")
        write(f, "\"$a1, $a2\"")
        write(f, "\n")
        write(f, "\"$instancename\";")
        write(f, "\"coeffb\";")
        write(f, "\"$b0, $b1, $b2\";")
        write(f, "\"$wordlength\";")
        write(f, "\"$varepsilon\";")
        write(f, "\"$shifta\";")
        write(f, "\"$shiftb\";")
        write(f, "\n")
    end
end


function read_rpag_outputs()
    lines = Vector{String}()
    filename = rpagcalls_base_name*"truncated_with_adder_graph.csv"
    open(filename, "r") do f
        lines = readlines(f)
    end
    lines_wl_eps = Vector{String}()
    filename = rpagcalls_base_name*"truncated.csv"
    open(filename, "r") do f
        lines_wl_eps = readlines(f)
    end

    open(flopococalls_base_name*".txt", "w") do f
        open(flopococalls_base_name*"_with_python.txt", "w") do fwp
            for i in 1:2:length(lines)
                wordlength = parse(Int, strip(split(lines_wl_eps[i+1], ";")[4], '\"'))
                shifta = parse(Int, strip(split(lines_wl_eps[i+1], ";")[6], '\"'))
                shiftb = parse(Int, strip(split(lines_wl_eps[i+1], ";")[7], '\"'))
                varepsilon = parse(Float64, strip(split(lines_wl_eps[i+1], ";")[5], '\"'))
                instancename = strip(split(lines[i], ";")[1], '\"')
                solution_a = strip(split(lines[i], ";")[5], '\"')
                solution_b = strip(split(lines[i+1], ";")[5], '\"')
                nbAa = parse(Int, strip(split(lines[i], ";")[4], '\"'))
                nbAb = parse(Int, strip(split(lines[i+1], ";")[4], '\"'))
                coeffa = parse.(Int, strip.(split(split(lines[i], ";")[3], ","), '\"'))
                coeffb = parse.(Int, strip.(split(split(lines[i+1], ";")[3], ","), '\"'))
                nbcoefsA = length(filter(x -> x != 0, coeffa))
                nbcoefsB = length(filter(x -> x != 0, coeffb))

                for lsbin in -16:4:-8
                    write(f, "flopoco FixIIRShiftAdd msbIn=-1 lsbIn=$(lsbin) lsbOut=$(lsbin) ")
                    write(f, "guardBits=-1 ")
                    write(f, "coeffa=\"$(coeffa[1]):$(coeffa[2])\" ")
                    write(f, "coeffb=\"$(coeffb[1]):$(coeffb[2]):$(coeffb[3])\" ")
                    write(f, "shifta=$(shifta) shiftb=$(shiftb) ")
                    write(f, "graphb=\"$(solution_b)\" ")
                    write(f, "grapha=\"$(solution_a)\" ")
                    write(f, "method=\"multiplierless\" ")
                    write(f, "outputFile=$(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd.vhd\n")

                    write(fwp, "flopoco FixIIRShiftAdd msbIn=-1 lsbIn=$(lsbin) lsbOut=$(lsbin) ")
                    write(fwp, "guardBits=-1 ")
                    write(fwp, "coeffa=\"$(coeffa[1]):$(coeffa[2])\" ")
                    write(fwp, "coeffb=\"$(coeffb[1]):$(coeffb[2]):$(coeffb[3])\" ")
                    write(fwp, "shifta=$(shifta) shiftb=$(shiftb) ")
                    write(fwp, "graphb=\"$(solution_b)\" ")
                    write(fwp, "grapha=\"$(solution_a)\" ")
                    write(fwp, "method=\"multiplierless\" ")
                    write(fwp, "outputFile=$(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd.vhd\n")
                    write(fwp, "python ../tools/vivado-runsyn.py --implement --vhdl $(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd.vhd\n")

                    write(f, "flopoco FixIIRShiftAdd msbIn=-1 lsbIn=$(lsbin) lsbOut=$(lsbin) ")
                    write(f, "guardBits=-1 ")
                    write(f, "coeffa=\"$(coeffa[1]):$(coeffa[2])\" ")
                    write(f, "coeffb=\"$(coeffb[1]):$(coeffb[2]):$(coeffb[3])\" ")
                    write(f, "shifta=$(shifta) shiftb=$(shiftb) ")
                    write(f, "method=\"plain\" ")
                    write(f, "outputFile=$(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain.vhd\n")

                    write(fwp, "flopoco FixIIRShiftAdd msbIn=-1 lsbIn=$(lsbin) lsbOut=$(lsbin) ")
                    write(fwp, "guardBits=-1 ")
                    write(fwp, "coeffa=\"$(coeffa[1]):$(coeffa[2])\" ")
                    write(fwp, "coeffb=\"$(coeffb[1]):$(coeffb[2]):$(coeffb[3])\" ")
                    write(fwp, "shifta=$(shifta) shiftb=$(shiftb) ")
                    write(fwp, "method=\"plain\" ")
                    write(fwp, "outputFile=$(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain.vhd\n")
                    write(fwp, "python ../tools/vivado-runsyn.py --implement --vhdl $(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain.vhd\n")

                    write(f, "flopoco FixIIRShiftAdd msbIn=-1 lsbIn=$(lsbin) lsbOut=$(lsbin) ")
                    write(f, "guardBits=-1 ")
                    write(f, "coeffa=\"$(coeffa[1]):$(coeffa[2])\" ")
                    write(f, "coeffb=\"$(coeffb[1]):$(coeffb[2]):$(coeffb[3])\" ")
                    write(f, "shifta=$(shifta) shiftb=$(shiftb) ")
                    write(f, "method=\"plain\" ")
                    write(f, "outputFile=$(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain0dsp.vhd\n")

                    write(fwp, "flopoco FixIIRShiftAdd msbIn=-1 lsbIn=$(lsbin) lsbOut=$(lsbin) ")
                    write(fwp, "guardBits=-1 ")
                    write(fwp, "coeffa=\"$(coeffa[1]):$(coeffa[2])\" ")
                    write(fwp, "coeffb=\"$(coeffb[1]):$(coeffb[2]):$(coeffb[3])\" ")
                    write(fwp, "shifta=$(shifta) shiftb=$(shiftb) ")
                    write(fwp, "method=\"plain\" ")
                    write(fwp, "outputFile=$(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain0dsp.vhd\n")
                    write(fwp, "python ../tools/vivado-runsyn.py --implement --vhdl $(instancename)_dw$(abs(lsbin) < 10 ? "0" : "")$(abs(lsbin))_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain0dsp.vhd --maxdsp 0\n")
                end

                write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedshiftandadd", nbAa+nbAb+nbcoefsA+nbcoefsB-1, varepsilon, wordlength)
                write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain", 0, varepsilon, wordlength)
                write_table_epsilon("$(instancename)_dw00_cw$(wordlength < 10 ? "0" : "")$(wordlength)_truncatedplain0dsp", 0, varepsilon, wordlength)
            end
        end
    end
end







################################################################################
function read_and_rewrite()
    filename_fixiir = flopococalls_base_name*"_fixiir.txt"
    filename_onestep = flopococalls_base_name*"_shiftandadd.txt"
    filename_truncated_plain = flopococalls_base_name*"_truncated_plain.txt"
    filename_threesteps = flopococalls_base_name*"_truncated_shiftandadd.txt"

    lines_fixiir = Vector{String}()
    lines_onestep = Vector{String}()
    lines_truncated_plain = Vector{String}()
    lines_threesteps = Vector{String}()

    open(filename_fixiir, "r") do f
        lines_fixiir = readlines(f)
    end
    open(filename_onestep, "r") do f
        lines_onestep = readlines(f)
    end
    open(filename_truncated_plain, "r") do f
        lines_truncated_plain = readlines(f)
    end
    open(filename_threesteps, "r") do f
        lines_threesteps = readlines(f)
    end

    open(flopococalls_base_name*".txt", "a") do f
        open(flopococalls_base_name*"_with_python.txt", "a") do fpython
            for line_fixiir in lines_fixiir
                write(f, line_fixiir*"\n")
                line_fixiir_vec = split(line_fixiir)
                outputname = ""
                for val in line_fixiir_vec
                    if startswith(val, "outputFile=")
                        outputname = val[(length("outputFile=")+1):end]
                    end
                end
                write(fpython, line_fixiir*"\n")
                write(fpython, "python ../tools/vivado-runsyn.py --implement --vhdl $outputname"*"\n")

                for wordlength in 4:32
                    for line_truncated_plain in lines_truncated_plain
                        line_vec = split(line_truncated_plain)
                        current_outputname = ""
                        for val in line_vec
                            if startswith(val, "outputFile=")
                                current_outputname = val[(length("outputFile=")+1):end]
                            end
                        end

                        if occursin("cw$(wordlength < 10 ? "0" : "")$(wordlength)", current_outputname) && current_outputname[1:(first(findfirst("_cw", current_outputname))-1)] == outputname[1:(first(findfirst("_cw00_fixiir", outputname))-1)]
                            write(f, line_truncated_plain*"\n")
                            write(fpython, line_truncated_plain*"\n")
                            write(fpython, "python ../tools/vivado-runsyn.py --implement --vhdl $current_outputname"*"\n")
                            write(f, "$(replace(line_truncated_plain, ".vhd" => "0dsp.vhd"))\n")
                            write(fpython, "$(replace(line_truncated_plain, ".vhd" => "0dsp.vhd"))\n")
                            write(fpython, "python ../tools/vivado-runsyn.py --implement --vhdl $(current_outputname[1:(end-4)])0dsp.vhd --maxdsp 0\n")

                        end
                    end

                    for line_threesteps in lines_threesteps
                        line_vec = split(line_threesteps)
                        current_outputname = ""
                        for val in line_vec
                            if startswith(val, "outputFile=")
                                current_outputname = val[(length("outputFile=")+1):end]
                            end
                        end
                        if occursin("cw$(wordlength < 10 ? "0" : "")$(wordlength)", current_outputname) && current_outputname[1:(first(findfirst("_cw", current_outputname))-1)] == outputname[1:(first(findfirst("_cw00_fixiir", outputname))-1)]
                            write(f, line_threesteps*"\n")
                            write(fpython, line_threesteps*"\n")
                            write(fpython, "python ../tools/vivado-runsyn.py --implement --vhdl $current_outputname"*"\n")
                        end
                    end

                    for line_onestep in lines_onestep
                        line_vec = split(line_onestep)
                        current_outputname = ""
                        for val in line_vec
                            if startswith(val, "outputFile=")
                                current_outputname = val[(length("outputFile=")+1):end]
                            end
                        end
                        if occursin("cw$(wordlength < 10 ? "0" : "")$(wordlength)", current_outputname) && current_outputname[1:(first(findfirst("_cw", current_outputname))-1)] == outputname[1:(first(findfirst("_cw00_fixiir", outputname))-1)]
                            write(f, line_onestep*"\n")
                            write(fpython, line_onestep*"\n")
                            write(fpython, "python ../tools/vivado-runsyn.py --implement --vhdl $current_outputname"*"\n")
                        end
                    end
                end
            end
        end
    end
    return nothing
end


function write_special_cases()
    filename_onestep = flopococalls_base_name*"_shiftandadd.txt"
    lines_onestep = Vector{String}()
    open(filename_onestep, "r") do f
        lines_onestep = readlines(f)
    end

    open(flopococalls_base_name*".txt", "a") do f
        open(flopococalls_base_name*"_with_python.txt", "a") do fpython
            for outputname in ["lp1x5_dw16_cw00_fixiir.vhd" , "lp1x5_dw12_cw00_fixiir.vhd", "lp1x5_dw08_cw00_fixiir.vhd", "lp2x3_dw16_cw00_fixiir.vhd" , "lp2x3_dw12_cw00_fixiir.vhd", "lp2x3_dw08_cw00_fixiir.vhd", "lp3x3_dw16_cw00_fixiir.vhd" , "lp3x3_dw12_cw00_fixiir.vhd", "lp3x3_dw08_cw00_fixiir.vhd"]
                for wordlength in 4:32
                    for line_onestep in lines_onestep
                        line_vec = split(line_onestep)
                        current_outputname = ""
                        for val in line_vec
                            if startswith(val, "outputFile=")
                                current_outputname = val[(length("outputFile=")+1):end]
                            end
                        end
                        if occursin("cw$(wordlength < 10 ? "0" : "")$(wordlength)", current_outputname) && current_outputname[1:(first(findfirst("_cw", current_outputname))-1)] == outputname[1:(first(findfirst("_cw00_fixiir", outputname))-1)]
                            write(f, line_onestep*"\n")
                            write(fpython, line_onestep*"\n")
                            write(fpython, "python ../tools/vivado-runsyn.py --implement --vhdl $current_outputname"*"\n")
                        end
                    end
                end
            end
        end
    end
    return nothing
end
