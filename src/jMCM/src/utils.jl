# function odd(number::Int)
#     if number == 0
#         return 0
#     end
#     while mod(number, 2) == 0
#         number = div(number, 2)
#     end
#     return number
# end

function get_min_wordlength(number::Int)
    return round(Int, max(log2(odd(abs(number))), 1), RoundUp)
end
