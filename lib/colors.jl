# import Base.println
# import Base.print
#
# export println
# export print
# export stylize
# export multistylize
#
#
# function multistylize(x::AbstractString, colors::Union{Vector, Tuple})
#   res = x
#   for color in colors
#     res = stylize(res, color)
#   end
#   return res
# end
#
# function stylize(x::AbstractString, color::Symbol)
#   colors = Dict(
#     :reset => (0, 0),
#
#     :bold => (1, 22),
#     :dim => (2, 22),
#     :italic => (3, 23),
#     :underline => (4, 24),
#     :inverse => (7, 27),
#     :hidden => (8, 28),
#     :strikethrough => (9, 29),
#
#     :black => (30, 39),
#     :red => (31, 39),
#     :green => (32, 39),
#     :yellow => (33, 39),
#     :blue => (34, 39),
#     :magenta => (35, 39),
#     :cyan => (36, 39),
#     :white => (37, 39),
#     :gray => (90, 39),
#     :grey => (90, 39),
#
#     :bgBlack => (40, 49),
#     :bgRed => (41, 49),
#     :bgGreen => (42, 49),
#     :bgYellow => (43, 49),
#     :bgBlue => (44, 49),
#     :bgMagenta => (45, 49),
#     :bgCyan => (46, 49),
#     :bgWhite => (47, 49),
#
#     # =>l(gacy styles for colors pre v1.0)0
#     :blackBG => (40, 49),
#     :redBG => (41, 49),
#     :greenBG => (42, 49),
#     :yellowBG => (43, 49),
#     :blueBG => (44, 49),
#     :magentaBG => (45, 49),
#     :cyanBG => (46, 49),
#     :whiteBG => (47, 49)
#   )
#   return string("\u001b[", colors[color][1] ,"m", x, "\u001b[", colors[color][2] ,"m")
# end
#
# println(x::AbstractString; color::Symbol) = println(stylize(x, color))
# print(x::AbstractString; color::Symbol) = print(stylize(x, color))
#
# println(x::AbstractString; colors::Union{Vector, Tuple}) = println(multistylize(x, colors))
# print(x::AbstractString; colors::Union{Vector, Tuple}) = print(multistylize(x, colors))
