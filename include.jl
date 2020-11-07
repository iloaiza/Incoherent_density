using LinearAlgebra
using FastGaussQuadrature
using SpecialFunctions
using HDF5
using FFTW

include("config.jl")
include("saving.jl")
include("general_functions.jl")
include("sanity.jl")
include("LVC_pot.jl")
include("rhodopsin.jl")
include("heaviside.jl")
include("krylov.jl")
include("lindblad.jl")
include("tipt.jl")
include("shift.jl")
include("PLOT.jl")

using Plots
using LaTeXStrings
using Plots.PlotMeasures
FONT = font(44)
SIZE = [1720,880]
pyplot()


function fdebug()
	include("fdebug.jl")
end