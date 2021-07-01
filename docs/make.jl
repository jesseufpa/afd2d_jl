using Yeppp
using Plots 
using BenchmarkTools
using DelimitedFiles
using Printf

include("../src/SU_mod.jl")
include("../src/KaiserWindow_mod.jl")
include("../src/LanczosInterpolation_mod.jl")
include("../src/InvSmooth2d_mod.jl")
include("../src/lsq_cgsolver.jl")

using .LSQ_CGSolver
using .SU_mod
using .KaiserWindow_mod
using .LanczosInterpolation_mod
using .InvSmooth2d_mod

using Documenter, AFD2d_mod

makedocs(sitename="AFD2d_docs",modules=[Documenter,AFD2d_mod])
