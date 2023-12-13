using DrWatson;
@quickactivate "LatticeModels";
include(srcdir("LatticeModels.jl"))
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
include(srcdir("../parameters.jl"));
using BenchmarkTools, CUDA

## ----------------- Test NRXY Model on GPU ----------------- ##
## ----------------- Test NRXY Model on GPU ----------------- ##
include(srcdir("../parameters.jl"));
tmax = Tf(1000)
t = Tf(0)
lattice = SquareLattice(L)
thetas = CUDA.rand(Tf, L, L, R) * Tf(2pi)
thetas_new = CUDA.zeros(Tf, L, L, R) 

times = collect(0:tmax/20:tmax)

CUDA.@time evolve!(thetas, thetas_new, lattice, params, tmax)

## ----------------- Test HydroXY Model on GPU ----------------- ##
## ----------------- Test HydroXY Model on GPU ----------------- ##
include(srcdir("../parameters.jl"));
tmax = Tf(1000)




