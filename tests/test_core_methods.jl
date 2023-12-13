using DrWatson;
@quickactivate "LatticeModels";
include(srcdir("LatticeModels.jl"))
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();
include(srcdir("../parameters.jl"));
using BenchmarkTools

## Tests get neighbours
# Physical Parameters
include(srcdir("../parameters.jl"));
params["symmetry"] = "polar"
params["algo"] = "Langevin"

## Benchmark update
model = NRXY(params)
model = XY(params)
lattice = SquareLattice(L, periodic=true)
# model = DiscreteNRXY(params)
lattice = TriangularLattice(L, periodic=true)
thetas = init_thetas(model, lattice, params)
update!(thetas, model, lattice)
@btime update!($thetas, $model, $lattice)
@time evolve!(thetas, model, lattice, 2000)
ndef = number_defects(thetas, model, lattice)
plot_thetas(thetas, model, lattice,title="t = $(round(model.t,digits=1)), n = $ndef")
prinz(10000*100) 


##
# 53 s for tmax = 1000 
#= Runtimes
WN = Wrapped Normal ; VM = VonMises ; AF = antiferro
Avec WN, AF = false update takes 11 ms
Avec VM, AF = false update takes 23 ms
Avec WN, AF = true update takes 11 ms
Avec VM, AF = true update takes 25 ms
=#

## Tests update!()
model = XY(params)
@btime update!(thetas, model, lattice)

model = ForcedXY(params)
update!(thetas, model, lattice)
@btime update!(thetas, model, lattice)

model = VisionConeXY(params)
update!(thetas, model, lattice)
@btime update!(thetas, model, lattice)

model = DiscreteNRXY(params)
update!(thetas, model, lattice)
@btime update!(thetas, model, lattice)

model = NRXY(params)
update!(thetas, model, lattice)
@btime update!(thetas, model, lattice)

model = SPP(params)
thetas = init_thetas(model, lattice, init="2pair", q=1, r0=60, type=["source", "divergent"])
update!(thetas, model, lattice)
@btime update!(thetas, model, lattice)
plot_theta(thetas, model, lattice)

# Tests update!(....,tmax)
tmax = 1
model = VisionConeXY(params)
model.t
update!(thetas, model, lattice, tmax)
model.t

## Visual verification
include(srcdir("../parameters.jl"));
model = NRXY(params)
lattice = TriangularLattice(L)
# lattice = SquareLattice(L)
thetas = init_thetas(model, lattice, params_init=params_init)
plot_thetas(thetas, model, lattice)
tmax = 1000
z = @elapsed evolve!(thetas, model, lattice, tmax=tmax)
plot_thetas(thetas, model, lattice)
## ok tout semble correspondre à mes attentes
prinz(z)
plot_thetas(thetas, model, lattice, defects=true)

## Check chebychev Metric
model = ForceXY(params_phys, params_num)
lattice = SquareLattice(L, periodic=true, metric="chebychev")
thetas = init_thetas(model, lattice, init="hightemp", q=1, r0=60, float_type=float_type, type=["source", "divergent"])
plot_theta(thetas, model, lattice)
tmax = 100
z = @elapsed update!(thetas, model, lattice, tmax)
plot_theta(thetas, model, lattice)
include("../src/defects_methods.jl")

plot_theta(thetas, model, lattice, defects=true)

## Checks whether the number of NN is plausible
i = 10;
j = 10;
model = XY(params_phys, params_num)
@btime get_neighbours(thetas, model, lattice, i, j)

model = ForcedXY(params_phys, params_num)
@btime get_neighbours(thetas, model, lattice, i, j)

model = SPP(params_phys, params_num)
@btime get_neighbours(thetas, model, lattice, i, j)

model = VisionConeXY(params_phys, params_num)
angle_neighbours = get_neighbours(thetas, model, lattice, i, j, true)
sum_influence_neighbours(thetas[i, j], angle_neighbours, model, lattice)
@btime get_neighbours(thetas, model, lattice, i, j, true)
@btime sum_influence_neighbours(thetas[i, j], angle_neighbours, model, lattice)
