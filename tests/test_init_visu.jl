cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson;
@quickactivate "LatticeModels";
include(srcdir("LatticeModels.jl"))
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();

## Test pair on Torus with theta0 
include(srcdir("../parameters.jl"));
params["init"] = "pair"
params["theta0"] = -0.1
params["r0"] = round(Int, L / 2)
params["mu_plus"] = 0
model = NRXY(params)
lattice = TorusLattice(L)
# lattice = SquareLattice(L)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice, defects=true,colorbar=false,quiver=false)



## Test init single 
include(srcdir("../parameters.jl"));
L = 28
params["L"] = L
params["init"] = "single_twisted"
params["q"] = 1
params["mu0"] = 0
params["mu1"] = pi/2
params["decay_function"] = x->exp(-x^2/50)

model = NRXY(params)
# model = XY(params)
lattice = TriangularLattice(L,periodic=false)
# lattice = SquareLattice(L)
thetas = init_thetas(model, lattice, params)
evolve!(thetas, model, lattice, 500)
plot_thetas(thetas, model, lattice, defects=true,colorbar=false,quiver=true,size=(800,800))   


## Test init single 
include(srcdir("../parameters.jl"));
L = 14
params["L"] = L
params["init"] = "single"
params["q"] = 1
params["mu0"] = 0
model = NRXY(params)
lattice = TriangularLattice(L,periodic=false)
# lattice = SquareLattice(L)
thetas = init_thetas(model, lattice, params)
# evolve!(thetas, model, lattice, 50)
plot_thetas(thetas, model, lattice, defects=true,colorbar=false,quiver=true,size=(800,800))   
xlims!(0.5,L+1)
ylims!(0.,L+1)

## Test init spinwave
include(srcdir("../parameters.jl"));
params["init"] = "spinwave"
nb_horizontal_waves = 1; params["nb_horizontal_waves"] = nb_horizontal_waves
nb_vertical_waves = 1; params["nb_vertical_waves"] = nb_vertical_waves
model = NRXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice, defects=false,colorbar=false,quiver=false)

cr = corr(thetas, model, lattice)
plot(1:100,remove_negative(cr))
evolve!(thetas, model, lattice,model.t + 100)
cr = corr(thetas, model, lattice)
plot(1:100,remove_negative(cr),axis=:log)

## Theoretical predictions
θ_(x) = 2*atanh(tan(π*(x-1/2)))
θ_(-0.2)
θ(x,t) = 2*atan(tanh(θ_(x)+ t))
x = 0:0.01:1
tan(-pi/2)
plot(x,θ.(x,0),label="t=0")
## Test flag (sharp and smooth) initialisations
include(srcdir("../parameters.jl"));
model = XY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice)

## Tests Movies
# gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5)
model = VisionConeXY(params)
# model = NRXY(params)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice)
saving_times = 0:2:200;
transients = saving_times[end] / 2;
# update!(thetas,model,lattice)
# plot_theta(thetas,model,lattice,defects=true)
z = @elapsed anim_extensile = movies(thetas, model, lattice, defects=true, saving_times=saving_times, transients=transients)
prinz(z)
mp4(anim_extensile, "NRI/NRXY/investigation/system/smooth_flag_VC.mp4")

## Basic Tests
thetas = init_thetas(model, lattice, init="isolated", q=1, type="source")
plot_theta(thetas, model, lattice)
&

thetas = init_thetas(model, lattice, init="isolated", q=-1, type="convergent")
plot_theta(thetas, model, lattice)
&

thetas = init_thetas(model, lattice, init="isolated", q=1 / 2, type="source")
plot_theta(thetas, model, lattice)
&

thetas = init_thetas(model, lattice, init="pair", q=1, r0=60, type=["source", "divergent"])
plot_theta(thetas, model, lattice)
&

thetas = init_thetas(model, lattice, init="2pair", q=1, r0=60, type=["source", "divergent"])
plot_theta(thetas, model, lattice)
&

thetas = init_thetas(model, lattice, init="pair", q=1 / 2, r0=60, type=["source", "divergent"])
plot_theta(thetas, model, lattice)
&

thetas = init_thetas(model, lattice, init="2pair", q=1 / 2, r0=60, type=["source", "divergent"])
plot_theta(thetas, model, lattice)


## Test init pair with angle phi
include(srcdir("../parameters.jl"));
model = XY(params)
lattice = SquareLattice(L)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice, defects=true)
dft = DefectTracker(thetas, model, lattice, find_type=true)
if number_active_defects(dft) > 0
    tp = string(round(last_type(dft.defectsP[1]), digits=2))
    tm = string(round(last_type(dft.defectsN[1]), digits=2))
    titre = L"µ_{+}=" * tp * L" ; µ_{-}=" * tm
else
    titre = ""
end
zoom_quiver(thetas, model, lattice, 50, 50, 12)
title!(titre)


## Test the zoom with periodic lattice
include(srcdir("../parameters.jl"));

model = XY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model, lattice, params)
p = plot_thetas(thetas, model, lattice, colorbar=false)
evolve!(thetas, model, lattice, 50)
p = plot_thetas(thetas, model, lattice, colorbar=false)

# window = 12
# thetas_zoom = zoom(thetas, lattice, round(Int, L / 2), 55, window)[2]
# p = plot_thetas(thetas_zoom, model, lattice, colorbar=false)

window = 10
zoom_quiver(thetas, model, lattice, round(Int, L / 2), round(Int, L / 2), window)

## Tests plotting 1D chains
lattice = Chain1D(L)
model = PropagationForcedXY(params)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice)

## Test plotting routines with quiver integrated
include(srcdir("../parameters.jl"));
lattice = TriangularLattice(L)
model = DiscreteNRXY(params)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice, quiver=true, force_quiver=false)

dft = DefectTracker(thetas, model, lattice, find_type=true)
dft.defectsP[1].type[1]
dft.defectsN[1].type[1]
tmax = 20;
every = 0.25;
update_and_track!(thetas, model, lattice, dft, tmax, every, find_type=false, verbose=false)
plot_thetas(thetas, model, lattice, quiver=true, force_quiver=false)
plot_mus(dft, lattice.L, vel=true, acc=true, over=10)

##
include("../parameters.jl");
model = DiscreteNRXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice)



## For a single defect
include(srcdir("../parameters.jl"));
params_init["init"] = "single";
params_init["mu0"] = -pi / 2
model = NRXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice)

## For pairs of defects, one has to 
include(srcdir("../parameters.jl"));
params_init["init"] = "pair";
params_init["r0"] = round(Int, L / 2)
muplus, muminus, phii =  pi/2, nothing, pi/2 # one of the three MUST be nothing because only two of them are linearly independent
params_init["mu_plus"] = muplus;
params_init["mu_minus"] = muminus;
params_init["phi"] = phii

model = NRXY(params)
lattice = SquareLattice(L)
thetas = init_thetas(model, lattice, params)
# dft = DefectTracker(thetas, model, lattice, find_type=true)
# titre = L"µ_{+}=" * string(round(last_type(dft.defectsP[1]), digits=2)) * L" ; µ_{-}=" * string(round(last_type(dft.defectsN[1]), digits=2))
plot_thetas(thetas, model, lattice, colorbar=false)
# savefig("test.png")
# title!(titre)

## Visualisation arrows single defect
include(srcdir("../parameters.jl"));
params_init["init"] = "single";
params_init["mu0"] = -0pi / 2
params_init["q"] = 1 
model = NRXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice, quiver=true, force_quiver=false)



## Visualisation arrows pair defect
include(srcdir("../parameters.jl"));
params_init["init"] = "pair";
params_init["q"] = 1 
mup, mum, fi = nothing, -pi/2, 0pi/2# one of the three MUST be nothing because only two of them are linearly independent
params_init["mu_plus"] = mup 
params_init["mu_minus"] = mum
params_init["phi"] = fi

model = NRXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model, lattice, params)
# dft = DefectTracker(thetas, model, lattice, find_type=true)
# titre = L"µ_{+}=" * string(round(last_type(dft.defectsP[1]), digits=2)) * L" ; µ_{-}=" * string(round(last_type(dft.defectsN[1]), digits=2))
plot_thetas(thetas, model, lattice, quiver=true, force_quiver=false, colorbar=true, defects=true)
# title!(titre)
# round.(thetas, digits=1)
# aa
# thetas[8,:]


## Visualisation pair init on Torus (TorusHydroLattice)
include(srcdir("../parameters.jl"));
params["init"] = "pair";
params["eta"] = 2;
eta = Tf(params["eta"])
params["q"] = 1
params["r0"] = round(Int, L / 2)
params["mu_plus"] = 0
params["theta0"] = 0.045pi
params["sigma"] = 0

model = HydroTest(params)
lattice = TorusHydroLattice(L, eta)
thetas, phis = init_thetas(model, lattice, params)
evolve!(thetas, phis, model, lattice, 3E-2)
plot_thetas(thetas, phis, model, lattice)
# psis, S2 =get_psis_S2(thetas, phis, lattice)
# dummy_lattice = SquareLattice(L)
# plot_thetas(psis, model, dummy_lattice, defects=false, colorbar=false, quiver=false)