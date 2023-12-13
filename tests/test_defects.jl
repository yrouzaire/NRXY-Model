using DrWatson;
@quickactivate "LatticeModels";
include(srcdir("LatticeModels.jl"));
using BenchmarkTools, Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();

## ----------- Tests for defect spotting without dynamics with PBC ----------- ##
include(srcdir("../parameters.jl"));
params["symmetry"] = "polar"
params["init"] = "pair"
params["q"] = 1
params["r0"] = 40
params["phi"] = rand()*2pi
params["mu_plus"] = nothing
params["mu_minus"] = 1

model = XY(params)
lattice = SquareLattice(L)
lattice = TriangularLattice(L)

thetas = init_thetas(model, lattice, params)

evolve!(thetas, model, lattice, 1)

ndf = number_defects(thetas, model, lattice)
dfp = spot_defects(thetas, model, lattice)[1][1]
dfm = spot_defects(thetas, model, lattice)[2][1]
plot_thetas(thetas, model, lattice, defects=true,colorbar=false,quiver=false)
title!("n=$ndf, µ = $(round(dfp[4],digits=2)), $(round(dfm[4],digits=2))")

## ----------- Tests for defect spotting without dynamics without PBC ----------- ##
include(srcdir("../parameters.jl"));
params["symmetry"] = "polar"
params["init"] = "single"
params["q"] = 1
params["mu0"] = 5
# params["mu0"] = 

model = XY(params)
# lattice = SquareLattice(L)
lattice = TriangularLattice(L)

thetas = init_thetas(model, lattice, params)

evolve!(thetas, model, lattice, 10)

ndf = number_defects(thetas, model, lattice)
p=plot_thetas(thetas, model, lattice, quiver=true,defects=true,colorbar=false)

try
	df = spot_defects(thetas, model, lattice)[params["q"] == 1 ? 1 : 2][1]
	title!("$(df[1:2]),$(df[4])")
catch ; end 
p

# xyqµ_plaquette(thetas,model,lattice,5,4,2pi)
# spot_defects(thetas, model, lattice)

## 
include(srcdir("../parameters.jl"));
model = NRXY(params)
lattice = TriangularLattice(L)
lattice = SquareLattice(L)
params_init["init"] = "single"
params_init["q"] = 1
params_init["mu0"] = -0.2
thetas = init_thetas(model, lattice, params_init=params_init)
plot_thetas(thetas, model, lattice, defects=false)
zoom_quiver(thetas, model, lattice, round(Int, L / 2), round(Int, L / 2), 4)

dft = DefectTracker(thetas, model, lattice)
dft.defectsP[1].type[1]
# dft.defectsN[1].type[1]
# dft = DefectTracker(thetas, model, lattice)
# dft.defectsP[1].type[1]

# infer_mu(zoom(thetas, lattice, 50,50)[2],q=1)
# thetaszoom = zoom(thetas, lattice, 50,50)[2]
# X_noisy_reshaped = reshape(thetaszoom,(W21,W21,1,:))
# recon = NN_test(provide_div_rot(X_noisy_reshaped))[:,:,1,1]
# infer_mu(recon,q=1, decay=false)

## Calibration of the inferrence procedure
include(srcdir("../parameters.jl"));
params["L"] = 53
params["symmetry"] = "polar"
params["init"] = "single"
params["q"] = 1
params["sigma"] = 0
mus = 0:pi/16:2pi-pi/16

R = 2
inferred = zeros(length(mus),R)
for i in each(mus)
	params["mu0"] = mus[i]
	Threads.@threads for r in 1:R
		model = XY(params)
		# model = NRXY(params)
		# lattice = SquareLattice(L)
		lattice = TriangularLattice(L)
		thetas = init_thetas(model, lattice, params)
		evolve!(thetas, model, lattice, 10)
		
		vp, vm = spot_defects(thetas, model, lattice)
		inferred[i,r] = params["q"] == 1 ? vp[1][4] : vm[1][4]
	end
end
inferred_avg = mean(inferred,dims=2)[:,1]
inferred_std = std(inferred,dims=2)[:,1]/sqrt(R)

plot()
plot!(mus,x->x,lw=2)
plot!(mus, inferred_avg, lw=1,line=:solid,c=:black)
# plot!(mus, mod.(inferred_avg.+0.02,2pi),rib=inferred_std, lw=1,line=:solid)
# plot!(mus,x->x+0.15)
# plot!(mus,x->x-0.15)
xlims!(0,2pi)
ylims!(0,2pi)

## Test acceleration and velocity of defects
include(srcdir("../parameters.jl"));
model = NRXY(params)
lattice = TriangularLattice(L, periodic=true)
thetas = init_thetas(model, lattice, params)
# plot_thetas(thetas, model, lattice, defects=false)

p=plot(xlabel="µ",ylabel="P(µ)")
times = 0:10:50
for tt in each(times)
	evolve!(thetas, model, lattice, times[tt])
	dft = DefectTracker(thetas, model, lattice)
	histogram!(get_musP(dft),normalize=true,bins=50,yaxis=:log)
	# plot_thetas(thetas, model, lattice, defects=false,title="t = $(model.t)")
end 
p 
ylims!(1e-2,1.1)

## 
include(srcdir("../parameters.jl"));
params["init"] = "pair"
model = NRXY(params)
model = VisionConeXY(params)
lattice = TriangularLattice(L, periodic=true)
thetas = init_thetas(model, lattice, params)
plot_thetas(thetas, model, lattice, defects=false)
##
evolve!(thetas, model, lattice, model.t + 100)
plot_thetas(thetas, model, lattice, defects=false,title="t = $(model.t)")


## 
using DifferentialEquations
f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

using Plots
plot(sol, linewidth = 5, title = "Solution to the linear ODE with a thick line",
     xaxis = "Time (t)", yaxis = "u(t) (in μm)", label = "My Thick Line!") # legend=false
plot!(sol.t, t -> 0.5 * exp(1.01t), lw = 3, ls = :dash, label = "True Solution!")