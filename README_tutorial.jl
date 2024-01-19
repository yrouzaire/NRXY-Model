using DrWatson;
@quickactivate "LatticeModels"; # tells Julia to load the correct package versions
include(srcdir("LatticeModels.jl")) # loads the entire code, might precompile some packages

using Plots, ColorSchemes, LaTeXStrings
# Load the plotting package. Can be very slow first time it is run, could take up to 30s - 1min
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot();

#= For the reader to do: Modify parameters.jl with the desired values.
For this tuto, take
    L = 100
    T = 0.1
    symmetry = "polar"
    rho = 1
    dt = 0.1
    init = "hightemp"
    The rest doesn't matter for this tutorial.
=#
include(srcdir("../parameters.jl")) # Load the parameters

## Declare your lattice: SquareLattice(L) or TriangularLattice(L).
lattice = SquareLattice(L)

## Declare your model: LangevinXY(params), MonteCarloXY(params), VisionConeXY(params), SPP(params) etc
model = LangevinXY(params)

## Initialisation of the theta field:
thetas = init_thetas(model, lattice, params)

## Temporal evolution
model.t # current time of the model = 0
# Update the model:
update!(thetas, model, lattice) # For one time step only
model.t # Now the current time of the model = dt
#= Indeed, the `update!` function also updates time.
Thus, don't forget to re-instantiate the model if you want to
restart a fresh simulation =#

duration = 15
z = @elapsed evolve!(thetas, model, lattice, duration) # For a given duration
model.t # Now the current time of the model = dt + duration
prinz(z) # prints the rounded runtime of this function call
tmax = 50
evolve!(thetas, model, lattice, tmax) # Until a given time
model.t # Now the current time of the model = tmax, no matter what was done before
# equivalently, if one wants to perform measurements over time
tmax2 = 80
while model.t < tmax2
    update!(thetas, model, lattice)
    # perform measurements
end

# Visualise the system
plot_thetas(thetas, model, lattice, defects=false)
plot_thetas(thetas, model, lattice, defects=true) # circle = +1 defect, triangle = -1 defect
#= WARNING: plotting the defects can take enormous time if they are to many (~ > 100).
It is recommended to first plot with `defects=false` (the default value) to evaluate visually
the number of defects. =#

## Some measurement on the system
polarOP, nematicOp = OP(thetas) # polar and nematic order parameters
C = corr(thetas, model, lattice) # the spatial correlation function C(r)
threshold = exp(-1)
ξ = corr_length(C, seuil=threshold) # the correlation length ξ (xi)
plot(xlabel="r", ylabel="C(r)")
plot!(1:L/2, remove_negative(C), axis=:log) # remove_negative
hline!([threshold], c=:grey, lw=0.5)
vline!([ξ], c=:grey, lw=0.5, line=:dash)
annotate!((1.5, 1.15threshold, text("y = threshold", 7, :left)))
annotate!((0.85ξ, 0.1, text("r = ξ", 7, :center, 90.0)))

## Create a movie
lattice = SquareLattice(L)
model = LangevinXY(params)
thetas = init_thetas(model, lattice, params_init)
every = 1;
tmax = 200;
transients = 50; # defects are not plotted before t ≥ transients (if defects=true)
saving_times = every:every:tmax
z = @elapsed animation = movies(thetas, model, lattice, defects=true, saving_times=saving_times, transients=transients)
prinz(z) # takes approximately 1 minute
filename = "films/my_first_movie_of_critical_XY_model.mp4"
mp4(animation, filename, fps=20) # creates the file in the given directory
# you have to open manually the .mp4 file

## Visualise the theta field at current (final) time
plot_thetas(thetas, model, lattice)

# Now if you want to visualize the arrows, you have to zoom in
# At random loc
xlocation, ylocation = (50, 50)
half_width_of_window = 10 # not too big because plotting the arrows can be damn slow
zoom_quiver(thetas, model, lattice, xlocation, ylocation, half_width_of_window)

# Where are the defects ?
vortices, antivortices = spot_defects(thetas, model, lattice)
nb_defects = length(vortices)

# Around a randomly chosen +1 defect
xlocation, ylocation = vortices[rand(1:nb_defects)][1:2]
zoom_quiver(thetas, model, lattice, xlocation, ylocation, half_width_of_window)

# Around a randomly chosen -1 defect
xlocation, ylocation = antivortices[rand(1:nb_defects)][1:2]
zoom_quiver(thetas, model, lattice, xlocation, ylocation, half_width_of_window)

## Explore the different initialisations
params_init["init"] = "hightemp" # equivalent of "disordered"
thetas = init_thetas(model, lattice, params_init)
plot_thetas(thetas, model, lattice)

params_init["init"] = "lowtemp" # equivalent of "ordered"
thetas = init_thetas(model, lattice, params_init)
plot_thetas(thetas, model, lattice) # everything is black, because the system is ordered (theta=0)

# For a +1 defect, one should provide the initial µ ∈ [0,2π] "mu0" (µ = 0 is a source, µ = ± π/2 is a vortex, µ = π is a sink)
params_init["init"] = "single";
params_init["q"] = 1;
mu0s = [0, π / 3, π, 5π / 3]
pp = []
for mu0 in mu0s
    params_init["mu0"] = mu0
    thetas = init_thetas(model, lattice, params_init)
    plott = zoom_quiver(thetas, model, lattice, 50, 50)
    title!(plott, "µ = $(round(mu0, digits=2))")
    push!(pp, plott)
end
plot(pp..., layout=(2, 2), size=(1000, 1000))

# The same goes for -1 defects
params_init["init"] = "single";
params_init["q"] = -1;
mu0s = [0, π / 3, π, 5π / 3]
pp = []
for mu0 in mu0s
    params_init["mu0"] = mu0
    thetas = init_thetas(model, lattice, params_init)
    plott = zoom_quiver(thetas, model, lattice, 50, 50)
    title!(plott, "µ = $(round(mu0, digits=2))")
    push!(pp, plott)
end
plot(pp..., layout=(2, 2), size=(1000, 1000))


## For pairs of defects, one has to 
params_init["init"] = "pair";
params_init["r0"] = round(Int, L / 2);
muplus, muminus, phii = pi / 2, nothing, pi/8 # one of the three MUST be nothing because only two of them are linearly independent
params_init["mu_plus"] = muplus;
params_init["mu_minus"] = muminus;
params_init["phi"] = phii
thetas = init_thetas(model, lattice, params_init)
plot_thetas(thetas, model, lattice)
