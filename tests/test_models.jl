include("../src/models.jl");
cols = cgrad([:black, :blue, :green, :orange, :red, :black])

# Physical Parameters
include(srcdir("../parameters.jl"));

model = ForcedXY(params)
heatmap(model.thetas')

model = XY(params)
model = MonteCarloXY(params)
model = MonteCarloNematicXY(params)
model = LangevinXY(params)
model = LangevinNematicXY(params)
model = VisionConeXY(params)
model = DiscreteNRXY(params)
model = NRXY(params)
model = SPP(params)


## Test NRXY
include(srcdir("../parameters.jl"));
# using Random #Random.seed!(1)
params["dt"] = 1E-1
params["T"] = 0
params["sigma"] = 3E-1
model = NRXY(params)
thetas = init_thetas(model, lattice, params)
# thetas[1, :] .= Float32(NaN)
# thetas[:, 1] .= Float32(NaN)
# thetas[end, :] .= Float32(NaN)
# thetas[:, end] .= Float32(NaN)
lattice = TriangularLattice(L)
evolve!(thetas, model, lattice, 1E1)
evolve!(thetas, model, lattice, 2E2)
# z = @elapsed evolve!(thetas, model, lattice, 1E3);
# prinz(z);
plot_thetas(thetas, model, lattice)
# number_defects(thetas, model, lattice)

## Test HydroVector
include(srcdir("../parameters.jl"));
# using Random #Random.seed!(1)
model = HydroVector(params)
lattice = TorusHydroLattice(L, eta)
thetas, phis = init_thetas(model, lattice, params)
evolve!(thetas, phis, model, lattice, 1E-4)
evolve!(thetas, phis, model, lattice, 1E-3)
z = @elapsed evolve!(thetas, phis, model, lattice, 1E-1);
prinz(z);
plot_thetas(thetas, phis, model, lattice);title!("dt = $dt, α = $alpha, η = $eta")


## Test HydroTest
include(srcdir("../parameters.jl"));
# using Random #Random.seed!(1)
model = HydroTest(params)
lattice = TorusHydroLattice(L, eta)
thetas, phis = init_thetas(model, lattice, params)
evolve!(thetas, phis, model, lattice, 1E-4)
evolve!(thetas, phis, model, lattice, 1E-3)
plot_thetas(thetas, phis, model, lattice);
title!("dt = $dt, α = $alpha, η = $eta")

##
z = @elapsed evolve!(thetas, phis, model, lattice, 4E-2);prinz(z);
plot_thetas(thetas, phis, model, lattice); title!("dt = $dt, α = $alpha, η = $eta")
