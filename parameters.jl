Tf = float_type = Float32
xpu = "cpu" # "gpu" or "cpu"

## ------------- Lattice Parameters ------------- ## 
L =  100
b = Tf(0.5) # for TorusLattice (here η = 1)
eta = Tf(2) # for TorusHydroLattice (here b = 1). Set to 0 or Inf to knock out curvature and make it a cylinder

## ------------- Physics Parameters ------------- ## 
# Generic or multi-use
T = Tf(0.03) # angular temperature
rho = 1 # density of particles	
symmetry = "polar"
algo = "Langevin" # rule for collision!() for model = SPP : algo = "A", or type of XY model : algo = "MonteCarlo"/"MC" or "Langevin"

# Kuramoto 
Var = Tf(0.1) # for ForcedXY
omegas_corr = 0 # for ForcedXY ; has to be a diviser of L

# SPP 
A = Tf(0) # for SPP, probability to self propel in the direction of the velocity instead of isotropically 
propulsion = "polar"

# For Hydrodynamic Model, α is the soft constraint to enforce the unit length
alpha = Tf(100)

# Non-Reciprocal
sigma = Tf(0.2) # for DiscreteNRXY and NRXY. Recover XY model if sigma = 0
vision_cone = Tf(2π - (4pi*sigma)/(2sigma+pi)) # for VisionConeXY. Recover XY model if vision_cone = 2π
sigma_VC(sigma_VM) = 4pi*sigma_VM/(2sigma_VM+pi)
sigma_VM(sigma_VC) = pi*sigma_VC/(4pi-2sigma_VC)
function g(theta, uij, sigma) 
	return exp(sigma * cos(theta - uij))
	# return 1.0
end


## ------------- Numerical Parameters ------------- ## 
#dt = Tf(1E-5)
dt = Tf(0.1)
wrapsT = 16  # for GPU
block3D = (wrapsT,wrapsT,1) # for GPU
R = 10
grid3D = (Int(ceil(L/wrapsT)),Int(ceil(L/wrapsT)),R)

## ------------- Initialisation Parameters ------------- ## 
init = "spinwave_horizontal" # lowtemp, hightemp, single, pair, twist, spinwave, 
init = "hightemp" # lowtemp, hightemp, single, pair, twist, spinwave, 
q = 1
r0 = round(Int, L / 2)
theta0 = pi/2 # 0 < location < 2π of the defect + in the θ direction, for the TorusLattice
mu0 = 0.1 # for one defect only
mu1 = 0 # for single twisted defects only
decay_function = x->exp(-x/(L/20))
mu_plus, mu_minus, phi = pi/2, nothing, 0pi/4 # one of them has to be Nothing
nb_horizontal_waves = 1
nb_vertical_waves = 0

# ------------- Containers ------------- ## 
params_phys = Dict("L" => L,"b" => b, "eta" => eta, "T" => T,"alpha" => alpha, "Var" => Var, "A" => A, "rho" => rho, "vision_cone" => vision_cone, "g" => g, "sigma" => sigma, "symmetry" => symmetry, "propulsion" => propulsion, "algo" => algo, "omegas_corr" => omegas_corr)
params_num  = Dict("dt" => dt, "float_type" => float_type, "R" => R, "wrapsT" => wrapsT, "block3D" => block3D, "xpu" => xpu)
params_init = Dict("init" => init, "q" => q, "r0" => r0, "theta0" => theta0, "phi" => phi, "mu0" => mu0,"mu1" => mu1, "decay_function" => decay_function, "mu_minus" => mu_minus, "mu_plus" => mu_plus, "nb_horizontal_waves" => nb_horizontal_waves, "nb_vertical_waves" => nb_vertical_waves)

params = merge(params_num, params_phys, params_init)



