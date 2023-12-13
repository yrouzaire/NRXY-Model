# LatticeModels

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> LatticeModels

It is authored by Ylann Rouzaire.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

------------------------------------------------------------------------------
# General philosophy of this package
This code aims at simulating the temporal evolution of a 2D vector-field,
parametrized in space by a scalar $\theta(x,y) \in [0,2\pi[$.

Throughout the entire project, we discretize this $\theta$ field onto a discrete 2D lattice.
  - Pros: the numerical simulations are much faster and easily controlled, in particular
the topological defect detection/tracking.
  - Cons: discrete lattices tend to induce undesired discretization artifacts, especially the
usual `SquareLattice`.



### Step-by-Step Tutorial
While reading this document, run `README_tutorial.jl` in parallel for a hands-on
discovery of the main tools of the code.

# Important files to get started

Here we comment the most important methods in each of these files.

## parameters.jl
Contains physical, numerical and initialisation parameters.

### Physical Parameters
  - $L$ is the (integer) size of the 2D $L\times L$ lattice. It thus contains $L^2$ agents/spins/particles/rotors/oscillators (all synonyms),
  such that the runtime is of minimum complexity $\mathcal{O}(L^2)$.
  - $T\ge 0$ is the temperature of the bath related to the angle diffusion.
  - `symmetry` is the symmetry of the interaction and can be either "nematic" or "polar".
  Polar means that 2 interacting particles $i$ and $j$ want to point in the same direction, and the force $\propto \sin(\theta_j - \theta_i)$
  Nematic means that 2 interacting particles $i$ and $j$ want to point in the same direction but without the notion of head/tail, and the force $\propto \sin(2(\theta_j - \theta_i))$
  - `propulsion` is the symmetry of propulsion of Self Propelled Particles (SPP model) and can only be "polar" for now. Polar propulsion: self-propelled rods. Nematic propulsion (reverses from times to times): active nematics
  - `Var` $\ge 0$ is the variance of the distribution of intrinsic/internal frequency for the model ForcedXY.
  - $A\ge 0$ is the value of the coupling between orientation and propulsion. $A = 0$: no coupling, $A \to \infty$: high coupling, i.e you are sure you move in the direction of your arrow (if polar propulsion)
  - $0 <\rho \le 1$ is the density. If $\rho = 0.9$, 10% of the sites will be empty.
  - `algo` is used independently (but never at the same time) for two cases.
    - If the model is of type XY, `algo = "Langevin"` will make its temporal evolution follow a Langevin dynamics.
    - If the model is of type XY, `algo = "MonteCarlo"` (alias `"MC"`)  will make its temporal evolution follow a Metropolis MonteCarlo dynamics.
    - If the model is `SPP`, set `algo = "A"` (only one functional)

### Numerical Parameters
  - `dt` is the timestep used in the evolution of Langevin-like models (ie. `ForcedXY`,  `LangevinXY`, `LangevinNematicXY`, `ContinuousNRXY`)
  - `float_type = Float32` is the recommanded Float type. Float64 takes more memory. Float16 is not a native Float type so operations are usually slower.
  - `tmax` is the maximum time for temporal evolution.
  - `transients` is the transient time, before which measurements on the system are not conducted.
  - `every` is the interval at which measurements are made, in the case you want linearly spaced data points.

### Initialisation Parameters
  - `init` is the type of initialisation. Can be:
    - "ordered" or "lowtemp" for an ordered initial configuration ( $T=0$ ).
    - "disordered" or "hightemp" for a disordered initial configuration ( $T=\infty$ ).
    - "single" or "isolated" for a manually created vortex (in this case, the charge $q$ and its type `type1defect` should be provided)
    - "pair" for a manually created vortex pair (in this case, the absolute value of the charge $q$, the interdefect separation `r0` (integer) and its type `type2defect` should be provided)
    - "2pair" or "2pairs" for TWO manually created vortex pairs (in cross shape).

## lattice.jl
All the methods to create and handle $L\times L$ 2D lattices: `SquareLattice` (4 nearest neighbours) and `TriangularLattice` (6 nearest neighbours) are implemented.
Motivation for implementing lattices different from the standard square one:
  - Discrete lattices may induce artifacts. The `TriangularLattice` is the closer to a continuum description in the thermodynamic limit $N\to\infty$.
  - Geometric frustration is usually observed on triangular or honeycomb lattices.

### Attributes
  - $L$ (integer) the linear size of the lattice. For now, only $L_x = L_y$ lattices are supported.
  - `periodic=true` stands for Periodic Boundary Conditions. `periodic=false` stands for Free Boundary Conditions.
  - `single=true` is unused for now, don't care about it.

### Implementation
Nothing exotic to say about the `SquareLattice`, the implementation is straightforward.
On the other hand, he TriangularLattice is somewhat more technical, so that it could be represented by a usual Matrix.
The trick is to modify the neighbourhood depending on the parity of the row of the considered site, cf. `get_neighbours(thetas,model,lattice::TriangularLattice,i,j)` in the file `core_methods.jl`

## models.jl
The file where all the models are defined.
Each instantiable model (i.e. `ForcedXY`) is a subtype of `AbstractModel`.
Note on the XY model, the base model from which I depart from equilibrium:
Two dynamics are available:
  - the overdamped Langevin dynamics (`LangevinXY`), the default.
  - the Metropolis Monte Carlo dynamics (`MonteCarloXY`)
You choose one or the other specifying `algo = "MonteCarlo"` or `algo = "Langevin"` in *parameters.jl* and then running `XY(params)`,
or by specifically running `LangevinXY(params)` / `MonteCarloXY(params)`
Both are subtypes of `AbstractXYModel`, which itself is a subtype of `AbstractModel`.
## init_visu.jl

## core_methods.jl
Basically implements the temporal evolution of the system, making extensive use of multiple dispatch.
The `get_neighbours` method is central, and dynamically adapts to the lattice type.
The `update!` method is defined for every single model.

Note: `update!(thetas,model,lattice)` updates each spin once. This is by definition what I define as *one* timestep $dt$.

# Other files
`auxiliary.jl` contains useful functions unrelated to the specifics of `LatticeModels`

`measurements.jl` contains functions to measure physical quantities in the system: polar/nematic order (=magnetisation), correlation functions and correlation lengths etc

# About topological defects (detection, tracking in-vivo)
The code is located in the defects_methods.jl file.
I will for now only comment on the first 3 methods, used to detect defects.
  - `arclength` returns the (signed) distance, in radians (hence the name) between two angles.
  - `get_vorticity` returns the vorticity of a given site $i,j$ on the lattice, following the standard procedure (well explained in Fig. 2.34 of *Advanced Statistical Physics: 2. Phase Transitions* of Leticia Cugliandolo)
  - `spot_defects` scans the system site by site and applies `get_vorticity` so localize defects. Often, the algo detects 2/3 (also depends on the lattice type) numerical defects for a single physical defect.
  This issue seems to be only dealt with *a posteriori*, this is what the routine `merge_duplicates` does. I won't explain the `find_types` routine for now.
