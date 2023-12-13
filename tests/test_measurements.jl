cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## With TriangularLattice
include(srcdir("../parameters.jl"));
    model = XY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
    thetas = init_thetas(model, lattice,params)
# plot_thetas(thetas,model,lattice)

# OP(thetas)

evolve!(thetas,model,lattice,1000)
ctri = corr(thetas,model,lattice)
plot(ctri)
corr_length(ctri)

## With SquareLattice
include(srcdir("../parameters.jl"));
    model = XY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
# plot_thetas(thetas,model,lattice)

# OP(thetas)
update!(thetas,model,lattice,10)
csqu = corr(thetas,model,lattice)
plot!(csqu)

## Test Energy on single defects
include(srcdir("../parameters.jl"));
params_init["init"] = "single"
params_init["mu0"] = 0
params["sigma"] = 1
model = NRXY(params)
# lattice = TriangularLattice(L)
lattice = SquareLattice(L)
thetas = init_thetas(model, lattice, params_init=params_init)
plot_thetas(thetas, model, lattice)
energy(thetas, model, lattice)
energy(thetas, model, lattice)

## Test tmax = 1E4
include(srcdir("../parameters.jl"));

## Energy on TorusLattice
include(srcdir("../parameters.jl"));
model = NRXY(params)
lattice = TorusLattice(L)
thetas = init_thetas(model, lattice, params)
energy(thetas, model, lattice)
energy(thetas, model, lattice)
@btime energy($thetas, $model, $lattice)

## -------------- Correlation length with FFT -------------- ##
## -------------- Correlation length with FFT -------------- ##
## -------------- Correlation length with FFT -------------- ##
## -------------- Correlation length with FFT -------------- ##
using FFTW
include(srcdir("../parameters.jl"));
model = XY(params)
lattice = SquareLattice(L)
lattice = TriangularLattice(L)
thetas = init_thetas(model, lattice, params)
evolve!(thetas, model, lattice, 2E3)
plot_thetas(thetas, model, lattice)


@time c_classique = corr(thetas, model, lattice)
@btime c_fft = corr_fft(thetas)
plot(c_classique, uaxis=:log)
plot!(c_fft)

Ls = collect(50:50:1000)
R = 10
zs = zeros(length(Ls), R)
for i in each(Ls)
    L = Ls[i]
    thetas = rand(L,L)
    for j in 1:R
        z = @elapsed corr_fft(thetas)
        zs[i,j] = z
    end
end # takes 7s for R=10 and Lmax = 1000
prinz(sum(zs))
zs_avg = mean(zs,dims=2)[:,1]
##
plot(Ls,zs_avg,m=true, axis=:log )
plot!(Ls,x->1.5E-8x^2*log(x))