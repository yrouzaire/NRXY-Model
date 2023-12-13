cd("D:/Documents/Research/projects/LatticeModels")
    using DrWatson ; @quickactivate "LatticeModels"
    include(srcdir("LatticeModels.jl"));
    using Plots,ColorSchemes,LaTeXStrings
    gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    include(srcdir("../parameters.jl"));
##
CHARGE = 1
@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP1.jld2")
base_datasetP = base_dataset
@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µN1.jld2")
base_datasetN = base_dataset;

model = XY(params)
lattice = SquareLattice(L)
lattice = TriangularLattice(L)
## Visually check the dataset
ind = rand(1:length(mus))
    p=plot_thetas(base_datasetP[:,:,ind],model,lattice)
    display_quiver!(p,base_datasetP[:,:,ind], lattice,window=WINDOW)
    title!("µ = $(mus[ind])")

##
ind = rand(1:length(mus))
    p=plot_thetas(base_datasetN[:,:,ind],model,lattice)
    display_quiver!(p,base_datasetN[:,:,ind], lattice,window=WINDOW)
    title!("µ = $(mus[ind])")

##

 
## Test noiseless
inferredP  = zeros(length(mus))
    inferredN = zeros(length(mus))
    decay = true
    for ind in 1:64
        inferredP[ind] = infer_mu(base_datasetP[:,:,ind],model,lattice,8,8,CHARGE)
        inferredN[ind] = infer_mu(base_datasetN[:,:,ind],model,lattice,8,8,-CHARGE)
    end
    p1 = plot(xlabel="True µ",ylabel="Inferred µ",legend=:topleft,title="Noiseless")
    plot!(mus,inferredP,m=true,ms=2,c=:red,label="q=1")
    plot!(mus,inferredN,m=true,ms=2,c=:green,label="q=-1")
    plot!(x->x,c=:black)

## Test with noise
lattice = TriangularLattice(W21,periodic=false)
lattice = SquareLattice(W21,periodic=false)
model = XY(params)
model.T = 0.2
inferredP  = zeros(length(mus))
    inferredN = zeros(length(mus))
    for ind in 1:64
        thetas_P = base_datasetP[:,:,ind]
        evolve!(thetas_P, model, lattice,1)
        inferredP[ind] = infer_mu(thetas_P,model,lattice,8,8,CHARGE)
        
        thetas_N = base_datasetN[:,:,ind]
        evolve!(thetas_N, model, lattice,1)
        inferredN[ind] = infer_mu(thetas_N,model,lattice,8,8,-CHARGE)
    end
    p2 = plot(xlabel="True µ",ylabel="Inferred µ",legend=false,title="Noisy dyn (T = $(model.T))")
    plot!(mus,inferredP,line=false,m=true,ms=2,c=:red,label="q=1")
    plot!(mus,inferredN,line=false,m=true,ms=2,c=:green,label="q=-1")
    plot!(x->x,c=:black)

# plot(p1,p2,size=(400,800),layout=(2,1))
# savefig("figures/procedure_infer_mu_polar.png")
