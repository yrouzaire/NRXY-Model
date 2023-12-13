# path_to_directory = "D:/Documents/Research/projects/LatticeModels"
path_to_directory = "Users/yrouzaire/Documents/Recherche/GitHub/Lattice-Models"
# cd(path_to_directory)
using DrWatson;
@quickactivate "LatticeModels";
include(srcdir("LatticeModels.jl"))
# using CUDA
using Plots, ColorSchemes, LaTeXStrings
gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot(rand([sin, cos, log, exp]))
include(srcdir("../parameters.jl"));
fn = "NRI/NRXY/investigation/"
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);

## ----------------- Data Generation (not for Hydrodyn on Torus) ----------------- ##
# ----------------- Data Generation (not for Hydrodyn on Torus) ----------------- ##
# ----------------- Data Generation (not for Hydrodyn on Torus) ----------------- ##
# ----------------- Data Generation (not for Hydrodyn on Torus) ----------------- ##

include(srcdir("../parameters.jl"));
lattice = SquareLattice(L)
model = NRXY_ExtraTerm(params)

# model = NRXY(params)
# lattice = TorusHydroLattice(L,eta)
# thetas,phis = init_thetas(model, lattice, params)
# saving_times = 0:1E-4:1E-2;

# lattice = SquareLattice(L)
# lattice = TorusLattice(L)
# model = NRXY(params)
# model = OUForcedXY(params)
# model = VisionConeXY(params)


update!(thetas, model, lattice)
thetas = init_thetas(model, lattice, params)
# thetas,phis = init_thetas(model, lattice, params)
saving_times = 0:10:1000;

thetas_saves = zeros(Float16, L, L, length(saving_times))

# thetas[1,:] .= Float32(NaN)
# thetas[:,1] .= Float32(NaN)
# thetas[end,:] .= Float32(NaN)
# thetas[:,end] .= Float32(NaN)

z = @elapsed for tt in each(saving_times)
    # if saving_times[tt] > saving_times[end]/6
    #     model.sigma = 0.0
    # end
    evolve!(thetas, model, lattice, saving_times[tt])
    # evolve!(thetas, phis, model, lattice, saving_times[tt])
    println("$(round(saving_times[tt]*100/saving_times[end],digits=1)) %")
    thetas_saves[:, :, tt] = thetas
end
prinz(z)
plot_thetas(thetas_saves[:,:,end], model, lattice, defects=false,colorbar=false,quiver=false)
# @save "data/etude_films/NRI_vision$(vision)_T$(T)_seed$(my_seed).jld2" vision T L seed=my_seed thetas_saves model lattice saving_times
# @load "data/etude_films/NRI_vision$(vision)_T$(T)_seed$(my_seed).jld2" vision T L seed thetas_saves model lattice saving_times

## ----------------- Movies ----------------- ##
# ----------------- Movies ----------------- ##

gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot_defects = false
find_mu = false
zm = @elapsed anim = @animate for tt in 1:round(Int,0.6*length(saving_times))
    println("$(round(tt*100/length(saving_times),digits=1)) %")
    thetass = Float32.(thetas_saves[:, :, tt])
    titre = "t = $(saving_times[tt])"
    if find_mu
        dft = DefectTracker(thetass, model, lattice)
        if (3 > number_active_defectsP(dft) > 0) 
            mup = mu2string(last_shape(dft.defectsP[1]))
            titre *= " µ+ = $mup"
        end
        if (3 > number_active_defectsN(dft) > 0) 
            mum = mu2string(last_shape(dft.defectsN[1]))
            titre *= " µ- = $mum"
        end
    end
        plot_thetas(thetass, model, lattice, size=(512, 512), 
            title=titre, defects=plot_defects)
end
prinz(zm)
prinz(zm+z)
#mp4(anim, path_to_directory*"/NRI/danpearce/PolarFluid_NR/films/" * filename * ".mp4", fps=30)
mp4(anim)
mp4(anim,pwd()*"/NRXY/films_paper/transverse_annihilation.mp4", fps=30)
## 

# ----------------- Data Generation for Hydrodyn on Torus ----------------- ##
# ----------------- Data Generation for Hydrodyn on Torus ----------------- ##
# ----------------- Data Generation for Hydrodyn on Torus ----------------- ##
# ----------------- Data Generation for Hydrodyn on Torus ----------------- ##

include(srcdir("../parameters.jl"));
model = HydroTest(params)
lattice = TorusHydroLattice(L,eta)
thetas,phis = init_thetas(model, lattice, params)
saving_times = 0:2E-3:2E-1;
thetas_saves = zeros(Float16, L, L, length(saving_times))
phis_saves = zeros(Float16, L, L, length(saving_times))

z = @elapsed for tt in each(saving_times)
    thetas_saves[:, :, tt] = thetas
    phis_saves[:, :, tt] = phis
    evolve!(thetas,phis, model, lattice, saving_times[tt])
    println("$(round(saving_times[tt]*100/saving_times[end],digits=1)) %")
end
prinz(z)
plot_thetas(Float32.(thetas_saves[:,:,end]),Float32.(phis_saves[:,:,end]), model, lattice)
# @save "data/etude_films/NRI_vision$(vision)_T$(T)_seed$(my_seed).jld2" vision T L seed=my_seed thetas_saves model lattice saving_times
# @load "data/etude_films/NRI_vision$(vision)_T$(T)_seed$(my_seed).jld2" vision T L seed thetas_saves model lattice saving_times

# ----------------- Movies ----------------- ##
# ----------------- Movies ----------------- ##

gr(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
zm = @elapsed anim = @animate for tt in 1:round(Int,length(saving_times))
    println("$(round(tt*100/length(saving_times),digits=1)) %")
    thetass = Float32.(thetas_saves[:, :, tt])
    phiss = Float32.(phis_saves[:, :, tt])
    titre = "t = $(saving_times[tt])"
    
    plot_thetas(thetass,phiss, model, lattice,size=(512,512), vertical=false)
end
prinz(zm)
prinz(zm+z)

filename = "$(init)_eta$(eta)_L$(L)_T$(T)_tmax$(saving_times[end])_dt$(dt)_α$(alpha)"
mp4(anim, path_to_directory*"/NRI/danpearce/PolarFluid_NR/films/C" * filename * ".mp4", fps=30)
pwd()
& 
## -------------- Old scripts -------------- ##
## -------------- Old scripts -------------- ##
## -------------- Old scripts -------------- ##
## -------------- Old scripts -------------- ##

## movie zoomed
# locs = [(130,0)]
# z = @elapsed for (locx, locy) in locs
#     println("locx = $locx, locy = $locy")
#     anim_zoom = @animate for tt in 1:length(saving_times)
#         zoom_quiver(thetas_saves[:, :, tt], model, lattice, locx, locy, 20, size=(1024, 1024))
#     end
#     # mp4(anim_zoom, fn * "one_defect/mutation-1_zoom($locx,$locy).mp4", fps=20)
#     mp4(anim_zoom, fn * "system/FXY_Triangular_zoom.mp4", fps=20)
#     # mp4(anim_zoom, fn * "system/FXY_Square_zoom.mp4", fps=20)
# end
# prinz(z)

## ------------------ Frame by Frame Movies ------------------
# pyplot(box=true, fontfamily="sans-serif",fmt=:png, label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
# include(srcdir("../parameters.jl"));
# lattice = SquareLattice(L)
# lattice = TriangularLattice(L)
# model = XY(params)
# thetas = init_thetas(model, lattice, params_init=params_init)

# z = @elapsed anim_fbf = @animate for i in 1:500
#     evolve!(thetas, model, lattice,5model.dt)
#     params_init["q"] == +1 ? mm = 1 : mm = 2 
#     i,j = 25,25#spot_defects(thetas, model, lattice)[mm][1][1:2] # location of the defect
#     zoom_quiver(thetas, model, lattice, i, j, WINDOW,size=(512, 512))
#     # plot_thetas(thetas, model, lattice, defects=false,colorbar=false,size=(512, 512))
# end
# filepath = "NRI/NRXY/investigation/one_defect/movies_frame_by_frame/negative_defect_XY_Squ.mp4"
# mp4(anim_fbf, filepath, fps=30)


# zoom_quiver(thetas, model, lattice, i, j, WINDOW,size=(600, 600))
# zoom_quiver(thetas, model, lattice, i, j, WINDOW,size=(400, 400))

