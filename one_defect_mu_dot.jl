# path_to_directory = "D:/Documents/Research/projects/LatticeModels"
path_to_directory = "Users/yrouzaire/Documents/Recherche/GitHub/Lattice-Models"
# cd(path_to_directory)
using DrWatson;
@quickactivate "LatticeModels";
include(srcdir("LatticeModels.jl"))
using Plots, ColorSchemes, LaTeXStrings
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5);
plot(rand([sin, cos, log, exp]))
include(srcdir("../parameters.jl"));
fn = "NRXY/investigation/one_defect/"

## ----------------- µ(t) during a torsion (+) ----------------- ##
include(srcdir("../parameters.jl"));
L = 32
params["L"] = L
params["dt"] = 1E-1
params["T"] = 0E-5
params["init"] = "single"
params["q"] = 1
params["r0"] = Int(round(L/2))
sigmas = collect(0.05:0.05:0.5)
# sigmas = [0.25, 0.3, 0.35, 0.4]
# sigmas = collect(0.1:0.1:0.5)
# sigmas = [0.1, 0.2, 0.5]
# mu0s = [0.15, 2pi-0.15]
mu0s = vcat(collect(0.1:0.05:0.8),collect(2pi-0.8:0.05:2pi-0.1))
tmax = 3000 ; times = 0:1:tmax
R = 1
M = length(sigmas)*length(mu0s)
m = 0

mus = NaN*zeros(length(sigmas), length(mu0s), length(times),R) 
z = @elapsed for i in each(sigmas)
    params["sigma"] = sigmas[i]
    params["vision_cone"] = 2pi - 4sigmas[i]
    for j in each(mu0s)
        m += 1
        params["mu0"] = mu0s[j] #+ 0.005randn()
        # params["mu_plus"] = mu0s[j] + randn() * 0.01
        println("σ = $(sigmas[i]), µ0 = $(round(mu0s[j],digits=2))" * " m/M = $m/$M")
        # println("VC = = $(2pi - sigmas[i]), µ0 = $(mu0s[j])")
        # for r in 1:R
        Threads.@threads for r in 1:R
            try 
                model = NRXY(params)
                # model = HydrodynNRXY(params)
                # model = XY(params)
                # model = VisionConeXY(params)
                # lattice = SquareLattice(L)
                lattice = TriangularLattice(L)
                thetas = init_thetas(model, lattice, params)
                
                # thetas[1,:] .= Float32(NaN)
                # thetas[:,1] .= Float32(NaN)
                # thetas[end,:] .= Float32(NaN)
                # thetas[:,end] .= Float32(NaN)

                dft = DefectTracker(thetas, model, lattice)
                update_and_track!(thetas, model, lattice,dft,times,plot=false,verbose=false)
                # for tt in each(times)
                #     evolve!(thetas,model,lattice,times[tt])
                #     update_DefectTracker!(dft, thetas, model, lattice)
                #     df = params["q"] > 0 ? dft.defectsP[1] : dft.defectsN[1]
                #     mus[i, j, tt, r] = mod(df.shape[end],2pi)
                #     # stop if mu == pi or close to pi
                #     # if abs(df.shape[end] - pi) < 0.02
                #     #     mus[i, j, tt+1:end, r] .= pi
                #     #     break
                #     # end
                # end
                df = params["q"] > 0 ? dft.defectsP[1] : dft.defectsN[1]
                mus[i,j,:,r] = df.shape
            catch e
                println(e)
                p=plot_thetas(thetas, model, lattice, defects=false,colorbar=false,quiver=false)
                display(p)
            end
        end
    end
end
prinz(z)
mus_avg = nanmean(mus,4)[:,:,:,1]
dmus_avg = diff(mus_avg, dims=3)/(times[2]-times[1])
# mus_avg_circular = zeros(length(sigmas), length(mu0s), length(times))
# for i in each(sigmas)
#     for j in each(mu0s)
#         for k in 1:length(times)
#             data = filter!(x->!isnan(x),mus[i,j,k,:])
#             mus_avg_circular[i,j,k] = mod(angle(sum(exp.(1im*data))),2π)
#         end
#     end
# end
prinz(z)
# @save "$(fn)one_defect_mu_dot_VC_many_mu0s.jld2" mus mus_avg sigmas mu0s times params runtime=z
# @save "$(fn)one_defect_mu_mu_dot_extension_sigma.jld2" mus mus_avg dmus_avg sigmas mu0s times params runtime = z
comments = "TriangularLattice, T=0"
# @save "$(fn)one_defect_mu_mu_dot_full_sigma_full_mu0s_Triangular.jld2" mus mus_avg dmus_avg sigmas mu0s times params comments runtime = z
@load "$(fn)one_defect_mu_mu_dot_full_sigma_full_mu0s_Triangular.jld2" mus mus_avg dmus_avg sigmas mu0s times params comments runtime
##
p_main = plot(legend=false)#, xlabel=L"µ", ylabel=L"\dot{µ}/\sigma")
plot!([pi, pi], [-0.5, 0.5], c=:black, ls=:dash, lw=0.5)
plot!([-0.5, 3pi], [0, 0], c=:black, ls=:dash, lw=0.5)
for i in 1:length(sigmas)
    ll = length(times) - 1
    for j in 1:1:Int(length(mu0s) / 2)
        plot!(mus_avg[i, j, 1:2:ll-1], dmus_avg[i, j, 1:2:ll] / sigmas[i], c=j, alpha=sqrt(i / length(sigmas)), label="μ = $(mu0s[j])", rib=0)
    end
    for j in reverse(Int(length(mu0s) / 2)+1:1:length(mu0s))
        plot!(mus_avg[i, j, 1:2:ll-1], dmus_avg[i, j, 1:2:ll] / sigmas[i], c=length(mu0s) + 1 - j, alpha=sqrt(i / length(sigmas)))
    end
end
xlims!((-0.1, 2pi + 0.1))
ylims!((-0.4, 0.42))
scatter!(mu0s[1:Int(end / 2)], dmus_avg[1, 1:Int(end / 2), 1] / sigmas[1], c=:black, ms=4, msw=0, label="", rib=0)
scatter!(mu0s[Int(end / 2)+1:end], dmus_avg[1, Int(end / 2)+1:end, 1] / sigmas[1], c=:black, ms=5, msw=0, label="", rib=0, m=:utriangle)
xticks!(0:pi/2:2pi, ["0", "π/2", "π", "3π/2", "2π"])
scatter!((pi, 0), c=:black, ms=12, msw=0, label="", rib=0, marker=:star5)
p_main

## --------------------
##
p1=plot(legend=:left,xlabel=L"t",ylabel=L"µ_{+}(t)")
for i in each(sigmas)
    plot!(times, mus_avg[i,j,:], label="$(sigmas[i])",rib=0,c=i)
    plot!(times, mus_avg[i,2,:],c=i)
end
hline!([pi],c=:black, ls=:solid,lw=0.5)
ylims!(0,2pi)
# xlims!(-0.5,301)
p1
# fpath = "D:/Documents/Research/projects/LatticeModels/NRI/NRXY/plotting_scripts_report/mu_dot/"
# savefig(p1,fpath*"mu_dot_VM.svg")
#
p2=plot(legend=false,xlabel=L"σt")
for i in 1:length(sigmas)
    # plot!(sinh(sigmas[i])*times, mus_avg[i,1,:], label="σ = $(sigmas[i])",c=i)
    plot!(sigmas[i]*times , mus_avg[i,1,:], label="σ = $(sigmas[i])",c=i)
end
hline!([pi],c=:black, ls=:solid,lw=0.5)
ylims!((0,2pi))
xlims!(-0.5,101)
cst = log(cot(mu0s[1]/2))
r = 3.5
plot!(p2,times*sigmas[end], x->2acot(exp(cst - 1/2/r*x)),ls=:dash,c=:black) # fits the data for 1NN 
p2
#
p3=plot(legend=false,xlabel=L"σt^{1+σ}")
for i in 2:length(sigmas)
    plot!((sigmas[i]*times).^(1+1sigmas[i]), mus_avg[i,1,:], label="σ = $(sigmas[i])",c=i)
    # plot!((sigmas[i])^(5/4)*times, mus_avg[i,1,:], label="σ = $(sigmas[i])",c=i)
end
hline!([pi],c=:black, ls=:solid,lw=0.5)
ylims!((0,2pi))
xlims!(-0.5,100)
cst = log(cot(mu0s[1]/2))
r = 1/sqrt(3)
# plot!(p3,times.^(1+sigmas[end])*sigmas[end], x->mu0s[1] + (pi-mu0s[1])*(1-exp(-1/14*x)),ls=:dash,c=:black)
p3

##
p4=plot(legend=false,xlabel=L"µ",ylabel=L"\dot{µ}/\sigma")
for i in 1:length(sigmas)
    for j in 1:length(mu0s)
        plot!(mus_avg[i,j,:][1:end-10-1],smooth(dmus_avg[i,j,:],over=1)[1:end-10]/sigmas[i], label="σ = $(sigmas[i])",c=i)
        # plot!(mus_avg[i,2,:][1:end-10-1],smooth(dmus_avg[i,2,:],over=1)[1:end-10], label="σ = $(sigmas[i])",c=i)
    end
end
xlims!((0,2pi))
# plot!(0:0.1:2pi, x->0.015*sin(x),c=:black,ls=:dash)
p4
##

# plot(p1,p2,p3,layout=(1,3),size=(1200,400)) # horizontal
plot(p1,p2,p4,layout=(3,1),size=(400,1200)) # vertical
# savefig("figures/issue_collapse.png")

##
# pplus=plot(legend=false,ylabel=L"µ_{+}(t)",xlabel=L"t")
# for i in each(mu0s)
#     plot!(times, mus_avg[1,i,:],c=i)
# end
# pplus
# ylims!(pplus, (0,2pi))
# pplus
##
pminus=plot(legend=false,ylabel=L"µ_{-}(t)",xlabel=L"t")
for i in 2:length(mu0s)-1
    plot!(times, smooth(mus_avg_circular[1,i,:],over=3),c=i)
end
pminus
ylims!(pminus, (0,2pi))
pminus

##
p=plot(pplus,pminus,layout=(1,2),size=(800,400))
fpath = "D:/Documents/Research/projects/LatticeModels/NRI/NRXY/plotting_scripts_report/mu_dot/"
# savefig(p,fpath*"mu_dot_many_µ0_VM.svg")



## -------------------- Inferring the f(r) function for the twisted defect for different times -------------------- 
include(srcdir("../parameters.jl"));
L = 50
params["L"] = L
params["T"] = 0.0
params["init"] = "single"
params["mu0"] = 1
params["q"] = 1
params["r0"] = Int(round(L/2))
params["sigma"] = 0.3

model = NRXY(params)
lattice = SquareLattice(L)
thetas = init_thetas(model, lattice, params)


# thetas[1,:] .= Float32(NaN)
# thetas[:,1] .= Float32(NaN)
# thetas[end,:] .= Float32(NaN)
# thetas[:,end] .= Float32(NaN)



times = 0:1:100
profiles = zeros(Float32, length(times), round(Int,L/2)+1)
mus = zeros(Float32, length(times))

z = @elapsed for i in each(times)
    evolve!(thetas, model, lattice, times[i])
    profiles[i,:] = thetas[round(Int,L/2),round(Int,L/2):end]# + thetas[round(Int,L/2)+1,round(Int,L/2):end]
    dft = DefectTracker(thetas, model, lattice)
    mus[i] = dft.defectsP[1].shape[1]
end
prinz(z) 

plot_thetas(thetas, model, lattice, defects=false,colorbar=false,quiver=false)
##
p=plot(axis=:log,legend=false,xlabel=L"r",ylabel=L"f(r)")#,ylims=(1E-3,1E-1))
rr = collect(2:round(Int,L/2)).+1
for i in 8:length(times)
    # plot!(rr, (remove_negative((profiles[i,2:end] - profiles[1,2:end])/mus[i]^0)),c=i)
    plot!(rr/times[i]^0, 1/times[i]*(remove_negative((profiles[i,2:end] - profiles[1,2:end])/mus[i]^0)),c=i)
end
# plot!(0.1:0.1:10, x->0.04x^(-1),c=:black)
p



## -------------------- Inferring the f(r) function for the twisted defect for different µ -------------------- 
include(srcdir("../parameters.jl"));
L = 50
params["L"] = L
params["T"] = 0.0
params["init"] = "single"
params["mu0"] = 0.3
params["q"] = 1
params["r0"] = Int(round(L/2))
params["sigma"] = 0.3

model = NRXY(params)
lattice = SquareLattice(L)
thetas = init_thetas(model, lattice, params)

thetas[1,:] .= Float32(NaN)
thetas[:,1] .= Float32(NaN)
thetas[end,:] .= Float32(NaN)
thetas[:,end] .= Float32(NaN)

dµ = 0.2
mus = collect(params["mu0"]:dµ:π-3dµ)
profiles = zeros(Float32, length(mus), round(Int,L/2)+1)

token_µ = 1
z = @elapsed while token_µ ≤ length(mus)
    update!(thetas, model, lattice)
    dft = DefectTracker(thetas, model, lattice)
    µ = dft.defectsP[1].shape[1]
    if µ ≥ mus[token_µ]
        half_arc = 0.5*arclength.(thetas[round(Int,L/2),round(Int,L/2):end],thetas[round(Int,L/2)+1,round(Int,L/2):end],2π)
        profiles[token_µ,:] = (thetas[round(Int,L/2),round(Int,L/2):end] .-params["mu0"] + half_arc)
        token_µ += 1
    end    
end
prinz(z) 

##
p=plot(axis=:log,legend=false,xlabel=L"r",ylabel=L"f(r)")#,ylims=(1E-3,1E-1))
firstind = 2 ; rr = collect(firstind:round(Int,L/2)+1)
for i in 1:length(mus)-2
    plot!(rr, 1/mus[i]*(remove_negative((profiles[i,firstind:end]))),c=i)
end
# plot!(rr, x->exp(-(x-2)/10),c=:black)
plot!(rr, x->1/x,c=:black)
p

## Figure for the PRL
# prinz(runtime)
# @load "$(fn)one_defect_mu_mu_dot_full_sigma_full_mu0s_Triangular.jld2" mus mus_avg dmus_avg sigmas mu0s times params runtime

p_main = plot(legend=false, background_color=:transparent)#, xlabel=L"µ", ylabel=L"\dot{µ}/\sigma")
plot!([pi, pi], [-0.5, 0], c=:black, ls=:dash, lw=0.5)
plot!([-0.5, 3pi], [0, 0], c=:black, ls=:dash, lw=0.5)
for i in length(sigmas)
    ll = length(times) - 1
    couleur = 2
    for j in 2:2:Int(length(mu0s) / 2)
        plot!(mus_avg[i, j, 1:2:ll-1], dmus_avg[i, j, 1:2:ll] / sigmas[i], c=couleur, alpha=sqrt(i / length(sigmas)), label="μ = $(mu0s[j])", rib=0)
        couleur += 1
    end
    couleur = 2
    for j in reverse(Int(length(mu0s) / 2)+2:2:length(mu0s))
        plot!(mus_avg[i, j, 1:2:ll-1], dmus_avg[i, j, 1:2:ll] / sigmas[i], c=couleur, alpha=sqrt(i / length(sigmas)))
        couleur += 1
    end
end
xlims!((-0.1, 2pi + 0.1))
ylims!((-0.4, 0.42))
scatter!(mu0s[2:2:Int(end / 2)], dmus_avg[1, 2:2:Int(end / 2), 1] / sigmas[1], c=:black, ms=4, msw=0, label="", rib=0)
# scatter!(mu0s[Int(end / 2)+1:end], dmus_avg[1, Int(end / 2)+1:end, 1] / sigmas[1], c=:black, ms=5, msw=0, label="", rib=0, m=:utriangle)
scatter!(mu0s[Int(end / 2)+2:2:end], dmus_avg[1, Int(end / 2)+2:2:end, 1] / sigmas[1], c=:black, ms=5, msw=0, label="", rib=0, m=:utriangle)
xticks!(0:pi/2:2pi, ["0", "π/2", "π", "3π/2", "2π"])
scatter!((pi, 0), c=:black, ms=12, msw=0, label="", rib=0, marker=:star5)
p_main
annotate!((0.95, 0.06), text(L"µ\!_+", :black, 15))
annotate!((0.12, 0.945), text(L"\dot{µ}\!_+\!/\sigma", :black, 15))
plot!([1.05, 1.75], [0.03, 0.32], arrow=true, color=:black, linewidth=1)

annotate!((0.335, 0.93), text(L"\mu_0", :black, 12))
annotate!((0.665, 0.04), text(L"\mu_0", :black, 12))
plot!([2pi - 1.85, 2pi - 1.1], [-0.35, -0.04], arrow=true, color=:black, linewidth=1)
scatter!([(mu0s[2], dmus_avg[1, 2, 1] / sigmas[1])], c=2, ms=8, msw=0, label="", rib=0, m=:circle)
scatter!([(mu0s[2], dmus_avg[1, 2, 1] / sigmas[1])], c=:black, ms=3, msw=0, label="", rib=0, m=:xcross)
scatter!([(mu0s[end-1], dmus_avg[1, end-1, 1] / sigmas[1])], c=2, ms=9, msw=0, label="", rib=0, m=:utriangle)
scatter!([(mu0s[end-1], -0.006 + dmus_avg[1, end-1, 1] / sigmas[1])], c=:black, ms=3, msw=0, label="", rib=0, m=:xcross)
# savefig(p_main, "$(fn)mu_dot_mu.svg")

## Now the inset representing various sigmas
p_inset = plot(legend=false, size=(300, 250), xlabel=L"µ\!_+", ylabel=L"\dot{µ}\!_+/\sigma")
plot!([-0.5, 3pi], [0, 0], c=:grey, ls=:dash, lw=0.5)
plot!([pi, pi], [-0.5, 0.5], c=:grey, ls=:dash, lw=0.5)
ind = 2
for i in 1:length(sigmas)
    ll = length(times) - 1
    for j in ind
        plot!(mus_avg[i, j, 1:2:ll-1], dmus_avg[i, j, 1:2:ll] / sigmas[i], c=ind, label="μ = $(mu0s[j])", rib=0, lw=1, alpha=max(0.4, i / length(sigmas)))
    end
    for j in length(mu0s) + 1 - ind
        plot!(mus_avg[i, j, 1:2:ll-1], dmus_avg[i, j, 1:2:ll] / sigmas[i], c=ind, lw=1, alpha=max(0.4, i / length(sigmas)))
    end
end
xlims!((-0.1, 2pi + 0.1))
ymax = 1.1 * maximum(dmus_avg[:, ind, 1:2:length(times)-1-1] ./ sigmas)
ylims!((-ymax, ymax))
# add an arrow from 0.1 , 0.2 to 0.3 , 0.4
plot!([0.6, 3], [0.015, 0.1], arrow=true, color=:black, linewidth=1)
annotate!((0.55, 0.9), text(L"\sigma", :black, 12))
annotate!((0.45, 0.1), text(L"\sigma", :black, 12))
plot!([2pi - 0.6, 2pi - 3], [-0.015, -0.1], arrow=true, color=:black, linewidth=1)
xticks!(0:pi:2pi, ["0", "π", "2π"])
yticks!(-0.1:0.05:0.1, ["-0.1", "", "0", "", "0.1"])
p_inset
# savefig(p_inset, "$(fn)mu_dot_mu_inset.svg")

