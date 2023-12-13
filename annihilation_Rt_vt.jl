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
fn = "NRXY/figures_paper/"


## ----------------- Load data ----------------- ##
filename = "data/Rt_vt_annihilation_stabilisation_track_defects.jld2"
filename = "data/Rt_L200_annihilation_stabilisation.jld2"
@load filename Ts dfts_fusion runtimes models R inits_mu sigmas_VM times tmax params  
hrun(runtimes)

sigmas_VM
params["L"]
inits_mu
Reff = length(dfts_fusion)
Ts
times
sigmas_VM
# sigmas_VC = round.(sigmas_VC,digits=2)
models

# # ----------------- Analysis R(t) ----------------- ##
rts = zeros(length(sigmas_VM), length(inits_mu), length(Ts), length(models), length(times), R)

latt = TriangularLattice(params["L"])
for i in each(sigmas_VM), j in each(inits_mu), k in each(Ts), m in each(models), r in 1:Reff
	dft = dfts_fusion[r][i,j,k,m]
	rt = interdefect_distance(dft.defectsP[1],dft.defectsN[1],latt)
	rts[i, j, k, m, 1:length(rt), r] = rt
end
rts_avg = nanmean(rts,6)[:,:,:,:,:,1] 
rts_std = nanstd(rts,6)[:,:,:,:,:,1]/sqrt(Reff)

rts_plus = NaN*zeros(length(sigmas_VM), length(inits_mu), length(Ts), length(models), length(times), R)
rts_minus = NaN*zeros(length(sigmas_VM), length(inits_mu), length(Ts), length(models), length(times), R)
for i in each(sigmas_VM), j in each(inits_mu), k in each(Ts), m in each(models), r in 1:Reff
	dft = dfts_fusion[r][i,j,k,m]
    for tt in 1:length(dft.defectsP[1].pos)
        rts_plus[i, j, k, m, tt, r] = dist(latt, dft.defectsP[1].pos[tt], dft.defectsP[1].pos[1])
    end
    for tt in 1:length(dft.defectsN[1].pos)
        rts_minus[i, j, k, m, tt, r] = dist(latt, dft.defectsN[1].pos[tt], dft.defectsN[1].pos[1])
    end
end
rts_plus_avg = nanmean(rts_plus,6)[:,:,:,:,:,1]
rts_minus_avg = nanmean(rts_minus,6)[:,:,:,:,:,1]
# data = filter(!isnan, rts_plus_avg[1, 1, 1, 1, 2:end])
# plot(times[1:length(data)], data)

##----------------- Plotting R(t) ----------------- ##
##----------------- Plotting R(t) ----------------- ##
##----------------- Plotting R(t) ----------------- ##
##----------------- Plotting R(t) ----------------- ##

plot_every = 1
styles = [:dashdot, :solid, :dot]
p1 = plot(xaxis=:log, legend=:false, dpi=1000, size=(300, 300))#, background_color=:transparent)#,title="Von Mises")
for i in 2:length(sigmas_VM)
    for j in each(inits_mu)
        for k in 1 #each(Ts)
            for m in 1 #each(models)
                j == 1 ? lwdt = 1.3 : lwdt = 0.9
                # plot!(times[2:length(rts_avg[i, j, k, m, :])], remove_negative(rts_avg[i, j, k, m, 2:end] ./ rts_avg[i, j, k, m, 1]),
                    # rib=0rts_std[i, j, k, m, :], c=i, ls=styles[j], lw=lwdt)
                plot!(times[2:plot_every:length(rts_avg[i, j, k, m, :])], remove_negative(rts_avg[i, j, k, m, 2:plot_every:end]),
                    rib=0rts_std[i, j, k, m, :], c=i, ls=styles[j], lw=lwdt)
            end
        end
    end
end
plot!(times[2:plot_every:length(rts_avg[1, 1, 1, 1, :])], remove_negative(rts_avg[1, 1, 1, 1, 2:plot_every:end]), rib=0, c=1, label="σ=0", ls=:solid, lw=2)
for i in 2:length(sigmas_VM)
    plot!([NaN, NaN], rib=0, c=i, label="$(round(sigmas_VM[i], digits=2))")
end
# plot!(times[1:length(rts_avg[1, 1, 1,1, :])], rts_avg[1, 1, 1,1, 1:end], rib=0, c=1, ls=:solid, lw=1.3)
# yticks!(10.0 .^ (-1:2), [L"10^{-1}", L"10^{0}", L"10^1", L"10^2"])
# ylims!(.01, 200)
ylims!(0.000, 115)
xlims!(10, 20E5)
xticks!(10.0 .^ (1:1:6), [L"10^1", L"10^2", L"10^3", L"10^4", L"10^5", L"10^6"])
# yticks!(10.0 .^ (-3:1:0), [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])

# xticks!(0:1000:6000,["0", "1", "2", "3", "4", "5", "6"])
annotate!((0.03, 0.935), text(L"R(t)", :left))
annotate!((0.94, 0.05), text(L"t", :left))
# annotate!((0.93, 0.03), text(L"t [\times 10^3]", :left))
# annotate!((0.03, 0.95), text(L"R_0\,(1-\alpha t^2)", :left))
p1

# savefig(fn * "figure3_Rt.png")

##
p1_collapse = plot(uaxis=:log, legend=:bottomleft, size=(400, 400))#,xlabel=L"t/τ", ylabel=L"R/R_0")#, background_color=:transparent)#,title="Von Mises")
for i in 1:length(sigmas_VM)
    for j in 1#each(inits_mu)
        for k in 1 #each(Ts)
            for m in 1 #each(models)
                j == 1 ? lwdt = 1.3 : lwdt = 0.9
                index_time_end = findfirst(iszero, rts_avg[i, j, k, 1, :])
                sig = sigmas_VM[i]
                # tt = times[2:index_time_end] .* (1 + (60sigmas_VM[i])^1.25)
                tt = times[2:index_time_end] ./ times[index_time_end]
                data = remove_negative((rts_avg[i, j, k, m, 2:end]) ./ rts_avg[i, j, k, m, 1])
                plot!(tt, data,
                    rib=0rts_std[i, j, k, m, :], c=i, ls=styles[2], lw=lwdt)

                # plot!(times[2:length(rts_avg[i, j, k, m, :])], rts_avg[i, j, k, m, 2:end],
                #     rib=rts_std[i, j, k, m, :], c=i)
            end
        end
    end
end
p1_collapse
# ylims!(0.000, 1.15)
annotate!((0.03, 0.925), text(L"R(t)/R_0", :left))
# xlims!(0,1.15)
annotate!((0.851, 0.06), text(L"t/τ", :left))
# savefig(fn * "figure3_Rt_collapse.svg")
# fpath = "D:/Documents/Research/projects/LatticeModels/NRI/NRXY/plotting_scripts_report/"
# savefig(p1,fn*"annihilation_Rt.png")

##----------------- Plotting R(t) for plus and minus differentiated ----------------- ##
##----------------- Plotting R(t) for plus and minus differentiated ----------------- ##
##----------------- Plotting R(t) for plus and minus differentiated ----------------- ##
##----------------- Plotting R(t) for plus and minus differentiated ----------------- ##


pdiff_pm = plot(uaxis=:log, legend=false, size=(300, 300))#, background_color=:transparent)#,title="Von Mises")
for i in length(sigmas_VM):-1:1
    for j in 1#each(inits_mu)
        for k in 1 #each(Ts)
            for m in 1 #each(models)
                j == 1 ? lwdt = 1.3 : lwdt = 0.9
                data = 2*filter(!isnan, rts_plus_avg[i, j, k, m, 2:end])/ params["L"]
                tA = times[findfirst(iszero, rts_avg[i, j, k, m, :])]
                tt = (times[1:length(data)] / tA)# .^(3/2)       
                plot!(tt, data, rib=0, c=i, ls=:solid, lw=lwdt)
                # plot!(times[1:length(data)], data, rib=0, c=i, ls=:solid, lw=lwdt)

                data = 2*filter(!isnan, rts_minus_avg[i, j, k, m, 2:end])/ params["L"]
                tA = times[findfirst(iszero, rts_avg[i, j, k, m, :])]
                tt = (times[1:length(data)] / tA) #.^(3/2)  
                plot!(tt, data, rib=0, c=i, ls=:dot, lw=lwdt)
                # plot!(times[1:length(data)], data, rib=0, c=i, ls=:dash, lw=lwdt)
                
            end
        end
    end
end
# for i in 2:length(sigmas_VM)
#     plot!([NaN, NaN], rib=0, c=i, label="$(round(sigmas_VM[i], digits=2))")
# end
# plot!(times[2:length(rts_avg[1, 1, 1,1, :])], rts_avg[1, 1, 1,1, 2:end], rib=0, c=1, label="σ = $(sigmas_VM[1])", ls=:solid, lw=1.3)
# yticks!(10.0 .^ (-1:2), [L"10^{-1}", L"10^{0}", L"10^1", L"10^2"])
# ylims!(.01, 200)
# ylims!(0.000, 1.1)
# xlims!(10, 12E4)
# xticks!(10.0 .^ (1:1:5), [L"10^1", L"10^2", L"10^3", L"10^4", L"10^5"])
# yticks!(10.0 .^ (-3:1:0), [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])

# xticks!(0:1000:6000,["0", "1", "2", "3", "4", "5", "6"])
annotate!((0.03, 0.92), text(L"d_\pm /R_0", :left))
annotate!((0.81, 0.07), text(L"t/\tau", :left, 15))
annotate!((0.65, 0.23), text("     -1 \n defects",10, :left))
annotate!((0.25, 0.6), text("     +1 \n defects", 10, :left))
# hline!([0.5], c=:grey, ls=:dash, lw=0.8)
# scatter!([0, 0], rib=0, line=:solid, c=:grey, label="a", lw = 2)
# scatter!([0, 0], rib=0, line=:dot, c=:grey, label="a", lw = 2)

# annotate!((0.03, 0.95), text(L"R_0\,(1-\alpha t^2)", :left))
pdiff_pm

# savefig(fn * "figure3_Rt_plus_minus.svg")

## ----------------- Analysis v(t) and a(t) ----------------- ##
using LinearAlgebra

overs = [1000, 400, 300, 100, 50, 25, 20, 15]
velP = NaN*zeros(length(sigmas_VM),length(inits_mu), length(Ts),length(models),length(times),R)
velN = NaN*zeros(length(sigmas_VM),length(inits_mu), length(Ts),length(models),length(times),R)

for i in each(sigmas_VM), j in each(inits_mu), k in each(Ts), m in each(models), r in 1:Reff
	dft = dfts_fusion[r][i,j,k,m]
	over = overs[i]
    if length(dft.defectsP[1].pos) > over && length(dft.defectsN[1].pos) > over
		vtP = compute_velocity_defect(dft.defectsP[1],times,over=over)
		vtN = compute_velocity_defect(dft.defectsN[1],times,over=over)
		ll = min(length(vtP),length(vtN))
		velP[i,j,k,m,1:ll,r] = vtP[1:ll]
		velN[i,j,k,m,1:ll,r] = vtN[1:ll]
	end
end
velP_avg = nanmean(velP,6)[:,:,:,:,:,1]
velN_avg = nanmean(velN,6)[:,:,:,:,:,1]

# over = 3
# accP = NaN*zeros(length(sigmas_VM),length(inits_mu), length(Ts),length(models),length(times),R)
# accN = NaN*zeros(length(sigmas_VM),length(inits_mu), length(Ts),length(models),length(times),R)

# for i in each(sigmas_VM), j in each(inits_mu), k in each(Ts), m in each(models), r in 1:Reff
# 	dft = dfts_fusion[r][i,j,k,m]
# 	atP = compute_acceleration_defect(dft.defectsP[1],times,over=over)
# 	atN = compute_acceleration_defect(dft.defectsN[1],times,over=over)
# 	ll = min(length(atP),length(atN))
# 	accP[i,j,k,m,1:ll,r] = atP[1:ll]
# 	accN[i,j,k,m,1:ll,r] = atN[1:ll]
# end
# accP_avg = nanmean(accP,6)[:,:,:,:,:,1]
# accN_avg = nanmean(accN,6)[:,:,:,:,:,1]

##
# ----------------- Plotting v(t) ----------------- ##
# ----------------- Plotting v(t) ----------------- ##
# ----------------- Plotting v(t) ----------------- ##
# ----------------- Plotting v(t) ----------------- ##
ind_to_cuts = [0, 200, 60, 20, 30, 10, 7, 7] # 
ind_to_cuts = [0,0,0,0,0,0,0,0] # 

p2=plot(legend=:topright, axis=:log)#,title="Von Mises")
for i in [1,2,3,4,5,6,7,8]#1:length(sigmas_VM)
	for j in 1
		for k in 1
			ll = findfirst(isnan.(velP_avg[i, j, k, 1, :]))
            plot!(times[2:ll-ind_to_cuts[i]], remove_negative(velP_avg[i, j, k, 1, 2:ll-ind_to_cuts[i]]), line=:solid, c=i, rib=0, m=false)
            plot!(times[2:ll-ind_to_cuts[i]], remove_negative(velN_avg[i, j, k, 1, 2:ll-ind_to_cuts[i]]), line=:dash, c=i)
        end
	end
end
xlims!(8, 5E4)
xticks!(10.0 .^ (1:1:5), [L"10^1", L"10^2", L"10^3", L"10^4", L"10^5"])
# ylims!(1E-3, 0.15)
annotate!((0.03, 0.93), text(L"v\!_\pm\!(t)", :left))
annotate!((0.94, 0.05), text(L"t", :left))
plot!([NaN,NaN], line=:solid,label=L"+1",c=:grey)
plot!([NaN,NaN], line=:dash, label=L"- 1",c=:grey)
p2
##

p2_collapse = plot(legend=:topright, uaxis=:log)#,title="Von Mises")
for i in [5, 6, 7, 8]#1:length(sigmas_VM)
    for j in 1
        for k in 1
            ll = findfirst(isnan.(velP_avg[i, j, k, 1, :]))
            index_time_end = findfirst(isnan, velP_avg[i, j, k, 1, :])
            sig = sigmas_VM[i]
            tA = 1/(1 + (60sigmas_VM[i])^1.25)
            tt = times[2:index_time_end] ./ tA
            plot!(tt.^(3/8), tA*remove_negative(velP_avg[i, j, k, 1, 2:index_time_end]), line=:solid, c=i, rib=0, m=false)
			plot!(tt.^(3/8), tA*remove_negative(velN_avg[i, j, k, 1, 2:index_time_end]), line=:dash, c=i)
            # tt = times[2:index_time_end] ./ times[index_time_end]
            # plot!(tt .^(3/10), remove_negative(velP_avg[i, j, k, 1, 2:index_time_end]) / sig^(5 / 4), line=:solid, c=i, rib=0, m=false)
            # plot!(tt .^(3/10), remove_negative(velN_avg[i, j, k, 1, 2:index_time_end]) / sig^(5 / 4), line=:dash, c=i)

        end
    end
end
# xlims!(0, 1)
# ylims!(0, 1.3)
annotate!((0.03, 0.93), text(L"v\!_\pm\!/σ^{5/4}", :left))
annotate!((0.54, 0.05), text(L"t/t_{annihilation}", :left))
p2_collapse

## Trend of v with sigma and t
# plot(sigmas_VM[2:end], velP_avg[2:end, 1, 1, 1, 100], marker=:circle, uaxis=:log)
p3=plot()
for ind_sig in 2:length(sigmas_VM)
	index_time_end = findfirst(isnan, velP_avg[ind_sig, 1,1,1, :])
	tt = times[2:index_time_end] ./ times[index_time_end]
	plot!(tt, velP_avg[ind_sig, 1, 1, 1, 2:end] / sigmas_VM[ind_sig]^3 / times[index_time_end], axis=:log)
end
p3

