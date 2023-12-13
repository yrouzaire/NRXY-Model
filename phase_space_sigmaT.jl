cd("D:/Documents/Research/projects/LatticeModels") 
using DrWatson; 
@quickactivate "LatticeModels";
include(srcdir("LatticeModels.jl"));
using Plots, ColorSchemes, LaTeXStrings 
pyplot(box=true, fontfamily="sans-serif", label=nothing, palette=ColorSchemes.tab10.colors[1:10], grid=false, markerstrokewidth=0, linewidth=1.3, size=(400, 400), thickness_scaling=1.5); plot(); 
include(srcdir("../parameters.jl"));
cols = cgrad([:black, :blue, :green, :orange, :red, :black]);
plot();
fn = "NRXY/investigation/system/"



## ---------------- Analysis Data Cluster ---------------- ##
## Data Processing 
filename = "data/NRXY_long1E5_big400_sigma.jld2"
# filename = "data/NRXY_square_discretisation_impact.jld2"
@load filename Ps Cs xis ns Es sigmas mu_poss mu_negs Ts times params comments runtimes
hrun(runtimes) 

times
mu_poss
params["L"]
params["init"]
sigmas
Ts
Ps
L = params["L"]
Ps_avg = nanmean(Ps,4)[:,:,:,1]
xis_avg = nanmean(xis,4)[:,:,:,1]
for i in each(xis_avg)
	if iszero(xis_avg[i]) 
		xis_avg[i] = NaN
	end
end
ns_avg = nanmean(ns,4)[:,:,:,1]
Es_avg = nanmean(Es,4)[:,:,:,1]

good_indices = []
R = 40
for r in 1:R 
	try 
		Cs[1,1,1,r]
		push!(good_indices,r)
		catch ;
	end
end

Cs_avg = Array{Vector}(undef, length(sigmas),length(Ts), length(times))
for i in 1:length(sigmas)
	for j in 1:length(Ts)
		for k in 1:length(times)
			Cs_avg[i,j,k] = mean([Cs[i,j,k,r] for r in good_indices])
		end
	end
end
## ----------- Heatmaps of Phase Space at final time ----------- ##
p1 = plot(xlabel="σ",ylabel="T",size=(480,400),colorbartitle="P")
heatmap!(sigmas,Ts,Ps_avg[:,:,end]',c=cgrad([:red,:orange,:green,:blue]),clims=(0,1))

p2 = plot(xlabel="σ",ylabel="T",size=(480,400),colorbartitle=L"\log_{10}(n+1)")
heatmap!(sigmas, Ts,log10.(ns_avg[:,:,end]' .+ 1),c=cgrad(reverse([:red,:orange,:green,:blue])),clims=(0,log10(maximum(ns_avg[:,:,end] .+ 1))))

p3 = plot(xlabel="σ",ylabel="T",size=(480,400),colorbartitle=L"ξ/L")
heatmap!(sigmas, Ts,xis_avg[:,:,end]'/L,c=cgrad(([:red,:orange,:green,:blue])),clims=(0,.5))

p4 = plot(xlabel="σ",ylabel="T",size=(480,400),colorbartitle=L"E - E_0")
heatmap!(sigmas, Ts,Es_avg[:,:,end]'.+ 1,c=cgrad(reverse([:red,:orange,:green,:blue])),clims=(0,.65))

plot(p1,p2,p3,p4,layout=(2,2),size=(960,800))
title!("t = $(Int(times[end]))")
# savefig("NRI/NRXY/investigation/system/phase_space_sigmaT_tmax$(Int(times[end]))_L$(params["L"]).png")

## ----------- Movies Heatmaps of Phase Space over time ----------- ##
anim = @animate for tt in each(times)
	p1 = plot(xlabel="σ",ylabel="T",size=(480,400),colorbartitle="P")
	heatmap!(sigmas, Ts,Ps_avg[:,:,tt]',c=cgrad([:red,:orange,:green,:blue]),clims=(0,1))

	p2 = plot(xlabel="σ",ylabel="T",size=(480,400),colorbartitle=L"\log_{10}(n+1)")
	heatmap!(sigmas, Ts,log10.(ns_avg[:,:,tt]' .+ 1),c=cgrad(reverse([:red,:orange,:green,:blue])),clims=(0,log10(maximum(ns_avg[:,:,end] .+ 1))))

	p3 = plot(xlabel="σ",ylabel="T",size=(480,400),colorbartitle=L"ξ/L")
	heatmap!(sigmas, Ts,xis_avg[:,:,tt]'/L,c=cgrad(([:red,:orange,:green,:blue])),clims=(0,.5))

	p4 = plot(xlabel="σ",ylabel="T",size=(480,400),colorbartitle=L"E - E_0")
	heatmap!(sigmas, Ts,Es_avg[:,:,tt]'.+ 1,c=cgrad(reverse([:red,:orange,:green,:blue])),clims=(0,.5))

plot(p1,p2,p3,p4,layout=(2,2),size=(960,800))
	plot(p1,p2,p3,p4,layout=(2,2),size=(960,800))
	title!("t = $((times[tt]))")
end
filename = "NRI/NRXY/investigation/system/films/phase_space_sigmaT_L$(params["L"])_A.mp4"
mp4(anim,filename,fps=10)

## ----------- Measurements over time ----------- ##
ind_T = 1

p1 = plot(xlabel="t",ylabel="P",axis=:log,legend=:topleft)
for i in 1:length(sigmas)
	for j in ind_T
		plot!(times,Ps_avg[i,j,:],c=i, rib=0, label="σ = $(sigmas[i])")
		# plot!((sigmas[i]*(times).^(1-sigmas[i])).^-0.5,Ps_avg[i,j,:],c=i)
	end
end
# plot!(times[8:end],x->1E-2*sqrt(x/log(x)),c=:black)
p1

##
p2 = plot(axis=:log, size=(350,350), legend=:bottomleft, legend_title="σ")
# p2 = plot(xlabel=L"σ\,t^{1-σ}",ylabel=L"n/L^2",axis=:log)
for i in 1:length(sigmas)
	for j in ind_T
		plot!(times,ns_avg[i,j,:]/L^2,c=i, label="$(sigmas[i])", rib=0)
		# plot!((sigmas[i]*times).^((1-sigmas[i])),ns_avg[i,j,:]/L^2,c=i)
		# plot!(sigmas[i]*(times).^((1-sigmas[i])),ns_avg[i,j,:]/L^2,c=i)
	end
end

plot!(times[15:end], x -> 5E-2log(10x) / x, c=:black, label=L"(\log t) / {t}")
plot!(times[15:end], x -> 7E-2 * x^(-0.45), c=:black, line=:dot, label=L"t^{-0.45}")
xticks!(10.0 .^ (0:5), [L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}", L"10^{5}"])
yticks!(10.0 .^ (-5:-1), [L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}", L"10^{-1}"])
xlims!(0.8,5E5)
annotate!((0.25,0.92),text(L"n/L^2", 15, :left, :bottom, :black, :hcenter, :vcenter))
annotate!((0.95,0.06),text(L"t", 15, :left, :bottom, :black, :hcenter, :vcenter))
p2
savefig(p2,"NRXY/figures_paper/nt.svg")

##

p3 = plot(xlabel="t",ylabel=L"ξ/L",axis=:log,legend=:topleft)
for i in 2:length(sigmas)
	for j in ind_T
		plot!(times,xis_avg[i,j,:]/L,c=i)
	end
end
# plot!(times[15:end],x->50E-3*sqrt(x/log(x)),c=:black,label=L"\sqrt{\log(x)/x}")
p3



p4 = plot(xlabel="t",ylabel=L"\sqrt{n}ξ/L",xaxis=:log)
for i in 2:length(sigmas)
	for j in ind_T
		plot!(sigmas[i]*(times).^((1-sigmas[i])),xis_avg[i,j,:]/params["L"].*(ns_avg[i,j,:]).^0.5,c=i)
	end
end
plot!(times[8:end], x -> 5E-3 * sqrt(x / log(x)), c=:black)
plot!(times[8:end], x -> 5E-3 * x^(-0.7), c=:black, line=:dash)
## ICI ICI ICI 
p4
plot(p1,p2,p3,layout=(1,3),size=(1200,400))
# savefig(fn*"figures/measurements_sigmaT.pdf")
##

p2_collapse = plot(xlabel=L"σ\,t^{1-σ}",ylabel=L"n/L^2",axis=:log)
for i in 2:length(sigmas)
	for j in ind_T
		# plot!(times,ns_avg[i,j,:]/L^2,c=i)
		plot!(sigmas[i]*(times).^((1-sigmas[i])),ns_avg[i,j,:]/L^2,c=i)
	end
end
plot!(times[2:end-10],x->50E-2 /x,c=:black, label=L"\sim 1/x")
# ylims!(1E0,1e3)
yticks!([1e-3,1e-2,1e-1,1e0],[L"10^{-3}",L"10^{-2}",L"10^{-1}",L"10^{0}"])
xticks!(10.0 .^(-1:1:2),[L"10^{-1}",L"10^{0}",L"10^{1}",L"10^{2}"])
p2_collapse

p3_collapse = plot(xlabel=L"σ\,t^{1-σ}",ylabel=L"ξ/L",axis=:log,legend=:topleft)
for i in 2:length(sigmas)
	for j in ind_T
		# plot!(times,xis_avg[i,j,:]/L,c=i)
		plot!((sigmas[i]*times).^(1-sigmas[i]),xis_avg[i,j,:]/L,c=i)
	end
end
# plot!(times[8:end],x->5E-3*sqrt(x/log(x)),c=:black)
plot!(times[1:20],x->70E-3*sqrt(x),c=:black,label=L"\sim x")
yticks!([1e-1,1e0],[L"10^{-1}",L"10^{0}"])
xticks!(10.0 .^(-1:1:2),[L"10^{-1}",L"10^{0}",L"10^{1}",L"10^{2}"])
p3_collapse

plot(p2_collapse,p3_collapse,layout=(1,2),size=(800,400))
# savefig(fn*"figures/measurements_sigmaT_collapse.pdf")

## ----------- Correlation functions for different times ----------- ##
ind_T= 1
ind_sig = 1
couleurs = cgrad([:red,:orange,:green,:blue,:purple])
ind_points_in_phase_space = [(1,1),(1,2),(1,3),(1,5)] # T, sigma
times_to_plot = each(times)
times_to_plot = [1,3,5,7,9,11,13,15,17,19,21,23,24,25,26,27,28,29,30]

pp = Vector(undef,length(ind_points_in_phase_space))
rr = 1:Int(L/2)
for i in each(ind_points_in_phase_space)
	ind_T, ind_sig = ind_points_in_phase_space[i]
	p=plot(xlabel="r",ylabel="C(r)",axis=:log)
	for tt in times_to_plot
		plot!(rr,remove_negative(Cs_avg[ind_sig,ind_T,tt]),c=couleurs[(times[tt]/times[end])^0.5])
	end
	title!("T = $(Ts[ind_T]), σ = $(sigmas[ind_sig])")
	ylims!(1E-2,1.05)
	pp[i] = p
end

plot(pp...,layout=(2,2),size=(960,800))


## ----------- P(µ) over time ----------- ##
ind_T = 1
ind_sig = 1#length(sigmas)
dµ = 0.25 # precision on µ ~ ± 0.2
nbbins = 0:dµ:2pi 
tts = ([1,6,8,9,10,11,12,13,16])
# tts = reverse([30])
# tts = reverse([29])

p_muplus=plot(legend=:topright)
# p_muplus=plot(xlabel=L"µ_+",ylabel=L"P(µ_+)",xlims=(0,2pi),legend=:topright)
coll = length(tts)+1
for tt in (tts)
	coll += -1
	data = []
	for ind in good_indices
		push!(data,mod.(mu_poss[ind_sig,ind_T,tt,ind],2pi)...)
	end
	h = fit(Histogram, data, nbbins)
	plot!(nbbins .+ dµ/2,(h.weights/length(data)/dµ),c=coll,label="t = $(round(Int,times[tt]))",rib=0)
end
plot!([0,2pi],x->1/2pi,c=:black,ls=:dash,lw=0.7,label="1/2π")
vline!(p_muplus,[pi],c=:black,lw=0.3)
xlims!(0,7.5)
annotate!(p_muplus, (0.935,0.06), text(L"µ_{+}", 15, :left, :bottom, :black, :hcenter, :vcenter))
annotate!(p_muplus, (0.15,0.93), text("P"*L"(µ_{+})", 15, :left, :bottom, :black, :hcenter, :vcenter))
# ylims!(-0.1,1)
xticks!(0:pi/2:2pi,["0","π/2","π","3π/2","2π"])
p_muplus
# savefig(fn*"figures/Pµ+.svg")
# savefig(fn*"figures/Pµ+.pdf")
#
p_muminus=plot(xlims=(0,2pi),uaxis=:log)
coll = length(tts)+1
for tt in tts
	coll += -1
	data = []
	for ind in good_indices
		push!(data,mod.(mu_negs[ind_sig,ind_T,tt,ind],2pi)...)
	end
	h = fit(Histogram, data, nbbins)
	plot!(nbbins.+ dµ/2,h.weights/length(data)/dµ,c=coll)
end
xticks!(0:pi/2:2pi,["0","π/2","π","3π/2","2π"])
vline!(p_muminus,[i*pi/2 for i in [1,3]],c=:black,lw=0.5)
# xticks!(0:pi/3:2pi,["0","π/3","2π/3","π","4π/3","5π/3","2π"])
# vline!(p_muminus,[i*pi/3 for i in [1,2,4,5]],c=:black,lw=0.5)
plot!(p_muminus,[pi,pi],[0.03,0.31],c=:black,lw=0.5)
hline!(p_muminus,[1/2pi],c=:black,ls=:dash)
annotate!(p_muminus, (0.92,0.06), text(L"µ_{-}", 15, :left, :bottom, :black, :hcenter, :vcenter))
annotate!(p_muminus, (0.5,0.93), text("P"*L"(µ_{-})", 15, :left, :bottom, :black, :hcenter, :vcenter))
ylims!(0.05,0.36)
plot(p_muplus,p_muminus,layout=(1,2),size=(800,400))
# savefig(fn*"figures/Pµ.pdf")
# savefig(p_muminus,fn*"figures/Pµ-.svg")
# savefig(fn*"figures/Pµ-.pdf")

## Scaling defect density
ind_T = 1
p = plot(axis=:log,legend=:topright)
for i in 1:length(sigmas)
	for j in ind_T
		plot!(times,ns_avg[i,j,:]/L^2,c=i, label="σ = $(sigmas[i])",rib=0)
	end
end
# plot!(times[8:end-10],x->3E-2/x,c=:black)
ylims!(7E-5,0.6)
yticks!(10.0 .^(-4:1:0))
annotate!(p, (0.2,0.93),text(L"n/L^2", 15, :left, :bottom, :black, :hcenter, :vcenter))
annotate!(p, (0.9,0.044),text(L"t", 15, :left, :bottom, :black, :hcenter, :vcenter))
p
#

p_collapse = plot(size=(300,300),axis=:log)
for i in 2:length(sigmas)
	for j in ind_T
		plot!(sigmas[i]*(times).^((1-sigmas[i])),ns_avg[i,j,:]/L^2,c=i)
	end
end
ylims!(5E-5,0.6)
yticks!(10.0 .^(-4:1:0))
xticks!(10.0 .^(-1:1:2))
annotate!(p_collapse, (0.5,0.9),text(L"n/L^2", 15, :left, :bottom, :black, :hcenter, :vcenter))
annotate!(p_collapse, (0.75,0.06),text(L"σt^{1-σ}", 15, :left, :bottom, :black, :hcenter, :vcenter))
annotate!(p_collapse, (0.47,0.25),text(L"\sim \frac{\log x}{x}", 15, :left, :bottom, :black, :hcenter, :vcenter))
plot!(times[2:end-13],x->4E-3(log(100x)/x),c=:black)
# ylims!(1E0,1e3)
p_collapse

# savefig(p,fn*"figures/n_scaling.png")
# savefig(p_collapse,fn*"figures/n_scaling_collapse.svg")
