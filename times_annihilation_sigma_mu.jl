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

## ------------ Annihilation times ------------ ##
## ------------ Annihilation times ------------ ##
## ------------ Annihilation times ------------ ##
## ------------ Annihilation times ------------ ##


filename = "data/times_annihilation_pair.jld2"
@load filename sigmas Ts mus_plus r0s times_annihilation_all times_annihilation_avg times tmax params comments runtimes R
hrun(runtimes)
hrun(runtimes/3*10*2)
mus_plus
sigmas
Ts
params

indT = 1
indRo = 3
times_annihilation_avg = nanmean(times_annihilation_all, 5)
mean_XY_annihilation_time = mean(times_annihilation_avg[1, :, :, indRo, 1])


# p1 = plot(xlabel=L"Initial $\mu_+$", ylabel=L"$\sigma$", size=(300, 300))
# data = log10.(times_annihilation_avg[:, :, indT, indRo, 1] / mean_XY_annihilation_time)
# heatmap!(p1, mus_plus, sigmas, data, c=:rainbow)
# xticks!(0:pi/2:pi, [L"$0$", L"$\pi/2$", L"$\pi$"], colorbar=true)
# title!(L"\log_{10}(\tau/\tau_{XY})")



p2 = plot(xlabel=L"Initial $\mu_+$", ylabel=L"$\sigma$", size=(300, 300))
data = log.(times_annihilation_avg[:, :, indT, indRo, 1] / mean_XY_annihilation_time)
heatmap!(p2, mus_plus, sigmas, data, c=:rainbow)
xticks!(0:pi/2:pi, [L"$0$", L"$\pi/2$", L"$\pi$"], colorbar=true)
title!(L"\log_{10}(\tau/\tau_{XY})")
# for i in 1:length(sigmas)
#     for j in 1:length(mus_plus)
#         if abs(data[i, j]) < 0.08
#             scatter!(p2, [mus_plus[j]], [sigmas[i]], color=:white, m=:circle, ms=1.5)
#         end
#     end
# end
xlims!(0, pi)
ylims!(0, 0.35)
p2
# savefig(p2, fn * "figure3_time_annihilation_phase_space.svg")


