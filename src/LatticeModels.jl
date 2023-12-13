include("lattices.jl")
include("models.jl")
include("init_visu.jl")

include("core_methods.jl")
include("defects_methods.jl")

include("auxiliary.jl")
include("measurements.jl")

using BenchmarkTools, JLD2
global const WINDOW = 7
global const W21 = 2WINDOW+1

# using Flux
# using Flux:onecold, Chain, Dense, softmax
# struct Reshape # has to be defined before Augmentor.jl is loaded to avoid conflicts
#     shape
# end
# Reshape(args...) = Reshape(args)
# (r::Reshape)(x) = reshape(x, r.shape)
# Flux.@functor Reshape ()

# include("infer_mu_src.jl")
# include("div_rot.jl")

# using BSON
# BSON.@load "NeuralNets/DAE_positive1___15_02_2023.bson" DAE_positive1
# global const DAE_polar = DAE_positive1;

# xx = provide_div_rot(ones(Float32,2WINDOW+1,2WINDOW+1,1))
# DAE_polar(xx)
# println("aaa  ",DAE_polar(xx)[1])


# # # using BSON
# global const DAE_positive1  = load("NeuralNets/DAE_positive1___15_02_2023.jld2","DAE_positive1");
# global const DAE_positive1  = load("NeuralNets/DAE_positive1___10_01_2023.jld2","DAE_positive1");
# global const DAE_negative1  = load("NeuralNets/DAE_negative1___16_12_2022.bson","DAE_negative1")
# global const DAE_positive12 = load("NeuralNets/DAE_positive12___15_12_2022.bson","DAE")
# global const DAE_negative12 = load("NeuralNets/DAE_negative12___15_12_2022.bson","DAE")

"All good !"