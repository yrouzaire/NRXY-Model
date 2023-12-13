include("lattices.jl");
include("models.jl");
using StatsBase, Distributions, SpecialFunctions

## ------------------------ Get Neighbours ------------------------ ##
## ------------------------ Get Neighbours ------------------------ ##
## ------------------------ Get Neighbours ------------------------ ##
## ------------------------ Get Neighbours ------------------------ ##
## ------------------------ Get Neighbours ------------------------ ##
function get_neighbours(thetas::Matrix{<:T}, model::AbstractModel{T}, lattice::SquareLattice, i::Int, j::Int, bulk::Bool=false)::Vector{T} where {T<:AbstractFloat}
    L = lattice.L
    # convention depuis la droite et sens trigo
    if bulk
        jm = j - 1
        jp = j + 1
        imm = i - 1
        ip = i + 1
    else
        jm = mod1(j - 1, L)
        jp = mod1(j + 1, L)
        imm = mod1(i - 1, L)
        ip = mod1(i + 1, L)
    end

    @inbounds angles = [
        thetas[i, jp],
        thetas[imm, j],
        thetas[i, jm],
        thetas[ip, j]]

    # never used so far, slows down the code so until useful, I comment it out
    # if lattice.metric in ["manhattan","euclidian"]
    #     #= Ordering of the nearest neighbours (pure convention, motivated in the file lattices.jl)
    #        2
    #     3  X  1
    #        4
    #     =#
    #     @inbounds angles = [
            # thetas[i, jp],
            # thetas[imm, j],
            # thetas[i, jm],
            # thetas[ip, j]]

    # elseif lattice.metric == "chebychev"
    #     #= Ordering of the nearest neighbours (pure convention, motivated in the file lattices.jl)
    #     4  3  2
    #     5  X  1
    #     6  7  8
    #     =#
    #     @inbounds angles =
    #        [thetas[i,jp],
    #         thetas[imm,jp],
    #         thetas[imm,j],
    #         thetas[imm,jm],
    #         thetas[i,jm],
    #         thetas[ip,jm],
    #         thetas[ip,j],
    #         thetas[ip,jp]]
    # end
    if model.rho < 1
        return filter!(!isnan, angles)
    else
        return angles
    end
end

function get_neighbours(thetas::Matrix{<:T}, model::AbstractModel{T}, lattice::TriangularLattice, i::Int, j::Int, bulk::Bool=false)::Vector{T} where {T<:AbstractFloat}
    L = lattice.L

    # convention depuis la droite et sens trigo
    if bulk
        jm = j - 1
        jp = j + 1
        imm = i - 1
        ip = i + 1
    else
        jm = mod1(j - 1, L)
        jp = mod1(j + 1, L)
        imm = mod1(i - 1, L)
        ip = mod1(i + 1, L)
    end

    if iseven(i) # E stands for even
        #= Ordering of the nearest neighbours (pure convention, motivated in the file lattices.jl)
           3  2 imagine this line is shifted 1/2 lattice spacing to the left to give :    3  2
        4  E  1                                                                           4  E  1
           5  6 imagine this line is shifted 1/2 lattice spacing to the left to give :    5  6
        =#
        @inbounds angles =
            [thetas[i,jp],
            thetas[imm,jp],
            thetas[imm,j],
            thetas[i,jm],
            thetas[ip,j],
            thetas[ip,jp]]
    else # O stands for odd
        #= Ordering of the nearest neighbours (pure convention, motivated in the file lattices.jl)
        3  2                                                                               3  2
        4  O  1 imagine this line is shifted 1/2 lattice spacing to the left to give :  4  O  1
        5  6                                                                               5  6
        =#
        @inbounds angles =
            [thetas[i,jp],
            thetas[imm,j],
            thetas[imm,jm],
            thetas[i,jm],
            thetas[ip,jm],
            thetas[ip,j]]

    end
    if model.rho < 1
        return filter!(!isnan, angles)
    else
        return angles
    end
end

function get_neighbours(thetas::Vector{T}, model::AbstractPropagationModel{T}, lattice::Chain1D, i::Int, bulk::Bool=false)::Vector{T} where {T<:AbstractFloat}
    L = lattice.L
    if bulk
        return [thetas[i-1], thetas[i+1]]
    else
        if lattice.periodic
            return [thetas[mod1(i - 1, L)], thetas[mod1(i + 1, L)]]
        else
            if i == 1
                return [thetas[2]]
            elseif i == L
                return [thetas[L-1]]
            else
                return [thetas[i-1], thetas[i+1]]
            end
        end
    end
end

function get_neighbours(thetas::Matrix{<:T}, model::AbstractModel{T}, lattice::TorusLattice, i::Int, j::Int, bulk::Bool=false)::Vector{T} where {T<:AbstractFloat}
    L = lattice.L
    # convention depuis la droite et sens trigo
    if bulk
        jm  = j - 1
        jp  = j + 1
        imm = i - 1
        ip  = i + 1
    else
        jm  = mod1(j - 1, L)
        jp  = mod1(j + 1, L)
        imm = mod1(i - 1, L)
        ip  = mod1(i + 1, L)
    end

    @inbounds angles = [
        thetas[i, jp],  # right 
        thetas[imm, j], # up 
        thetas[i, jm],  # left 
        thetas[ip, j]]  # down

    return angles
end



# default method
function sum_influence_neighbours(theta::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::AbstractModel{T}, lattice::AbstractLattice)::T where {T<:AbstractFloat}
    force = 0
    for angle in angles_neighbours
        force += sin(angle - theta)
    end
    return force
    # this is ~40% faster than sum(sin, angles_neighbours .- theta) 2.77 ms versus 4.5 ms for a 200x200 system (one update)
    # and waaay faster than for i in each(angles_neighbours) ... # 32 ms
end

# for LangevinNematicXY
function sum_influence_neighbours(theta::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::LangevinNematicXY{T}, lattice::AbstractLattice)::T where {T<:AbstractFloat}
    force = 0
    for angle in angles_neighbours
        force += sin(2(angle - theta)) # 3.45 ms
    end
    return force
    # return sum(sin, 2(angles_neighbours .- theta)) # 6.2 ms
    # return sum(x -> sin(2x), angles_neighbours .- theta) # 4.9 ms
end

# for VisionConeXY
function sum_influence_neighbours(theta::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::VisionConeXY{T}, lattice::Abstract2DLattice)::T where {T<:AbstractFloat}
    weights = zeros(T, length(angles_neighbours))
    theta0 = mod(theta, 2π) # VisionConeXY is for now only defined for polar symmetry only
    if lattice isa TriangularLattice
        constant = T(π / 3)
    elseif lattice isa SquareLattice
        constant = T(π / 2)
    end
    force = 0.0

    # straightforward implementation, but slow
    # for n in 1:length(angles_neighbours)
    #     dtheta     = theta0 - (n-1)*constant
    #     dtheta_abs = abs(dtheta)
    #     arclengt   = min(2π-dtheta_abs,dtheta_abs)

    #     if arclengt ≤ model.vision_cone/2
    #         @inbounds weights[n] = 1.0
    #     end
    # end
    # if isempty(angles_neighbours) return 0.0 # instead of sum(sin,...,init=0) because not possible in Julia 1.3.0 on the cluster I use
    # else return sum(sin.(sym(model)*(angles_neighbours .- theta)) .* weights)
    # end

    # optimised version
    for n in eachindex(angles_neighbours)
        dtheta = theta0 - (n - 1) * constant
        dtheta_abs = abs(dtheta)
        arclengt = min(2π - dtheta_abs, dtheta_abs)

        if abs(arclengt) ≤ model.vision_cone / 2
            force += sin(angles_neighbours[n] - theta)
        end
    end
    # normalisation = model.vision_cone / 2π

    nnn = number_nearest_neighbours(lattice)
    return force / nnn
end

# for DiscreteNRXY
# pistes pour optimiser : separation des cas polar et nematic, puis for angle in angles_neighbours pour calculer sum(sin,angles_neighbours .- theta)
# function sum_influence_neighbours(theta::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::DiscreteNRXY{T}, lattice::Abstract2DLattice)::T where {T<:AbstractFloat}
#     nnn = number_nearest_neighbours(lattice)
#     if isempty(angles_neighbours)
#         return 0.0 # instead of sum(sin,...,init=0) because not possible in Julia 1.3.0 on the cluster I use
#     else
#         ID = ID_projection_angle_onto_lattice(theta, i, j, lattice)
#         if model.symmetry == "polar"
#             force = sum(sin, angles_neighbours .- theta) * (1.0 - model.sigma)
#             correction = nnn * model.sigma * sin(angles_neighbours[ID] - theta)
#             return force + correction
#         elseif model.symmetry == "nematic"
#             nnn2 = Int(nnn / 2)
#             base = sum(sin, 2(angles_neighbours .- theta)) * (1.0 - model.vision)
#             correction_front = nnn2 * model.vision * sin(angles_neighbours[ID] - theta)
#             correction_back = nnn2 * model.vision * sin(angles_neighbours[mod1(ID + nnn2, nnn)] - theta)
#             return base + correction_front + correction_back
#         end
#     end
# end

# Different angular weight distributions
visioncone(ϕ, viscone) = mod(ϕ, 2π) ≤ viscone / 2
vonmises(ϕ, sigma; ϕ0=0.0) = exp(sigma * cos(ϕ-ϕ0))  # ϕ0 = 0 looking ahead, π/2 to the left, -π/2 to the right
left_and_right(ϕ, sigma; ϕ0=0.0) = exp(sigma * sin(ϕ-ϕ0)^2) - sigma/2 # exotic distribution, not used
faf(ϕ,sigma; ϕ0=0.0) = sigma*cos(2(ϕ-ϕ0))*exp(-(ϕ-ϕ0)^2*sigma) # Ferro/AntiFerro 
# faf(ϕ,sigma; ϕ0=0.0) = sigma*cos(0.5(ϕ-ϕ0))*exp(-(ϕ-ϕ0)^2*sigma) # Ferro/AntiFerro 
# for tests
# ff(x) = faf(x,1); plot(x->ff(x),xlims=(-pi,pi),ylims=(-0.5,1.5),title=string(round(mean([ff(x) for x in -pi:0.01:pi]),digits=3)))


function sum_influence_neighbours(theta::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::NRXY{T}, lattice::Abstract2DLattice)::T where {T<:AbstractFloat}
    force = 0.0
    nnn = number_nearest_neighbours(lattice)
    for n in 1:nnn
        diff_angle = angles_neighbours[n] - theta
        phi = (n - 1) / nnn * 2π
        weight = exp(model.sigma * cos(theta - phi))
        force += weight * sin(diff_angle)
    end
    return force / nnn

    #= Normalisation by the number of neighbours, very important for 
    the critical temperature, which also gets divided by nnn ! 
    Tkt Square = 0.89/4 = 0.2225, Tkt Triangular = 1.47/6 = 0.245 =#


    # # # IMPLEMENTATION of the ghost NaN at the boundary 
    # force = 0.0
    # nnn = number_nearest_neighbours(lattice)
    # no_nans = count(isnan, angles_neighbours) == 0
    # if no_nans
    #     for n in 1:nnn
    #         diff_angle = angles_neighbours[n] - theta
    #         phi = theta - (n - 1) / nnn * 2π
    #         weight = vonmises(phi, model.sigma, ϕ0=0)
    #         force += weight * sin(diff_angle)
    #     end
    # else
    #     for n in 1:nnn
    #         if !isnan(angles_neighbours[n])
    #             diff_angle = angles_neighbours[n] - theta
    #             phi = theta - (n - 1) / nnn * 2π
    #             weight = vonmises(phi, model.sigma, ϕ0=0)
    #             force += weight * sin(diff_angle)
    #         end
    #     end
    # end
    # return force / (nnn)
end

function sum_influence_neighbours(theta::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::NRXY_ExtraTerm{T}, lattice::Abstract2DLattice)::T where {T<:AbstractFloat}
    force = 0.0
    nnn = number_nearest_neighbours(lattice)
    for n in 1:nnn
        diff_angle = angles_neighbours[n] - theta
        phi = (n - 1) / nnn * 2π
        weight = exp(model.sigma * cos(theta - phi))
        force += weight * (sin(diff_angle) - model.sigma * cos(diff_angle)*sin(theta-phi)) # to obtain the same dynamics as that of Sarah Loos.
    end
    return force / nnn

    #= Normalisation by the number of neighbours, very important for 
    the critical temperature, which also gets divided by nnn ! 
    Tkt Square = 0.89/4 = 0.2225, Tkt Triangular = 1.47/6 = 0.245 =#
end

# special for NRXY on SquareLattice, looses generality but 2/3 time ratio
function sum_influence_neighbours(theta::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::NRXY{T}, lattice::SquareLattice)::T where {T<:AbstractFloat}
    expcos = exp(model.sigma * cos(theta))
    expsin = exp(model.sigma * sin(theta))
    force = 0
    force += sin(angles_neighbours[1] - theta) * expcos 
    force += sin(angles_neighbours[2] - theta) * expsin 
    force += sin(angles_neighbours[3] - theta) / expcos 
    force += sin(angles_neighbours[4] - theta) / expsin
    return force/4
end

# special for NRXY on TriangularLattice, , looses generality but 3/4 time ratio
function sum_influence_neighbours(theta::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::NRXY{T}, lattice::TriangularLattice)::T where {T<:AbstractFloat}
    expcos = exp(model.sigma * cos(theta))
    expcos1 = exp(model.sigma * cos(theta - π/3))
    expcos2 = exp(model.sigma * cos(theta - 2π/3))
    force = 0
    force += sin(angles_neighbours[1] - theta) * expcos 
    force += sin(angles_neighbours[2] - theta) * expcos1
    force += sin(angles_neighbours[3] - theta) * expcos2
    force += sin(angles_neighbours[4] - theta) / expcos 
    force += sin(angles_neighbours[5] - theta) / expcos1
    force += sin(angles_neighbours[6] - theta) / expcos2
    return force/6
end

function sum_influence_neighbours(theta::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::HydrodynNRXY{T}, lattice::SquareLattice)::T where {T<:AbstractFloat}
    # Original implementation
    
    dthetadx = (angles_neighbours[1] - angles_neighbours[3]) / 2
    dthetady = (angles_neighbours[2] - angles_neighbours[4]) / 2
    dthetadx2 = (angles_neighbours[1] - 2*theta + angles_neighbours[3])
    dthetady2 = (angles_neighbours[2] - 2*theta + angles_neighbours[4])

    force_EL = 2(cos(dthetadx)*sin(0.5*dthetadx2) + cos(dthetady)*sin(0.5*dthetady2))
    force_NR = 2* model.sigma * ( cos(theta)*sin(dthetadx)*cos(0.5*dthetadx2) + sin(theta)*sin(dthetady)*cos(0.5*dthetady2) )

    nnn = 4 # we know this method is only defined for SquareLattice 
    return (force_EL + force_NR)/nnn
end

# function sum_influence_neighbours(theta::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::ContinuousNRNematicXY{T}, lattice::Abstract2DLattice)::T where {T<:AbstractFloat}
#     nnn = number_nearest_neighbours(lattice)
#     force = 0.0
#     for n in 1:nnn
#         diff_angle = angles_neighbours[n] - theta
#         angle_with_neighbour = (n - 1) / nnn * 2π - theta
#         # weight = vonmises(2angle_with_neighbour, 2model.sigma)
#         force += weight * sin(2diff_angle)
#     end
#     return force / nnn
# end

## ------------------------ Update Methods ------------------------ ##
## ------------------------ Update Methods ------------------------ ##
## ------------------------ Update Methods ------------------------ ##
## ------------------------ Update Methods ------------------------ ##
## ------------------------ Update Methods ------------------------ ##

function evolve!(thetas::AbstractArray{T}, model::AbstractModel, lattice::AbstractLattice, tmax) where {T<:AbstractFloat}
    while model.t < tmax
        update!(thetas, model, lattice)
    end
    return thetas
end

function evolve!(thetas::AbstractArray{T},phis::AbstractArray{T}, model::AbstractModel, lattice::AbstractLattice, tmax) where {T<:AbstractFloat}
    while model.t < tmax
        update!(thetas, phis, model, lattice)
    end
    return thetas
end

#= Optimization tips
Fastest : 
for angle in angle_neighbours
    sumsin += sin(angle - θ)
end

Intermediary
sum(sin,angle_neighbours .- θ)
sum(x->sin(x- θ),angle_neighbours) is slightly better
NB : we you can avoid multiplying arrays, better. 
    So: sum(x->sin(2x - θ),angle_neighbours) is way better than 
    sum(sin,2angle_neighbours .- 2θ)

Very slow
for n in each(angle_neighbours)
    sumsin += sin(angle_neighbours[n] - θ)
end


To be honest, I don't really understand why, since the following two functions runtimes are the same
xx = rand(10000)
function f1(xx::Array{Float64})
    sumsin = 0.0
    for angle in xx
        sumsin += sin(angle - 1.0)
    end
    return sumsin
end

function f2(xx::Array{Float64})
    sumsin = 0.0
    for i in eachindex(xx)
        sumsin += sin(xx[i] - 1.0)
    end
    return sumsin
end

f1(xx)
f2(xx)
@btime f1(xx)
@btime f2(xx)

VERY IMPORTANT : it would seem that "each" is by faaaar slower than "eachindex" (at least by a factor 10 )
=#

## ------------------------ Langevin-like Dynamics ------------------------ ##
## ------------------------ Langevin-like Dynamics ------------------------ ##

function update!(thetas::Matrix{<:FT}, model::Union{LangevinXY{FT},LangevinNematicXY{FT},VisionConeXY{FT},DiscreteNRXY{FT},NRXY{FT},NRXY_ExtraTerm{FT},ContinuousNRNematicXY{FT},HydrodynNRXY{FT}}, lattice::Abstract2DLattice) where {FT<:AbstractFloat}
    thetas_old = copy(thetas)
    L = lattice.L
    dt = model.dt
    T = model.T

    bulk = true
    Threads.@threads for j in 2:L-1 
        for i in 2:L-1
            θ = thetas_old[i, j]
            angle_neighbours = get_neighbours(thetas_old, model, lattice, i, j, bulk)
            thetas[i, j] = θ + dt * sum_influence_neighbours(θ, i, j, angle_neighbours, model, lattice) + sqrt(2T * dt) * randn(FT)
        end
    end
    
    if lattice.periodic # update the borders
        bulk = false
        for j in [1,L] 
            for i in 2:L-1
                θ = thetas_old[i, j]
                angle_neighbours = get_neighbours(thetas_old, model, lattice, i, j, bulk)
                thetas[i, j] = θ + dt * sum_influence_neighbours(θ, i, j, angle_neighbours, model, lattice) + sqrt(2T * dt) * randn(FT)
            end
        end
        for j in 1:L
            for i in [1,L]
                θ = thetas_old[i, j]
                angle_neighbours = get_neighbours(thetas_old, model, lattice, i, j, bulk)
                thetas[i, j] = θ + dt * sum_influence_neighbours(θ, i, j, angle_neighbours, model, lattice) + sqrt(2T * dt) * randn(FT)
            end
        end
    end

    model.t += dt
    return thetas
end

function update!(thetas::Matrix{<:FT}, model::Union{ForcedXY{FT},PropagationForcedXY{FT}}, lattice::Abstract2DLattice) where {FT<:AbstractFloat}
    thetas_old = copy(thetas)
    L = lattice.L
    dt = model.dt
    T = model.T

    if lattice.periodic
        range_ij = 1:L
    else
        range_ij = 2:L-1
    end
    for j in range_ij
        Threads.@threads for i in range_ij
            θ = thetas_old[i, j]
            angle_neighbours = get_neighbours(thetas_old, model, lattice, i, j, is_in_bulk(i, j, L))
            sumsin = 0.0
            for angle in angle_neighbours
                sumsin += sin(angle - θ)
            end
            thetas[i, j] = θ + dt * (model.omegas[i, j] + sumsin) + sqrt(2T * dt) * randn(FT)
        end
    end
    model.t += dt
    return thetas
end

function update!(thetas::Matrix{<:FT}, model::OUForcedXY{FT}, lattice::Abstract2DLattice) where {FT<:AbstractFloat}
    thetas_old = copy(thetas)
    L = lattice.L
    dt = model.dt
    T = model.T

    if lattice.periodic
        range_ij = 1:L
    else
        range_ij = 2:L-1
    end
    for j in range_ij
        Threads.@threads for i in range_ij
            θ = thetas_old[i, j]
            angle_neighbours = get_neighbours(thetas_old, model, lattice, i, j, is_in_bulk(i, j, L))
            sumsin = 0.0
            for angle in angle_neighbours
                sumsin += sin(angle - θ)
            end
            thetas[i, j] = θ + dt * (model.omegas[i, j] + sumsin) + sqrt(2T * dt) * randn(FT)
        end
    end
    #= Update the omegas, according to the Ornstein-Uhlenbeck process 
    The equation is ̇ω = -κ ω + σ_ω ξ(t) where ξ(t) is a white noise.
    The variance of the process at long times is given by <ω(t)²> = σ_ω^2 / 2κ
    Since one wants the variance of the process to be equal to the initial variance 
    of the distribution, one can choose σ_ω = σ, thus fixing the value of κ to 1/2.
    =#
    kappa = 0.5
    model.omegas -= kappa * model.omegas * dt + sqrt(model.Var * dt) * randn(FT,L,L)

    model.t += dt
    return thetas
end

## ------------------------ MonteCarlo Methods ------------------------ ##
## ------------------------ MonteCarlo Methods ------------------------ ##

function update!(thetas::Matrix{<:FT}, model::MonteCarloXY{FT}, lattice::Abstract2DLattice) where {FT<:AbstractFloat}
    thetas_old = copy(thetas)
    L = lattice.L
    T = model.T

    if lattice.periodic
        range_ij = 1:L
    else
        range_ij = 2:L-1
    end
    for j in range_ij
        for i in range_ij
            theta = thetas_old[i, j]
            proposal = mod(2sqrt(T) * randn() + theta, 2pi)
            # proposal = 2pi * rand()
            angle_neighbours = get_neighbours(thetas_old, model, lattice, i, j, is_in_bulk(i, j, L))
            sumsin = 0.0
            for angle in angle_neighbours
                sumsin += sin((proposal + theta) / 2 - angle)
            end
            dE = 2sin((proposal - theta) / 2) * sumsin
            if rand() < exp(-dE / T)
                @inbounds thetas[i, j] = proposal
            end
        end
    end

    model.t += 1
    return thetas
end

function update!(thetas::Matrix{<:FT}, model::MonteCarloNematicXY{FT}, lattice::Abstract2DLattice) where {FT<:AbstractFloat}
    thetas_old = copy(thetas)

    L = lattice.L
    T = model.T

    for j in 1:L
        for i in 1:L
            theta = thetas[i, j]
            proposal = mod(2sqrt(T) * randn() + theta, 2pi)
            angle_neighbours = get_neighbours(thetas_old, model, lattice, i, j, is_in_bulk(i, j, L))
            sumsin = 0.0
            for angle in angle_neighbours
                sumsin += sin(proposal + theta - 2angle)
            end
            @fastmath dE = sin(proposal - theta) * sumsin
            if rand() < exp(-dE / T)
                @inbounds thetas[i, j] = proposal
            end
        end
    end

    model.t += 1
    return thetas
end

function update!(thetas::Matrix{<:FT}, model::SPP{FT}, lattice::Abstract2DLattice) where {FT<:AbstractFloat}
    L = lattice.L
    order_trials = StatsBase.shuffle(1:L^2)
    for n in order_trials
        i, j = linear_to_square_index(n, L)
        in_bulk = 2 < i < L - 1 && 2 < j < L - 1 # si "seulement" 1 < x,y < L , bug pour avoir le voisin d'un voisin

        θ = thetas[i, j]
        if !isnan(θ) && (in_bulk || lattice.periodic) # translation : if on the border AND lattice not periodic, do nothing

            ic, jc = angle2neighbour(θ, i, j, model, lattice)
            θc = thetas[ic, jc]

            if isnan(θc) # destination empty, let's go
                thetas[ic, jc] = θ
                thetas[i, j] = NaN
            else # destination occupied
                collision!(thetas, model, lattice, (i, j), θ, (ic, jc), θc, in_bulk)
            end
        end
    end
    model.t += 1 # here = number of MonteCarlo steps
    return thetas
end


## ------------------------ NRXY Model on the TorusLattice ------------------------ ##
## ------------------------ NRXY Model on the TorusLattice ------------------------ ##
#= Here ψ is the equivalent of θ in the other models, ie. the field of angles.
θ is the angle from the exterior to the interior of the torus.
ϕ is the angle round the torus viewed from above. =# 

function update!(psis::Matrix{<:FT}, model::NRXY{FT}, lattice::TorusLattice) where {FT<:AbstractFloat}
    psis_old = copy(psis)

    L = lattice.L
    b = lattice.b
    dt = model.dt
    T = model.T

    range_ij = 2:L-1 
    inbulk = true
    Threads.@threads for j in range_ij  
        for i in range_ij
            ψ = psis_old[i, j]
            angle_neighbours = get_neighbours(psis_old, model, lattice, i, j, inbulk)
            @inbounds psis[i, j] += dt * sum_influence_neighbours(ψ, i, j, angle_neighbours, model, lattice) + sqrt(2T * dt) * randn(FT)
        end
    end
    inbulk = false
    for j in 1:L
        for i in [1,L]
            ψ = psis_old[i, j]
            angle_neighbours = get_neighbours(psis_old, model, lattice, i, j, inbulk)
            @inbounds psis[i, j] += dt * sum_influence_neighbours(ψ, i, j, angle_neighbours, model, lattice) + sqrt(2T * dt) * randn(FT)
        end
    end
    for j in [1,L]
        for i in 2:L-1
            ψ = psis_old[i, j]
            angle_neighbours = get_neighbours(psis_old, model, lattice, i, j, inbulk)
            @inbounds psis[i, j] += dt * sum_influence_neighbours(ψ, i, j, angle_neighbours, model, lattice) + sqrt(2T * dt) * randn(FT)
        end
    end
    
    

    model.t += dt
    return psis
end

# Optimized versions with matrices A, B, C, D to store the values of ̂n ⋅ (gα x gβ)
function sum_influence_neighbours(ψ::T, i::Int, j::Int, angles_neighbours::Vector{<:T}, model::NRXY{T}, lattice::TorusLattice{T})::T where {T<:AbstractFloat}
    L  = lattice.L
    b  = lattice.b
    σ = model.sigma
    
    A  = lattice.A_cross
    B  = lattice.B_cross
    C  = lattice.C_cross
    D  = lattice.D_cross
    
    # This function has to be highly optimized, we thus declare everything 
    sinpsi, cospsi = sincos(ψ)
    expcos = exp(σ*cospsi)
    expsin = exp(σ*sinpsi)
    sin1, cos1 = sincos(angles_neighbours[1])
    sin2, cos2 = sincos(angles_neighbours[2])
    sin3, cos3 = sincos(angles_neighbours[3])
    sin4, cos4 = sincos(angles_neighbours[4])

    force = 0 

    @fastmath force += cospsi * cos1 * A[1,i,j] * expcos 
    @fastmath force += cospsi * cos2 * A[2,i,j] * expsin
    @fastmath force += cospsi * cos3 * A[3,i,j] / expcos
    @fastmath force += cospsi * cos4 * A[4,i,j] / expsin

    @fastmath force += cospsi * sin1 * B[1,i,j] * expcos 
    @fastmath force += cospsi * sin2 * B[2,i,j] * expsin
    @fastmath force += cospsi * sin3 * B[3,i,j] / expcos
    @fastmath force += cospsi * sin4 * B[4,i,j] / expsin

    @fastmath force += sinpsi * cos1 * C[1,i,j] * expcos
    @fastmath force += sinpsi * cos2 * C[2,i,j] * expsin
    @fastmath force += sinpsi * cos3 * C[3,i,j] / expcos
    @fastmath force += sinpsi * cos4 * C[4,i,j] / expsin

    @fastmath force += sinpsi * sin1 * D[1,i,j] * expcos
    @fastmath force += sinpsi * sin2 * D[2,i,j] * expsin
    @fastmath force += sinpsi * sin3 * D[3,i,j] / expcos
    @fastmath force += sinpsi * sin4 * D[4,i,j] / expsin
    
    nnn = 4
    return force/nnn
end

## ------------------------ Update Hydrodynamic Vector Field Model (ft. Daniel Pearce) ------------------------ ##
## ------------------------ Update Hydrodynamic Vector Field Model (ft. Daniel Pearce) ------------------------ ##
## ------------------------ Update Hydrodynamic Vector Field Model (ft. Daniel Pearce) ------------------------ ##
## ------------------------ Update Hydrodynamic Vector Field Model (ft. Daniel Pearce) ------------------------ ##

function update!(Pthetas::Matrix{T}, Pphis::Matrix{T}, model::HydroVector{T}, lattice::TorusHydroLattice)::Tuple{Matrix{T},Matrix{T}} where {T<:AbstractFloat}
    Pthetas_old = copy(Pthetas)
    Pphis_old = copy(Pphis)

    d_theta_P_theta, d_phi_P_theta = first_derivatives(Pthetas)
    d_theta_P_phi, d_phi_P_phi = first_derivatives(Pphis)
    d_theta_P_theta *= -1  # because the θ axis is reversed wrt to the i indexing
    d_theta_P_phi *= -1

    d2_theta_P_theta, d2_phi_P_theta = second_derivatives(Pthetas)
    d2_theta_P_phi, d2_phi_P_phi = second_derivatives(Pphis)

    S2 = get_S2(Pthetas, Pphis, lattice)

    alpha = model.alpha
    L = lattice.L
    dt = model.dt

    for j in 1:L # φ
        for i in 1:L # θ
            theta = Float32(2π * (L - i) / L)
            rhoo = lattice.rho[i] # ρ = η + cos(θ) 
            rhoo_prime = lattice.rhop[i] # ρ' = - sin(θ)
            rhoo_second = -cos(theta) # ρ" = - cos(θ)

            Ptheta = Pthetas_old[i, j] # maybe change to Pthetas[i,j] for efficiency ?
            Pphi = Pphis_old[i, j] # maybe change to Pphis[i,j] for efficiency ?

            H_theta = 2 * rhoo * d2_theta_P_theta[i, j] +
                      2 * d2_phi_P_theta[i, j] / rhoo +
                      2 * rhoo_prime * d_theta_P_theta[i, j] -
                      4 * rhoo_prime * d_phi_P_phi[i, j] -
                      2 * rhoo_prime^2 / rhoo * Ptheta +
                      alpha * rhoo * (1 - S2[i, j]) * Ptheta

            H_phi = 2 * rhoo^3 * d2_theta_P_phi[i, j] +
                    2 * rhoo * d2_phi_P_phi[i, j] +
                    6 * rhoo^2 * rhoo_prime * d_theta_P_phi[i, j] +
                    4 * rhoo_prime * d_phi_P_theta[i, j] +
                    2 * rhoo^2 * rhoo_second * Pphi +
                    alpha * rhoo^3 * (1 - S2[i, j]) * Pphi

            Pthetas[i, j] += dt * H_theta # + sqrt(2 * model.T * dt) * randn(T)
            Pphis[i, j] += dt * H_phi # + sqrt(2 * model.T * dt) * randn(T)
        end
    end

    ############ Under matrix form = faster ?? ############
    # rhoo = 1
    # rhoo_prime  = 0
    # rhoo_second = 0
    # H_theta = 2*rhoo * d2_theta_P_theta + 2/rhoo * d2_phi_P_theta + alpha * rhoo * (1 .-S2)* Pthetas_old 
    # Pthetas += dt * H_theta 

    # H_phi = 2*rhoo^3 * d2_theta_P_phi + 2*rhoo * d2_phi_P_phi + alpha * rhoo^3 * (1 .-S2) * Pphis_old
    # Pphis += dt * H_phi  

    ############ ############ ############ ############


    model.t += dt
    return Pthetas, Pphis
end

function update!(Pthetas::Matrix{T}, Pphis::Matrix{T}, model::HydroTest{T}, lattice::TorusHydroLattice)::Tuple{Matrix{T},Matrix{T}} where {T<:AbstractFloat}
    Pthetas_old = copy(Pthetas)
    Pphis_old = copy(Pphis)

    d_theta_P_theta, d_phi_P_theta = first_derivatives(Pthetas)
    d_theta_P_phi, d_phi_P_phi = first_derivatives(Pphis)
    # d_theta_P_theta *= -1  # because the θ axis is reversed wrt to the i indexing
    # d_theta_P_phi *= -1

    d2_theta_P_theta, d2_phi_P_theta = second_derivatives(Pthetas)
    d2_theta_P_phi, d2_phi_P_phi = second_derivatives(Pphis)

    S2 = get_S2(Pthetas, Pphis, lattice)

    alpha = model.alpha
    sigma = model.sigma
    L = lattice.L
    dt = model.dt

    for j in 1:L # φ
        for i in 1:L # θ
            theta = Float32(2π * (L - i) / L)
            rhoo = lattice.rho[i] # ρ = η + cos(θ) 
            rhoo_prime = lattice.rhop[i] # ρ' = - sin(θ)
            rhoo_second = -cos(theta) # ρ" = - cos(θ)

            Ptheta = Pthetas_old[i, j] # maybe change to Pthetas[i,j] for efficiency ?
            Pphi = Pphis_old[i, j] # maybe change to Pphis[i,j] for efficiency ?

        kappa_p = rhoo * d_theta_P_phi[i,j] - d_phi_P_theta[i,j] / rhoo + 2*rhoo_prime * Pphi
        H_theta = 2 * rhoo * d2_theta_P_theta[i, j] +
                2 * d2_phi_P_theta[i, j] / rhoo +
                2 * rhoo_prime * d_theta_P_theta[i, j] +
                -4 * rhoo_prime * d_phi_P_phi[i, j] +
                -2 * rhoo_prime^2 / rhoo * Ptheta +
                sigma * kappa_p * rhoo * Pphi +
                alpha * rhoo * (1 - S2[i, j]) * Ptheta

        H_phi = 2 * rhoo^3 * d2_theta_P_phi[i, j] +
                2 * rhoo * d2_phi_P_phi[i, j] +
                6 * rhoo^2 * rhoo_prime * d_theta_P_phi[i, j] +
                4 * rhoo_prime * d_phi_P_theta[i, j] +
                2 * rhoo^2 * rhoo_second * Pphi +
                -sigma * kappa_p / rhoo * Ptheta +
                alpha * rhoo^3 * (1 - S2[i, j]) * Pphi

            Pthetas[i, j] += dt * H_theta # + sqrt(2 * model.T * dt) * randn(T)
            Pphis[i, j] += dt * H_phi # + sqrt(2 * model.T * dt) * randn(T)
        end
    end

    model.t += dt
    return Pthetas, Pphis
end

## ------------------------ Update Propagation Models ------------------------ ##
## ------------------------ Update Propagation Models ------------------------ ##
## ------------------------ Update Propagation Models ------------------------ ##
## ------------------------ Update Propagation Models ------------------------ ##

function update!(thetas::Vector{FT}, model::AbstractPropagationModel{FT}, lattice::Abstract1DLattice)::Vector{FT} where {FT<:AbstractFloat}
    thetas_old = copy(thetas)
    L = lattice.L
    dt = model.dt
    T = model.T

    model.symmetry == "polar" ? symm = 1 : symm = 2

    in_bulk = true
    for i in 2:L-1
        θ = thetas_old[i]
        angle_neighbours = get_neighbours(thetas_old, model, lattice, i, in_bulk)
        thetas[i] = θ + dt * (model.omegas[i] + sum(sin, symm * (angle_neighbours .- θ))) + sqrt(2T * dt) * randn(FT)
    end
    thetas[1] = thetas_old[1] + dt * (model.omegas[1] + sum(sin, symm * (get_neighbours(thetas_old, model, lattice, 1, false) .- thetas_old[1]))) + sqrt(2T * dt) * randn()
    thetas[L] = thetas_old[L] + dt * (model.omegas[L] + sum(sin, symm * (get_neighbours(thetas_old, model, lattice, L, false) .- thetas_old[L]))) + sqrt(2T * dt) * randn()

    model.t += dt
    return thetas
end

## ------------------------ Other Evolution Methods ------------------------ ##
## ------------------------ Other Evolution Methods ------------------------ ##

function collision!(thetas::Matrix{<:FT}, model::SPP{FT}, lattice::AbstractLattice, pos1::Tuple{T,T}, theta1::FT, pos2::Tuple{T,T}, theta2::FT, bulk::Bool) where {T<:Int,FT<:AbstractFloat}
    # @assert model.symmetry == "nematic" "Energy is only coded for nematic interaction for now !"
    width_proposal = 2sqrt(model.T)
    proposal = width_proposal * randn(FT) + theta1
    i, j = pos1
    if model.algo == "A" # Model A. Align nematically with all NN
        neighbours = 2get_neighbours(thetas, model, lattice, i, j, bulk)
        # dE = -1/2 * ( sum(cos, neighbours .- 2proposal ) - sum(cos, neighbours .- 2theta1 ))
        # dE = sin(proposal - theta1) * sum(sin,proposal + theta1 .- neighbours ) # should be computationally faster
        model.symmetry == "polar" ? coeff_symmetry = 1.0 : coeff_symmetry = 2.0
        coeff_symmetry2 = coeff_symmetry / 2.0
        @fastmath dE = 1.0 / coeff_symmetry2 * sin(coeff_symmetry2 * (proposal - theta1)) * sum(sin, coeff_symmetry2 * (proposal + theta1 .- neighbours))

    elseif model.algo == "B" # align nematically wrt collided spin
        dE = -1 / 2 * (cos(2(theta2 - proposal)) - cos(2(theta2 - theta1)))

    elseif model.algo == "C" # align F/AF wrt collided spin
        ccos = cos(theta1 - theta2)
        if ccos > 0
            J = +1.0 # ferro
        else
            J = -1.0 # antiferro
        end
        dE = -J * (cos(theta1 - proposal) - ccos)

    elseif model.algo == "CA"
        ccos = cos(theta1 - theta2)
        if ccos > 0
            J = +1.0 # ferro
        else
            J = -1.0 # antiferro
        end
        dE = -J * (cos(theta1 - proposal) - ccos)
        if rand() < exp(-dE / model.T) # Metropolis Monte Carlo
            @inbounds thetas[i, j] = proposal
        end
        proposal = model.width_proposal * randn(FT) + thetas[i, j] # new proposal
        neighbours = 2get_neighbours(thetas, model, lattice, i, j, bulk)
        dE = -1 / 2 * (sum(cos, neighbours .- 2proposal) - sum(cos, neighbours .- 2theta1))
        # the last Monte Carlo step, corresponding to the case "A", is performed by the last step (common to all cases)
    else
        error("Unknown Model")
    end
    if rand() < exp(-dE / model.T) # Metropolis Monte Carlo
        @inbounds thetas[i, j] = proposal
    end
    return thetas
end

function angle2neighbour(theta::AbstractFloat, i::Int, j::Int, model::AbstractModel, lattice::Abstract2DLattice)
    theta, A = Float64.((theta, model.A)) # Float64.((1,1)) 50x faster than Float64.([1,1])
    direction_motion = direction_of_motion(theta, A)
    NN = project_angle_onto_lattice(direction_motion, i, j, lattice)
    # if model.propulsion == "nematic" && rand(Bool)
    #     NN = (-1) .* NN # choose opposite direction with probability 1/2
    # end
    return add_2_positions((i, j), NN, lattice.L, true) # TODO, check whether false could make it and what runtime gain it would yields
end

function direction_of_motion(theta::T, A::T) where {T<:AbstractFloat}
    if A == 0
        angle = 2π * rand()
    else
        # angle = 1.0/sqrt(A)*randn()+theta  # Wrapped Normal
        # angle = rand(VonMises(theta,A)) # Von Mises, indistinguishable from Wrapped Normal for A > 4
        angle = mod(rand(Cauchy(theta, one(T) / A)), 2π) # Wrapped Cauchy, contractile activity
    end
    return angle
end
# histogram(mod.(rand(Cauchy(.5,0.25),Int(1E5)),2pi),normalize=true,bins=100)


function project_angle_onto_lattice(angle::AbstractFloat, i::Int, j::Int, lattice::Abstract2DLattice)
    nearest_neighbours = offsets(lattice, iseven(i)) # can be of length 4, 6 or 8
    nb_nn = length(nearest_neighbours)
    index_nearest_neighbour = round(Int, mod(angle, 2π) / (2π / nb_nn), RoundNearest) + 1
    if index_nearest_neighbour == nb_nn + 1
        index_nearest_neighbour = 1
    end
    return nearest_neighbours[index_nearest_neighbour]
end

# same as above but only returns the index instead of the offset to reach the said projected neighbour
function ID_projection_angle_onto_lattice(angle::AbstractFloat, i::Int, j::Int, lattice::Abstract2DLattice)::Int
    nearest_neighbours = offsets(lattice, iseven(i)) # can be of length 4, 6 or 8
    nb_nn = length(nearest_neighbours)
    index_nearest_neighbour = round(Int, mod(angle, 2π) / (2π / nb_nn), RoundNearest) + 1 # +1 so that theta=0 points to the spin on the
    if index_nearest_neighbour == nb_nn + 1
        index_nearest_neighbour = 1
    end
    return index_nearest_neighbour
end

# Meant to relax reconstruction for spotting defects or to relax the system at initialisation
function relax(thetas::Matrix{FT}, model::AbstractModel, lattice::Abstract2DLattice; trelax=2, T=0 ) where {FT<:AbstractFloat}
    thetas_copy = copy(thetas)
    dummy_dt = FT(1E-2)
    dummy_model = LangevinXY{FT}(T, model.symmetry, dummy_dt, zero(FT), model.rho)
    # dummy_lattice = SquareLattice(size(thetas_copy, 1), true, true, "chebychev")
    evolve!(thetas_copy, dummy_model, lattice, trelax)
    return thetas_copy
end
