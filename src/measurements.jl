include("lattices.jl");
include("models.jl");
include("core_methods.jl"); # for vonmises 

function OP(thetas::Matrix{T}) where T<:AbstractFloat
    tmp_polar   = Complex(0)
    tmp_nematic = Complex(0)
    for theta in thetas
        if !isnan(theta)
            tmp_polar   += exp(im*theta)
            tmp_nematic += exp(2im*theta)
        end
    end
    return abs(tmp_polar/length(thetas)) , abs(tmp_nematic/length(thetas))
end

function local_mag(thetas::Matrix{T},L::Int,l::Int)::Matrix{T} where T<:AbstractFloat
    # L is the size of the system
    # l is the size of the local box for averaging the magnetisation
    @assert isinteger(L/l)
    c = Int(L/l)
    δ = zeros(2,c^2)
    for j in 1:c
        for i in 1:c
            δ[:,square_to_linear_index(i,j,c)] = get_Delta(thetas[1+l*(i-1):i*l, 1+l*(j-1):j*l])
        end
    end
    return δ
end

# For SquareLattice and TriangularLattice
function corr_fft(thetas) 
    C = cos.(thetas)
    S = sin.(thetas)
    L = size(thetas,1)
    c_fft_2d = abs.((ifft(abs.(fft(C) .^ 2 + fft(S) .^ 2))))
    c_fft = (c_fft_2d[2:Int(end / 2), 1] + c_fft_2d[1, 2:Int(end / 2)]) / 2c_fft_2d[1, 1]
    # does the mean of the first row and first column and renormalises by the first element
    return c_fft
end
function corr(thetas::Matrix{T},model::AbstractModel,lattice::Abstract2DLattice) where T<:AbstractFloat
    L = lattice.L ; Lover2 = round(Int,L/2,RoundDown)
    matrix_corr = Matrix{T}(undef,Lover2,L^2)

    if lattice.periodic nmax = Lover2
    else 
        #TODO i ,j not defined
        error()
        #nmax = distance_to_border(thetas,i,j)-1
    end

    for m in 1:L^2
        i,j = linear_to_square_index(m,L)
        for n in 1:nmax
            matrix_corr[n,m] = corr(thetas,model,lattice,i,j,n)
        end
    end
    return nanmean(matrix_corr,2)[:,1]
end

# For SquareLattice and TriangularLattice
function corr(thetas::Matrix{T},model::AbstractModel,lattice::Abstract2DLattice,i::Int,j::Int,n::Int) where T<:AbstractFloat
    angle_neighbours_at_distance_n = [thetas[mod1(i+n,L),j],
                                      thetas[mod1(i-n,L),j],
                                      thetas[i,mod1(j+n,L)],
                                      thetas[i,mod1(j-n,L)] ]
    filter!(!isnan,angle_neighbours_at_distance_n)
    model.symmetry == "polar" ? symm = 1 : symm = 2
    if length(angle_neighbours_at_distance_n) > 0
        return mean(cos,symm*(angle_neighbours_at_distance_n .- thetas[i,j]))
    else
        return NaN
    end
end

# Method with correct_sites way too slow, so let's go back to a flawed but tractable solution : use the same algo for TriangularLattice and SquareLattice
# function corr(thetas::Matrix{T},model::AbstractModel,lattice::TriangularLattice,i::Int,j::Int,n::Int) where T<:AbstractFloat
#     # correct_sites = []
#     # for a in -n:+n , b in -n:+n # scann the square with chebychev distance ≤ n
#     #     first_selection = abs(a) + abs(b) ≥ n # to exclude already some inner part of the square
#     #     correct_distance = (dist(lattice,true_position(lattice,a,b),true_position(lattice,i,j)) == n)
#     #     if first_selection && correct_distance
#     #         push!(correct_sites,(a,b))
#     #     end
#     # end
#     angle_neighbours_at_distance_n = [thetas[mod1(i+n,L),j],
#                                       thetas[mod1(i-n,L),j],
#                                       thetas[i,mod1(j+n,L)],
#                                       thetas[i,mod1(j-n,L)],
#                                       thetas[mod1(i-n,L),mod1(j-n,L)],
#                                       thetas[mod1(i+n,L),mod1(j+n,L)],
#                                       thetas[mod1(i-n,L),mod1(j+n,L)],
#                                       thetas[mod1(i+n,L),mod1(j-n,L)] ]
#     posij = true_position(lattice,i,j)
#     filter!(x->(dist(lattice,true_position(lattice,x...),posij) == n),angle_neighbours_at_distance_n)
#     filter!(!isnan,angle_neighbours_at_distance_n)
#     model.symmetry == "polar" ? symm = 1 : symm = 2
#     if length(angle_neighbours_at_distance_n) > 0
#         return mean(cos,symm*(angle_neighbours_at_distance_n .- thetas[i,j]))
#     else
#         return NaN
#     end
# end

# For TorusLattice, the ϕ and θ direction are not equivalent so the result cannot be averaged over both
function corr(thetas::Matrix{T},model::AbstractModel,lattice::AbstractTorusLattice) where T<:AbstractFloat
    Lover2 = round(Int, L/2,RoundDown)
    corr_along_theta = zeros(T,L,L,Lover2+1,2) # Last index : 1 for positive distance
    corr_along_phi = zeros(T,L,L,Lover2+1,2)   # Last index : 1 for positive distance rφ, 2 for rθ, 2 for negative distance

    for j in 1:L # φ direction
        for i in 1:L # θ direction
            theta = thetas[i,j]
            for k in 0:Lover2
                @inbounds corr_along_theta[i,j,1+k,1] = cos(thetas[mod1(i+k,L),j] - theta) 
                @inbounds corr_along_theta[i,j,1+k,2] = cos(thetas[mod1(i-k,L),j] - theta)
                @inbounds corr_along_phi[i,j,1+k,1]   = cos(thetas[i,mod1(j+k,L)] - theta)
                @inbounds corr_along_phi[i,j,1+k,2]   = cos(thetas[i,mod1(j-k,L)] - theta)
            end
        end
    end
    corr_along_theta_for_each_theta = reverse(mean(corr_along_theta, dims=2),dims=1)[:,1,:,:] # averaged over φ, reversed to obtain θ increasing from 0 to 2π  
    corr_along_phi_for_each_theta   = reverse(mean(corr_along_phi, dims=2),dims=1)[:,1,:,:]   # averaged over θ, reversed to obtain θ increasing from 0 to 2π 
    return corr_along_theta_for_each_theta, corr_along_phi_for_each_theta
end

function corr_length(C::Vector{T},rs=1:length(C);seuil=exp(-1))::T where T<:AbstractFloat # from a time series, returns the correlation length ()
    i_after = findfirst(x->x<seuil,C)
    if i_after ≠ nothing && i_after > 1
        # Linear interpolation
        i_befor = i_after - 1
        r_after = rs[i_after]
        r_befor = rs[i_befor]
        c_after = C[i_after]
        c_befor = C[i_befor]
        ξ = (seuil*(r_after-r_befor) -  (c_befor*r_after - r_befor*c_after))/(c_after-c_befor)
    else
    ξ = NaN
    end
    return ξ
end


function energy(thetas::Matrix{T},model::AbstractModel,lattice::Abstract2DLattice) where T<:AbstractFloat
    L = size(thetas,1)
    E = 0
    if lattice.periodic
        range_ij = 1:L
    else
        range_ij = 2:L-1
    end
    for j in range_ij
        for i in range_ij
            angle_neighbours = get_neighbours(thetas,model,lattice,i,j)
            E += energy(thetas[i,j],angle_neighbours,model,lattice)
        end
    end
    nnn = number_nearest_neighbours(lattice)
    return E/(L^2)/2 # because we counted each pair twice
end


function energy(theta::T,angle_neighbours::Vector{T},model::AbstractModel,lattice::Abstract2DLattice) where T<:AbstractFloat
    energy = 0.0
    symm_coeff = model.symmetry == "polar" ? 1 : 2
    for angle in angle_neighbours
        energy -= cos(symm_coeff*(angle - theta))
    end
    energy = 1/symm_coeff * energy
    return energy
end

function energy(theta::T,angle_neighbours::Vector{T},model::NRXY,lattice::Abstract2DLattice) where T<:AbstractFloat
    energy = 0.0
    nnn = length(angle_neighbours)
    for n in 1:nnn
        phi = theta - (n-1)*2π/nnn
        w = vonmises(phi, model.sigma) # exp(σ*cos(φ))
        energy -= w * cos(angle_neighbours[n] - theta)
    end
    return energy
end


## ---------------------------- On the TorusLattice ---------------------------- ##
## ---------------------------- On the TorusLattice ---------------------------- ##
## ---------------------------- On the TorusLattice ---------------------------- ##
## ---------------------------- On the TorusLattice ---------------------------- ##
function energy(psis::Matrix{T},model::AbstractModel,lattice::TorusLattice)::T where T<:AbstractFloat
    L = size(thetas,1)
    E = 0

    A = lattice.A_dot
    B = lattice.B_dot
    C = lattice.C_dot
    D = lattice.D_dot

    for j in 1:L
        for i in 1:L
            ψ = psis[i,j] ; sinpsi, cospsi = sin(ψ), cos(ψ)
            angle_neighbours = get_neighbours(psis, model, lattice, i,j)        

            ## right neighbour (θ, φ+2π/L)
            E -= cospsi*cos(angle_neighbours[1]) * A[1,i,j]
            E -= cospsi*sin(angle_neighbours[1]) * B[1,i,j]
            E -= sinpsi*cos(angle_neighbours[1]) * C[1,i,j]
            E -= sinpsi*sin(angle_neighbours[1]) * D[1,i,j]

            ## top neighbour (θ+2π/L, φ)
            E -= cospsi*cos(angle_neighbours[2]) * A[2,i,j]
            E -= cospsi*sin(angle_neighbours[2]) * B[2,i,j]
            E -= sinpsi*cos(angle_neighbours[2]) * C[2,i,j]
            E -= sinpsi*sin(angle_neighbours[2]) * D[2,i,j]

            ## left neighbour (θ, φ-2π/L)
            E -= cospsi*cos(angle_neighbours[3]) * A[3,i,j]
            E -= cospsi*sin(angle_neighbours[3]) * B[3,i,j]
            E -= sinpsi*cos(angle_neighbours[3]) * C[3,i,j]
            E -= sinpsi*sin(angle_neighbours[3]) * D[3,i,j]

            ## bottom neighbour (θ-2π/L, φ)
            E -= cospsi*cos(angle_neighbours[4]) * A[4,i,j]
            E -= cospsi*sin(angle_neighbours[4]) * B[4,i,j]
            E -= sinpsi*cos(angle_neighbours[4]) * C[4,i,j]
            E -= sinpsi*sin(angle_neighbours[4]) * D[4,i,j]

        end
    end
    return E/(L^2)/2 # because we counted each pair twice
end

