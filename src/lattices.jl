## ------------------------ Lattices ------------------------
abstract type AbstractLattice end
abstract type Abstract1DLattice <: AbstractLattice end
abstract type Abstract2DLattice <: AbstractLattice end
abstract type AbstractTorusLattice <: AbstractLattice end

mutable struct Chain1D <: Abstract1DLattice
    L::Int
    periodic::Bool
end
Chain1D(L::Int;periodic::Bool=true) = Chain1D(L,periodic)

mutable struct TriangularLattice <: Abstract2DLattice
    L::Int
    periodic::Bool
    single::Bool
    metric::String
end
TriangularLattice(L::Int;periodic::Bool=true,single::Bool=true,metric::String="euclidian") = TriangularLattice(L,periodic,single,metric)

mutable struct SquareLattice <: Abstract2DLattice
    L::Int
    periodic::Bool
    single::Bool
    metric::String
end
SquareLattice(L::Int;periodic::Bool=true,single::Bool=true,metric::String="euclidian") = SquareLattice(L,periodic,single,metric)

function number_nearest_neighbours(lattice::Abstract2DLattice)
    if isa(lattice,TriangularLattice) nnn = 6
    elseif isa(lattice,SquareLattice) nnn = 4
    end
    return nnn
end

# Finite difference derivatives with periodic boundary conditions
# function first_derivatives_order1(P::Matrix{T})::Tuple{Matrix{T},Matrix{T}} where {T<:AbstractFloat}
#     L = size(P,1)
#     dP_dx = similar(P) # zeros(T,L,L)
#     dP_dy = similar(P) # zeros(T,L,L)
#     for j in 1:L 
#         for i in 2:L-1
#             @inbounds dP_dx[i,j] = (P[i+1,j] - P[i-1,j])/2
#         end
#         @inbounds dP_dx[1,j] = (P[2,j] - P[L,j])/2
#         @inbounds dP_dx[L,j] = (P[1,j] - P[L-1,j])/2
#     end

#     # y direction
#     for j in 2:L-1
#         for i in 1:L 
#             @inbounds dP_dy[i,j] = (P[i,j+1] - P[i,j-1])/2
#         end
#     end
#     for i in 1:L 
#         @inbounds dP_dy[i,1] = (P[i,2] - P[i,L])/2
#         @inbounds dP_dy[i,L] = (P[i,1] - P[i,L-1])/2
#     end
#     Delta_x = 2π/L
#     Delta_y = 2π/L
#     return dP_dx/Delta_x, dP_dy/Delta_y 
# end
function first_derivatives(P::Matrix{T})::Tuple{Matrix{T},Matrix{T}} where {T<:AbstractFloat}
    L = size(P,1)
    dP_dx = similar(P) # zeros(T,L,L)
    dP_dy = similar(P) # zeros(T,L,L)
    for j in 1:L 
        for i in 3:L-2
            @inbounds dP_dx[i,j] = (1/12*P[i-2,j] - 2/3*P[i-1,j] +2/3*P[i+1,j] -1/12*P[i+2,j])/2
        end
        @inbounds dP_dx[1,j]   = (1/12*P[L-1,j] - 2/3*P[L,j]   +2/3*P[2,j] -1/12*P[3,j])/2
        @inbounds dP_dx[2,j]   = (1/12*P[L,j]   - 2/3*P[1,j]   +2/3*P[3,j] -1/12*P[4,j])/2
        @inbounds dP_dx[L,j]   = (1/12*P[L-2,j] - 2/3*P[L-1,j] +2/3*P[1,j] -1/12*P[2,j])/2
        @inbounds dP_dx[L-1,j] = (1/12*P[L-3,j] - 2/3*P[L-2,j] +2/3*P[L,j] -1/12*P[1,j])/2
    end

    # y direction
    for j in 3:L-2
        for i in 1:L 
            @inbounds dP_dy[i,j] = (1/12*P[i,j-2] - 2/3*P[i,j-1] +2/3*P[i,j+1] -1/12*P[i,j+2])/2
        end
    end
    for i in 1:L 
        @inbounds dP_dy[i,1]   = (1/12*P[i,L-1] - 2/3*P[i,L]   +2/3*P[i,2] -1/12*P[i,3])/2
        @inbounds dP_dy[i,2]   = (1/12*P[i,L]   - 2/3*P[i,1]   +2/3*P[i,3] -1/12*P[i,4])/2
        @inbounds dP_dy[i,L-1] = (1/12*P[i,L-3] - 2/3*P[i,L-2] +2/3*P[i,L] -1/12*P[i,1])/2
        @inbounds dP_dy[i,L]   = (1/12*P[i,L-2] - 2/3*P[i,L-1] +2/3*P[i,1] -1/12*P[i,2])/2
    end

    Delta_x = 2π/L
    Delta_y = 2π/L
    return dP_dx/Delta_x, dP_dy/Delta_y 
end

# function second_derivatives_order1(P::Matrix{T})::Tuple{Matrix{T},Matrix{T}} where {T<:AbstractFloat}
#     L = size(P,1)
#     ddP_x = similar(P)
#     ddP_y = similar(P)
#     # x direction
#     for j in 1:L 
#         for i in 2:L-1
#             @inbounds ddP_x[i,j] = (P[i+1,j] - 2*P[i,j] + P[i-1,j])
#         end
#         @inbounds ddP_x[1,j] = (P[2,j] - 2*P[1,j] + P[L,j])
#         @inbounds ddP_x[L,j] = (P[1,j] - 2*P[L,j] + P[L-1,j])
#     end

#     # y direction
#     for j in 2:L-1
#         for i in 1:L 
#             @inbounds ddP_y[i,j] = (P[i,j+1] - 2*P[i,j] + P[i,j-1])
#         end
#     end
#     for i in 1:L 
#         @inbounds ddP_y[i,1] = (P[i,2] - 2*P[i,1] + P[i,L])
#         @inbounds ddP_y[i,L] = (P[i,1] - 2*P[i,L] + P[i,L-1])
#     end
#     Delta_x = 2π/L
#     Delta_y = 2π/L
#     return ddP_x/Delta_x^2, ddP_y/Delta_y^2
#     # return ddP_x, ddP_y
# end

function second_derivatives(P::Matrix{T})::Tuple{Matrix{T},Matrix{T}} where {T<:AbstractFloat}
    L = size(P,1)
    ddP_x = similar(P)
    ddP_y = similar(P)
    # x direction
    for j in 1:L 
        for i in 3:L-2
            @inbounds ddP_x[i,j] = (-1/12*P[i+2,j] + 4/3*P[i+1,j] - 5/2*P[i,j] + 4/3*P[i-1,j] +  -1/12*P[i-2,j])
        end
        @inbounds ddP_x[1,j]   = (-1/12*P[3,j] + 4/3*P[2,j] - 5/2*P[1,j] + 4/3*P[L,j] +  -1/12*P[L-1,j])
        @inbounds ddP_x[2,j]   = (-1/12*P[4,j] + 4/3*P[3,j] - 5/2*P[2,j] + 4/3*P[1,j] +  -1/12*P[L,j])
        @inbounds ddP_x[L-1,j] = (-1/12*P[1,j] + 4/3*P[L,j] - 5/2*P[L-1,j] + 4/3*P[L-2,j] +  -1/12*P[L-3,j])
        @inbounds ddP_x[L,j]   = (-1/12*P[2,j] + 4/3*P[1,j] - 5/2*P[L,j] + 4/3*P[L-1,j] +  -1/12*P[L-2,j])
    end

    # y direction
    for j in 3:L-2
        for i in 1:L 
            @inbounds ddP_y[i,j] = (-1/12*P[i,j+2] + 4/3*P[i,j+1] - 5/2*P[i,j] + 4/3*P[i,j-1] +  -1/12*P[i,j-2])
        end
    end
    for i in 1:L 
        @inbounds ddP_y[i,1]   = (-1/12*P[i,3] + 4/3*P[i,2] - 5/2*P[i,1] + 4/3*P[i,L] +  -1/12*P[i,L-1])
        @inbounds ddP_y[i,2]   = (-1/12*P[i,4] + 4/3*P[i,3] - 5/2*P[i,2] + 4/3*P[i,1] +  -1/12*P[i,L])
        @inbounds ddP_y[i,L-1] = (-1/12*P[i,1] + 4/3*P[i,L] - 5/2*P[i,L-1] + 4/3*P[i,L-2] +  -1/12*P[i,L-3])
        @inbounds ddP_y[i,L]   = (-1/12*P[i,2] + 4/3*P[i,1] - 5/2*P[i,L] + 4/3*P[i,L-1] +  -1/12*P[i,L-2])
    end

   
    Delta_x = 2π/L
    Delta_y = 2π/L
    return ddP_x/Delta_x^2, ddP_y/Delta_y^2
    # return ddP_x, ddP_y
end

## ------------------------ Torus Lattices ------------------------ ##
## ------------------------ Torus Lattices ------------------------ ##
## ------------------------ Torus Lattices ------------------------ ##
## ------------------------ Torus Lattices ------------------------ ##


number_nearest_neighbours(lattice::AbstractTorusLattice) = 4
mutable struct TorusLattice{T} <: AbstractTorusLattice  where {T<:AbstractFloat}
    L::Int # number of nodes in the theta and phi directions (the lattice is LxL)
    b::T # the aspect ratio of the torus, distance between the center of the torus and the center of the tube = 1, b<1 is the radius of the tube. 
   
    # Precomputation of the dot product of the cross products gα x gβ with the normal : n ⋅ (gα x gβ)
    A_cross::Array{T,3}
    B_cross::Array{T,3}
    C_cross::Array{T,3}
    D_cross::Array{T,3}
    
    # Precomputation of the dot products gα ⋅ gβ
    A_dot::Array{T,3}
    B_dot::Array{T,3}
    C_dot::Array{T,3}
    D_dot::Array{T,3}
end
function TorusLattice(L::Int, b=0.5) 
    A_cross, B_cross, C_cross, D_cross = matricesABCD_cross(L)
    A_dot, B_dot, C_dot, D_dot = matricesABCD_dot(L)
    return TorusLattice(L,Float32(b), A_cross, B_cross, C_cross, D_cross, A_dot, B_dot, C_dot, D_dot)
end

function matricesABCD_cross(L::Int)::Tuple{Array{Float32,3}, Array{Float32,3}, Array{Float32,3}, Array{Float32,3}} 
    A = zeros(Float32,4,L,L) # dot product of the normal ̂n ⋅ gtheta_X_gtheta
    B = zeros(Float32,4,L,L) # dot product of the normal ̂n ⋅ gtheta_X_gphi
    C = zeros(Float32,4,L,L) # dot product of the normal ̂n ⋅ gphi_X_gtheta
    D = zeros(Float32,4,L,L) # dot product of the normal ̂n ⋅ gphi_X_gphi
      
    for z in 1:4 
        if z == 1 # right
            dtheta = 0
            dphi   = Float32(2π/L) 
        elseif z == 2 # top
            dtheta = Float32(2π/L)
            dphi   = 0
        elseif z == 3 # left
            dtheta = 0
            dphi   = -Float32(2π/L)
        elseif z == 4 # bottom
            dtheta = -Float32(2π/L)
            dphi   = 0
        end

        for j in 1:L
            for i in 1:L 
                phi = Float32(2π * (j - 1) / L)
                theta = Float32(2π * (L - i) / L)

                nx = -cos(theta)*cos(phi)
                ny = -cos(theta)*sin(phi)
                nz = -sin(theta)

                Ax = -sin(theta) * cos(theta + dtheta) * sin(phi) + sin(theta + dtheta) * cos(theta) * sin(phi + dphi)
                Ay = -cos(theta) * sin(theta + dtheta) * cos(phi + dphi) + cos(theta + dtheta) * sin(theta) * cos(phi)
                Az = sin(theta)*sin(theta + dtheta) * sin(dphi)
                A[z,i,j] = nx*Ax + ny*Ay + nz*Az

                Bx = - cos(theta) * cos(phi + dphi)
                By = - cos(theta) * sin(phi + dphi)
                Bz = - sin(theta) * cos(dphi)
                B[z,i,j] = nx*Bx + ny*By + nz*Bz

                Cx = cos(theta + dtheta) * cos(phi)
                Cy = cos(theta + dtheta) * sin(phi)
                Cz = sin(theta + dtheta) * cos(dphi)
                C[z,i,j] = nx*Cx + ny*Cy + nz*Cz

                D[z,i,j] = nz * sin(dphi) 
            end
        end
    end

    return A, B, C, D
end

# # Cross products gα x gβ
function gtheta_X_gtheta(theta1::T, theta2::T, phi1::T, phi2::T, rho1::T, rho2::T, b::T)::Vector{T} where {T<:AbstractFloat}
    g1 = -sin(theta1) * cos(theta2) * sin(phi1) + sin(theta2) * cos(theta1) * sin(phi2)
    g2 = -cos(theta1) * sin(theta2) * cos(phi2) + cos(theta2) * sin(theta1) * cos(phi1) 
    g3 = sin(theta1)*sin(theta2) * sin(phi2 - phi1) # optimised version with trigo identities
    return b^2 * [g1, g2, g3]
end

function gtheta_X_gphi(theta1::T, theta2::T, phi1::T, phi2::T, rho1::T, rho2::T, b::T)::Vector{T} where {T<:AbstractFloat}
    g1 = cos(theta1) * cos(phi2) 
    g2 = cos(theta1) * sin(phi2) 
    g3 = sin(theta1) * cos(phi2-phi1) # optimized version with trigo identities
    return -b*rho2*[g1, g2, g3]
end

function gphi_X_gtheta(theta1::T, theta2::T, phi1::T, phi2::T, rho1::T, rho2::T, b::T)::Vector{T} where {T<:AbstractFloat}
    g1 = cos(theta2) * cos(phi1) 
    g2 = cos(theta2) * sin(phi1) 
    g3 = sin(theta2) * cos(phi2-phi1) # optimized version with trigo identities
    return b*rho1*[g1, g2, g3]
end

function gphi_X_gphi(theta1::T, theta2::T, phi1::T, phi2::T, rho1::T, rho2::T, b::T)::Vector{T} where {T<:AbstractFloat}
    g1 = 0
    g2 = 0
    g3 = rho1 * rho2 * sin(phi1 - phi2) # optimized version with trigo identities
    # g3 = rho1 * rho2 * sin(phi1 - phi2) # optimized version with trigo identities
    return [g1, g2, g3]
end

function matricesABCD_dot(L::Int)::Tuple{Array{Float32,3}, Array{Float32,3}, Array{Float32,3}, Array{Float32,3}} 
    A = zeros(Float32,4,L,L) # dot product gtheta ⋅ gtheta
    B = zeros(Float32,4,L,L) # dot product gtheta ⋅ gphi
    C = zeros(Float32,4,L,L) # dot product gphi ⋅ gtheta
    D = zeros(Float32,4,L,L) # dot product gphi ⋅ gphi
      
    for z in 1:4 
        if z == 1 # right
            dtheta = 0
            dphi   = Float32(2π/L) 
        elseif z == 2 # top
            dtheta = Float32(2π/L)
            dphi   = 0
        elseif z == 3 # left
            dtheta = 0
            dphi   = -Float32(2π/L)
        elseif z == 4 # bottom
            dtheta = -Float32(2π/L)
            dphi   = 0
        end

        for j in 1:L
            for i in 1:L 
                phi = Float32(2π * (j - 1) / L)
                theta = Float32(2π * (L - i) / L)

                A[z,i,j] = sin(theta) * sin(theta + dtheta) * cos(phi) * cos(phi + dphi) + sin(theta) * sin(theta + dtheta) * sin(phi) * sin(phi + dphi) + cos(theta) * cos(theta + dtheta)
                B[z,i,j] = sin(theta) * cos(phi) * sin(phi + dphi) - sin(theta) * sin(phi) * cos(phi + dphi)
                C[z,i,j] = sin(theta + dtheta) * cos(phi + dphi) * sin(phi) - sin(theta + dtheta) * sin(phi + dphi) * cos(phi)
                D[z,i,j] = sin(phi) * sin(phi + dphi) + cos(phi) * cos(phi + dphi)  
            end
        end
    end

    return A, B, C, D
end

function gtheta_dot_gtheta(theta1::T, theta2::T, phi1::T, phi2::T, rho1::T, rho2::T, b::T)::T where {T<:AbstractFloat}
    return b^2 * (sin(theta1) * sin(theta2) * cos(theta2 - theta1) + cos(theta1) * cos(theta2) ) 
end

function gtheta_dot_gphi(theta1::T, theta2::T, phi1::T, phi2::T, rho1::T, rho2::T, b::T)::T where {T<:AbstractFloat}
    return b * rho2 * sin(theta1) * sin(phi2 - phi1)
end

function gphi_dot_gtheta(theta1::T, theta2::T, phi1::T, phi2::T, rho1::T, rho2::T, b::T)::T where {T<:AbstractFloat}
    return b * rho1 * sin(theta2) * sin(phi1 - phi2)
end

function gphi_dot_gphi(theta1::T, theta2::T, phi1::T, phi2::T, rho1::T, rho2::T, b::T)::T where {T<:AbstractFloat}
    return rho1 * rho2 * cos(phi1 - phi2)
end

# for the continuous hydrodynamic model ft. Daniel Pearce
mutable struct TorusHydroLattice{T} <: AbstractTorusLattice  where {T<:AbstractFloat}
    L::Int # number of nodes in the theta and phi directions (the lattice is LxL)
    eta::T # η > 1 is the distance from the center of the hole to the center of the annulus, where the annulus has radius 1. 
    rho    # ρ  = η + cos(θ)
    rhop   # ρ' = dρ/dθ= -sin(θ)
    rhopp  # ρ" = d²ρ/dθ²= -cos(θ)
end
function TorusHydroLattice(L::Int, eta) 
    rho, rhop, rhopp = rho_rhoprime_rhosecond_cpu(L, eta)
    return TorusHydroLattice(L, eta, rho, rhop, rhopp)
end

function rho_rhoprime_rhosecond_cpu(L,eta)
    if eta == 0 || eta == Inf return ones(Float32,L), zeros(Float32,L), zeros(Float32,L) end # to knock out curvature
    rho   = zeros(Float32,L,L)
    rhop  = zeros(Float32,L,L)
    rhopp = zeros(Float32,L,L)
    for i in 1:L
        theta   = Float32(2π * (L - i) / L)
        rho[i,:]   .= eta + cos(theta)
        rhop[i,:]  .= -sin(theta)
        rhopp[i,:] .= -cos(theta)
    end
    return rho, rhop, rhopp
end

"""
    get_psis_S2(thetas::Matrix{T}, phis::Matrix{T}, lattice::TorusHydroLattice) where {T<:AbstractFloat}

The theta and phi components of the vector field are defined as follows: 
    
        Pθ = S cos(ψ)
        Pφ = S sin(ψ) / ρ, where ρ = η + cos(θ) and ψ is the angle of the vector field.

        Therefore, to inverse the relations and obtain ψ,S² from Pθ and Pφ, we have:

        ψ = atan(Pφ,Pθ/ρ)
        S² = Pθ² + (Pφ/ρ)²

        Example of usage: psis, S2 = get_psis_S2(thetas, phis, lattice)
"""
function get_psis_S2(thetas::Matrix{T}, phis::Matrix{T}, lattice::TorusHydroLattice)::Tuple{Matrix{T},Matrix{T}} where {T<:AbstractFloat}
    L = size(thetas,1)
    S2 = zeros(T,L,L)
    psis = zeros(T,L,L)
    for j in 1:L 
        for i in 1:L
            ρ = lattice.rho[i]
            @inbounds psis[i,j] = atan(phis[i,j],thetas[i,j]/ρ)
            @inbounds S2[i,j] = thetas[i,j]^2 + (phis[i,j]/ρ)^2
        end
    end

    return mod.(psis,2π), S2
end

function get_psis(thetas::Matrix{T}, phis::Matrix{T}, lattice::TorusHydroLattice)::Matrix{T} where {T<:AbstractFloat}
    L = size(thetas,1)
    psis = zeros(T,L,L)
    for j in 1:L 
        for i in 1:L
            θ = Float32(2π * (L - i) / L)
            ρ = lattice.rho[i]
            @inbounds psis[i,j] = atan(phis[i,j],thetas[i,j]/ρ)
        end
    end

    return psis
end

function get_S2(thetas::Matrix{T}, phis::Matrix{T}, lattice::TorusHydroLattice)::Matrix{T} where {T<:AbstractFloat}
    L = size(thetas,1)
    S2 = zeros(T,L,L)
    for j in 1:L 
        for i in 1:L
            θ = Float32(2π * (L - i) / L)
            ρ = lattice.rho[i]
            @inbounds S2[i,j] = thetas[i,j]^2 + (phis[i,j]/ρ)^2
        end
    end

    return S2
end

## ------------------------ Functions ------------------------
function dist(lattice::Abstract2DLattice,pos1,pos2)
    a,b = pos1 ; x,y = pos2
    dx = abs(x-a)
    dy = abs(y-b)
    if lattice.periodic
        dx = min(dx,lattice.L-dx)
        dy = min(dy,lattice.L-dy)
    end
    if     lattice.metric == "euclidian" return sqrt(dx^2 + dy^2)
    elseif lattice.metric == "manhattan" return dx + dy
    elseif lattice.metric == "chebychev" return max(dx,dy)
    else error("Unknown metric !")
    end
end

function dist2(lattice::Abstract2DLattice,pos1,pos2)
    a,b = pos1 ; x,y = pos2
    dx = abs(x-a)
    dy = abs(y-b)
    if lattice.periodic
        dx = min(dx,lattice.L-dx)
        dy = min(dy,lattice.L-dy)
    end
    if     lattice.metric == "euclidian" return dx^2 + dy^2
    elseif lattice.metric == "manhattan" return (dx + dy)^2
    elseif lattice.metric == "chebychev" return max(dx,dy)^2
    else error("Unknown metric !")
    end
end

dist(lattice::TorusLattice,pos1,pos2) = sqrt(dist2(lattice,pos1,pos2))
function dist2(lattice::TorusLattice,pos1,pos2)
    a,b = pos1 ; x,y = pos2
    dx = abs(x-a)
    dy = abs(y-b)
    dx = min(dx,lattice.L-dx)
    dy = min(dy,lattice.L-dy)
    return dx^2 + dy^2
end


function distance_matrix(new,old,lattice::AbstractLattice)
    m_new,m_old = length(new),length(old)
    distance_matrix = zeros(m_new,m_old)
    for j in 1:m_old
        for i in 1:m_new
            distance_matrix[i,j] = dist(lattice,new[i],old[j])
        end
    end
    return distance_matrix
end

#= By convention, the first neighbour = the first element of the vector "offsets"
is the one on the right =#

function offsets(lattice::TriangularLattice,even::Bool)::Vector{Tuple{Int,Int}}
    if even return [(0,1)  , (-1,1)  , (-1,0) ,  (0,-1), (1,0) , (1,1) ]
    else    return [(0,1) , (-1,0) , (-1,-1) , (0,-1) , (1,-1), (1,0) ]
    end
end

function offsets(lattice::SquareLattice,even::Bool)::Vector{Tuple{Int,Int}}
    return [(0,1) , (-1,0) , (0,-1) , (1,0)]

    # the following is unused
    # if lattice.metric in ["manhattan","euclidian"]
    #     return [(1,0), (0,1) , (-1,0) , (0,-1) ]
    # elseif lattice.metric =="chebychev"
    #     return [(1,0) , (1,1) , (0,1) , (-1,1) , (-1,0) , (-1,-1) , (0,-1) , (1,-1)]
    # else error("Unknown metric !")
    # end
end

function linear_to_square_index(n::Int,L::Int)
    i = div(n-1,L) + 1  # 1 ≤ i ≤ L
    j = mod1(n,L)       # 1 ≤ j ≤ L
    #= formula for reversing the linear indexing. i is the quotient and j
    the reminder of the euclidian division of n by L and the corrections
    deal with the 1-indexing of Julia =#
    return i,j
end

function square_to_linear_index(i::Int,j::Int,L::Int)::Int
    return L*(i-1) + j # 1 ≤ n ≤ L^2
end

is_in_bulk(i::Int,j::Int,L::Int)::Bool = !ontheborder(i,j,L)
function ontheborder(i::Int,j::Int,L::Int)::Bool
    if i == 1 || j == 1 || i == L || j == L
        return true
    else
        return false
    end
end

function at_least_at_distance_X_from_border(i,j,L;X=2)::Bool
    if i < X || j < X || i > L-X+1 || j > L-X+1
        return false
    else
        return true
    end
end

true_position(lattice::SquareLattice,i,j) = (i,j)
function true_position(lattice::TriangularLattice,i,j)
    if iseven(i) return (i,j+0.5)
    else return (i,j)
    end
end

function distance_to_border(thetas::Matrix{<:AbstractFloat},i,j)
    L = size(thetas,1)
    distance_to_left   = i-1
    distance_to_right  = L-i
    distance_to_bottom = j-1
    distance_to_top    = L-j
    return minimum([distance_to_top,distance_to_left,distance_to_right,distance_to_bottom])
end

function add_2_positions(pos1::Tuple{T,T},pos2::Tuple{T,T},L::T,should_take_mod::Bool)::Tuple{T,T} where T<:Int
    if should_take_mod return mod1.(pos1 .+ pos2 ,L)
    else return pos1 .+ pos2
    end
end

function mean_2_positions(pos1,pos2,L,should_take_mod::Bool=true)
    a,b = pos1 ; x,y = pos2

    dx = (x - a) #; dx = min(dx,L-dx)
    dy = (y - b) #; dy = min(dy,L-dy)

    if should_take_mod
        if abs(L-dx) < abs(dx) dx = -(L-dx) end
        # if L-dx < dx dx = -(L-dx) end
        if abs(L-dy) < abs(dy) dy = -(L-dy) end
        # if L-dy < dy dy = -(L-dy) end
        estimate_loc_collision = mod1.((a,b) .+ 0.5.*(dx,dy),L)
        # estimate_loc_collision = mod1.(0.5.*(a,b) .+ 0.5.*(x,y),L)
        return estimate_loc_collision
    end
    return (a,b) .+ 0.5.*(dx,dy) # if should_take_mod == false
end
# l = 100
# mean_2_positions((50,50),(60,60),l) == (55,55)
# mean_2_positions((10,10),(90,90),l) == (100,100)
# mean_2_positions((49,66),(51,61),l) == (50.0, 63.5)

function mean_N_positions(vec_pos,L,should_take_mod::Bool=true)
    averaged_pos = vec_pos[1]
    for i in 2:length(vec_pos)
        averaged_pos = mean_2_positions(averaged_pos,vec_pos[i],L,should_take_mod)
    end
    return averaged_pos
end

# coordinate-free methods based on scalar products (does not force SquareLattice)
function get_div_rot(thetas::Matrix{T},lattice::Abstract2DLattice) where T<:AbstractFloat # coordinate-free
    L = size(thetas,1)
    divergence = NaN*zeros(L,L)
    rotational = NaN*zeros(L,L)
    for j in 1:L , i in 1:L
        divergence[i,j] , rotational[i,j] = get_div_rot(thetas,i,j,lattice)
    end
    return divergence , rotational
end

function get_div_rot(thetas::Matrix{T},i,j,lattice::Abstract2DLattice) where T<:AbstractFloat # coordinate-free
    L = size(thetas,1)

    dummy_rho = one(T) # so that the NaN do not get filtered out, I need them here
    rho1_model = LangevinXY{T}(zero(T),"polar",zero(T),zero(T),dummy_rho)

    if isa(lattice,TriangularLattice) cst = π/3 # ; surface_unit_cell = 3sqrt(3)/2
    elseif isa(lattice,SquareLattice) cst = π/2 # ; surface_unit_cell = 2
    end

    angles_neighbours = get_neighbours(thetas,rho1_model,lattice,i,j,is_in_bulk(i,j,L))
    div = 0.0
    rot = 0.0
    for k in 1:length(angles_neighbours)
        if !isnan(angles_neighbours[k]) # normally it should not be, since preconditionning! first, but one never knows.
            div += cos((k-1)*cst) * cos(angles_neighbours[k]) + sin((k-1)*cst) * sin(angles_neighbours[k])
            rot += cos(k*cst) * cos(angles_neighbours[k]) + sin(k*cst) * sin(angles_neighbours[k])
        end
    end
    return div,rot
    #= No division by surface_unit_cell because eventually
    we do not care about the actual numerical value,
    only its sign. =#
end

# methods based on derivatives (forces SquareLattice)
# function get_div_rot(thetas::Matrix{T}) where T<:AbstractFloat # derivative based
#     L = size(thetas,1)
#     divergence = NaN*zeros(L,L)
#     rotational = NaN*zeros(L,L)
#     for j in 1:L , i in 1:L
#         divergence[i,j] , rotational[i,j] = get_div_rot(thetas,i,j)
#     end
#     return divergence , rotational
# end
#
# function get_div_rot(thetas::Matrix{T},i,j) where T<:AbstractFloat # derivative based
#     L = size(thetas,1)
#
#     dummy_rho = one(T) # so that the NaN do not get filtered out, I need them here
#     rho1_model = XY{T}(zero(T),"polar",zero(T),zero(T),dummy_rho)
#     square_lattice = SquareLattice(L,periodic=false)
#
#     angles_neighbours = get_neighbours(thetas,rho1_model,square_lattice,i,j,is_in_bulk(i,j,L))
#     # TODO check the order of angles_neighbours[1 or 3] and whether 2pi systematiquement
#     dx_theta = arclength(angles_neighbours[3],angles_neighbours[1],2pi)/2
#     dy_theta = arclength(angles_neighbours[4],angles_neighbours[2],2pi)/2
#     theta = thetas[i,j]
#     div = -sin(theta)*dx_theta + cos(theta)*dy_theta
#     rot =  cos(theta)*dx_theta + sin(theta)*dy_theta
#     return div,rot
# end

function zoom(thetas::Matrix{<:Number},lattice::Abstract2DLattice,i,j,window=WINDOW)
    #= Returns a 2window+1 portion of thetas (centered on i,j),
    together with a bool to say whether the operation has worked =#
    L = lattice.L
    i = round(Int,i)
    j = round(Int,j)

    if at_least_at_distance_X_from_border(i,j,L,X=window+1)
        thetas_zoom = thetas[i-window:i+window,j-window:j+window]
        no_problem_go_ahead = true
    else
        if lattice.periodic
            shift = 20
            thetas_shifted = copy(thetas)
            for i in 1:L thetas_shifted[i,:] = circshift(thetas_shifted[i,:],shift) end
            for i in 1:L thetas_shifted[:,i] = circshift(thetas_shifted[:,i],shift) end
            return zoom(thetas_shifted,lattice,mod1(i+shift,L),mod1(j+shift,L),window)

        else
            no_problem_go_ahead = false
            thetas_zoom = NaN
        end
    end
    return no_problem_go_ahead,thetas_zoom
end

## Rotations
rotate_clockwise90(thetas::Matrix{<:Number}) = rotr90(thetas .- pi/2)
rotate_counterclockwise90(thetas::Matrix{<:Number}) = rotl90(thetas .+ pi/2)
rotate_180(thetas::Matrix{<:Number}) = rot180(thetas .+ pi)
function randomly_rotate(thetas::Matrix{T})::Matrix{T} where T<:Number
    u = rand()
    if     u < 0.25  return thetas
    elseif u < 0.50  return rotate_clockwise90(thetas)
    elseif u < 0.75  return rotate_counterclockwise90(thetas)
    else             return rotate_180(thetas)
    end
end
