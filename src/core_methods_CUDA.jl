include("lattices.jl");
include("models.jl");
using StatsBase, Distributions, SpecialFunctions, CUDA


## ----------------- For NRXY (Non-Reciprocal XY Model, GPU version) ----------------- ##
## ----------------- For NRXY (Non-Reciprocal XY Model, GPU version) ----------------- ##
## ----------------- For NRXY (Non-Reciprocal XY Model, GPU version) ----------------- ##
## ----------------- For NRXY (Non-Reciprocal XY Model, GPU version) ----------------- ##

# for SquareLattice
function kernel_update_Square!(thetas_new, thetas, L::Int, R::Int, g, T::Tf, sigma::Tf, dt::Tf) where {Tf}
	i = (blockIdx().x-1)*blockDim().x + threadIdx().x
	j = (blockIdx().y-1)*blockDim().y + threadIdx().y
	k = (blockIdx().z-1)*blockDim().z + threadIdx().z

	check = i <= L && j <= L && k <= R # if L is not a power of 2, we need to check that we are not out of bounds
	if check
		# not slower than differentiating the bulk and the boundaries to avoid unnecessary mod1 usage
		ip = mod1(i+1,L)
		i_ = mod1(i-1,L)
		jp = mod1(j+1,L)
		jm = mod1(j-1,L)
		theta = thetas[i,j,k]

		force = sin(thetas[i,jp,k]-theta) * g(theta,0, sigma)    +    
				sin(thetas[ip,j,k]-theta) * g(theta,pi/2, sigma) + 
				sin(thetas[i,jm,k]-theta) * g(theta,-pi, sigma)  +  
				sin(thetas[i_,j,k]-theta) * g(theta,-pi/2, sigma) 

		thetas_new[i,j,k] = theta + force*dt/4 + sqrt(2T*dt) * randn(Tf)
	end
	
	return nothing
end

# for TriangularLattice
function kernel_update_Triangular!(thetas_new, thetas, L::Int, R::Int, g, T::Tf, sigma::Tf, dt::Tf) where {Tf}
	i = (blockIdx().x-1)*blockDim().x + threadIdx().x
	j = (blockIdx().y-1)*blockDim().y + threadIdx().y
	k = (blockIdx().z-1)*blockDim().z + threadIdx().z

	check = i <= L && j <= L && k <= R # if L is not a power of 2, we need to check that we are not out of bounds
	if check
        ip = mod1(i+1,L)
        i_ = mod1(i-1,L)
        jp = mod1(j+1,L)
        jm = mod1(j-1,L)
        theta = thetas[i,j,k]
        if iseven(i)
            force = sin(thetas[i,jp,k]-theta)  * g(theta,0, sigma)     +    
                    sin(thetas[i_,j,k]-theta)  * g(theta,pi/3, sigma)  + 
                    sin(thetas[i_,jm,k]-theta) * g(theta,2pi/3, sigma) +  
                    sin(thetas[i,jm,k]-theta)  * g(theta,pi, sigma)    +
                    sin(thetas[ip,jm,k]-theta) * g(theta,4pi/3, sigma) +
                    sin(thetas[ip,j,k]-theta)  * g(theta,5pi/3, sigma)

        else
            force = sin(thetas[i,jp,k]-theta)  * g(theta,0, sigma)     +    
                    sin(thetas[i_,jp,k]-theta) * g(theta,pi/3, sigma)  + 
                    sin(thetas[i_,j,k]-theta)  * g(theta,2pi/3, sigma) +  
                    sin(thetas[i,jm,k]-theta)  * g(theta,pi, sigma)    +
                    sin(thetas[ip,j,k]-theta)  * g(theta,4pi/3, sigma) +
                    sin(thetas[ip,jp,k]-theta) * g(theta,5pi/3, sigma)
        end
        thetas_new[i,j,k] = theta + force*dt/6 + sqrt(2T*dt) * randn(Tf)
    end
	
	return nothing
end

function update_Square!(thetas, thetas_new, L, R, g, T, sigma, dt)
	@cuda threads=block3D blocks=grid3D kernel_update_Square!(thetas_new, thetas, L, R, g, T, sigma, dt)
	@. thetas = thetas_new # quand tous les calculs sur la ligne sont ".able" , vraiment plus rapide que  "thetas .= thetas_new" dixit Ludovic Dumoulin
end


function update_Triangular!(thetas, thetas_new, L, R, g, T, sigma, dt)
	@cuda threads=block3D blocks=grid3D kernel_update_Triangular!(thetas_new, thetas, L, R, g, T, sigma, dt)
	@. thetas = thetas_new # quand tous les calculs sur la ligne sont ".able" , vraiment plus rapide que  "thetas .= thetas_new" dixit Ludovic Dumoulin
end

function evolve!(thetas, thetas_new, lattice, params, tmax)
    @unpack L, R, g, T, sigma, dt = params
    t = 0.0
    if lattice isa SquareLattice
        while t < tmax
            update_Square!(thetas, thetas_new, L, R, g, T, sigma, dt)
        t += dt
        end
    elseif lattice isa TriangularLattice
        while t < tmax
            update_Triangular!(thetas, thetas_new, L, R, g, T, sigma, dt)
        t += dt
        end
    end

    return thetas
end


## ----------------- For HydroXY (Hydrodynamic XY Model from Dan Pearce, GPU version) ----------------- ##
## ----------------- For HydroXY (Hydrodynamic XY Model from Dan Pearce, GPU version) ----------------- ##

# function kernel_update_HydroXY!(Pthetas_new, Pthetas,Pphis_new, Pphis, L::Int, R::Int, g, T::Tf, sigma::Tf, dt::Tf) where {Tf}
# 	i = (blockIdx().x-1)*blockDim().x + threadIdx().x
# 	j = (blockIdx().y-1)*blockDim().y + threadIdx().y
# 	k = (blockIdx().z-1)*blockDim().z + threadIdx().z
    
#     delta = 1/L

# 	check = i <= L && j <= L && k <= R # if L is not a power of 2, we need to check that we are not out of bounds
# 	if check
#         ip = mod1(i+1,L)
#         i_ = mod1(i-1,L)
#         jp = mod1(j+1,L)
#         jm = mod1(j-1,L)
#         Ptheta = thetas[i,j,k]
#         Pphi = phis[i,j,k]
#         d2Pphi_dphi2 = (phis[ip,j,k] - 2Pphi + phis[i_,j,k])/delta^2
#         d2Pphi_dtheta2 = (phis[i,jp,k] - 2Pphi + phis[i,jm,k])/delta^2
        
#         dPphi_dphi = (phis[ip,j,k] - phis[i_,j,k])/2/delta
#         dPphi_dtheta = (phis[i,jp,k] - phis[i,jm,k])/2/delta

#         Htheta = 
        
#         thetas_new[i,j,k] = theta + Htheta*dt + sqrt(2T*dt) * randn(Tf)
#     end
	
# 	return nothing
# end