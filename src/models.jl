using Parameters, Interpolations

abstract type AbstractModel{AbstractFloat} end
abstract type AbstractPropagationModel{AbstractFloat} <: AbstractModel{AbstractFloat} end

# function modd(model::AbstractModel{FT}) where {FT<:AbstractFloat}
#     if model.symmetry == "polar"
#         return FT(2pi)
#     elseif model.symmetry == "nematic"
#         return FT(pi)
#     end
# end
# function sym(model::AbstractModel{FT})::FT where {FT<:AbstractFloat}
#     if model.symmetry == "polar"
#         return 1.0
#     elseif model.symmetry == "nematic"
#         return 2.0
#     end
# end

## ---------------------------- Classical XY Model ---------------------------- ##
## ---------------------------- Classical XY Model ---------------------------- ##
## ---------------------------- Classical XY Model ---------------------------- ##
## ---------------------------- Classical XY Model ---------------------------- ##

abstract type AbstractXYModel{AbstractFloat} <: AbstractModel{AbstractFloat} end
function XY(params) # by default, return LangevinXY
    @unpack T, symmetry, dt, float_type, rho, algo = params
    if symmetry == "polar"
        if algo in ["MC", "MonteCarlo"]
            return MonteCarloXY{float_type}(T, symmetry, zero(float_type), rho)
        elseif algo == "Langevin"
            return LangevinXY{float_type}(T, symmetry, dt, zero(float_type), rho)
        else
            println()
            println("WARNING ! Unknown algo provided for XY model: return XY with Langevin dynamics.")
            return LangevinXY{float_type}(T, symmetry, dt, zero(float_type), rho)
        end
    elseif symmetry == "nematic"
        if algo in ["MC", "MonteCarlo"]
            return MonteCarloNematicXY{float_type}(T, symmetry, zero(float_type), rho)
        elseif algo == "Langevin"
            return LangevinNematicXY{float_type}(T, symmetry, dt, zero(float_type), rho)
        else
            println()
            println("WARNING ! Unknown algo provided for NematicXY model: return XY with Langevin dynamics.")
            return LangevinXY{float_type}(T, symmetry, dt, zero(float_type), rho)
        end
    else
        error("Symmetry not recognized. Please choose between 'polar' and 'nematic'.")
    end
end

mutable struct LangevinXY{AbstractFloat} <: AbstractXYModel{AbstractFloat}
    T::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat

end
function LangevinXY(params)
    @unpack T, symmetry, dt, float_type, rho = params
    T, dt, rho = convert.(float_type, (T, dt, rho))
    @assert symmetry == "polar" "LangevinXY only works for polar symmetry."
    return LangevinXY{float_type}(T, symmetry, dt, zero(float_type), rho)
end

mutable struct LangevinNematicXY{AbstractFloat} <: AbstractXYModel{AbstractFloat}
    T::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat

end
function LangevinNematicXY(params)
    @unpack T, symmetry, dt, float_type, rho = params
    T, dt, rho = convert.(float_type, (T, dt, rho))
    @assert symmetry == "nematic" "LangevinNematicXY only works for nematic symmetry."

    return LangevinNematicXY{float_type}(T, symmetry, dt, zero(float_type), rho)
end

# MonteCarlo
mutable struct MonteCarloXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    symmetry::String
    t::AbstractFloat
    rho::AbstractFloat

end
function MonteCarloXY(params)
    @unpack T, symmetry, float_type, rho = params
    T, rho = convert.(float_type, (T, rho))
    @assert symmetry == "polar" "MonteCarloXY only works for polar symmetry."

    return MonteCarloXY{float_type}(T, symmetry, zero(float_type), rho)
end

mutable struct MonteCarloNematicXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    symmetry::String
    t::AbstractFloat
    rho::AbstractFloat

end
function MonteCarloNematicXY(params)
    @unpack T, symmetry, float_type, rho = params
    T, rho = convert.(float_type, (T, rho))
    @assert symmetry == "nematic" "MonteCarloNematicXY only works for nematic symmetry."

    return MonteCarloNematicXY{float_type}(T, symmetry, zero(float_type), rho)
end


## ------------------------- Kuramoto Models ------------------------- ##
## ------------------------- Kuramoto Models ------------------------- ##
## ------------------------- Kuramoto Models ------------------------- ##
## ------------------------- Kuramoto Models ------------------------- ##

mutable struct ForcedXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    Var::AbstractFloat
    omegas_corr::AbstractFloat
    symmetry::String
    omegas::Matrix{AbstractFloat}
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end


function ForcedXY(params)
    @unpack L, T, Var, symmetry, dt, float_type, rho, omegas_corr = params
    T, Var, dt, rho = convert.(float_type, (T, Var, dt, rho, omegas_corr))

    pbc = true
    @unpack init = params
    if init == "single" pbc = false end
    
    if omegas_corr in [0,1]
        omegas = sqrt(Var) * randn(float_type, L, L)
    else
        @assert isinteger(L/omegas_corr) "L must be a multiple of omegas_corr."
        omegas_coarse = sqrt(Var) * randn(float_type, Int(L/omegas_corr)+1, Int(L/omegas_corr)+1) # uncorrelated coarse grid of omegas
        if pbc 
            omegas_coarse[1,:] = omegas_coarse[end,:] # enforce periodicity
            omegas_coarse[:,1] = omegas_coarse[:,end]
        end
        
        #= Interpolate omegas (Cubic Spline) on fine grid to create a smooth field 
        with a correlation length given by the coarse grid spacing "omegas_corr"=#
        itp = interpolate(omegas_coarse, BSpline(Cubic(Periodic(OnGrid()))))
        omegas = [itp(1+a/omegas_corr,1+b/omegas_corr) for a in 1:L, b in 1:L]
    end
    return ForcedXY{float_type}(T, Var, omegas_corr, symmetry, omegas, dt, zero(float_type), rho)
end

mutable struct OUForcedXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    Var::AbstractFloat
    omegas_corr::AbstractFloat # the difference with ForcedXY is that ω change in time, according to the Ornstein-Uhlenbeck process
    symmetry::String
    omegas::Matrix{AbstractFloat}
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end


function OUForcedXY(params)
    @unpack L, T, Var, symmetry, dt, float_type, rho, omegas_corr = params
    T, Var, dt, rho = convert.(float_type, (T, Var, dt, rho, omegas_corr))

    pbc = true
    @unpack init = params
    if init == "single" pbc = false end
    
    if omegas_corr in [0,1] # the field is uncorrelated
        omegas = sqrt(Var) * randn(float_type, L, L)
    else
        @assert isinteger(L/omegas_corr) "L must be a multiple of omegas_corr."
        omegas_coarse = sqrt(Var) * randn(float_type, Int(L/omegas_corr)+1, Int(L/omegas_corr)+1) # uncorrelated coarse grid of omegas
        if pbc 
            omegas_coarse[1,:] = omegas_coarse[end,:] # enforce periodicity
            omegas_coarse[:,1] = omegas_coarse[:,end]
        end
        
        #= Interpolate omegas (Cubic Spline) on fine grid to create a smooth field 
        with a correlation length given by the coarse grid spacing "omegas_corr"=#
        itp = interpolate(omegas_coarse, BSpline(Cubic(Periodic(OnGrid()))))
        omegas = [itp(1+a/omegas_corr,1+b/omegas_corr) for a in 1:L, b in 1:L]
    end
    return OUForcedXY{float_type}(T, Var, omegas_corr, symmetry, omegas, dt, zero(float_type), rho)
end

## --------------------- Self Propelled Particles (SPP) --------------------- ## 
## --------------------- Self Propelled Particles (SPP) --------------------- ## 
## --------------------- Self Propelled Particles (SPP) --------------------- ## 
## --------------------- Self Propelled Particles (SPP) --------------------- ## 

mutable struct SPP{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    A::AbstractFloat
    symmetry::String
    propulsion::String
    t::AbstractFloat
    rho::AbstractFloat
    algo::String
end
function SPP(params)
    @unpack T, A, rho, symmetry, algo, propulsion, float_type = params
    T, A, rho = convert.(float_type, (T, A, rho))

    return SPP{float_type}(T, A, symmetry, propulsion, zero(float_type), rho, algo)
end


## ------------------- Non Reciprocal XY Models ------------------- ##
## ------------------- Non Reciprocal XY Models ------------------- ##
## ------------------- Non Reciprocal XY Models ------------------- ##
## ------------------- Non Reciprocal XY Models ------------------- ##
# (Vision Cone)
mutable struct VisionConeXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    vision_cone::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function VisionConeXY(params)
    @unpack T, vision_cone, symmetry, dt, float_type, rho = params
    T, vision_cone, dt, rho = convert.(float_type, (T, vision_cone, dt, rho))

    if vision_cone ≠ 2π
        @assert symmetry == "polar" "VisionConeXY model only defined for polar symmetry"
    end

    return VisionConeXY{float_type}(T, vision_cone, symmetry, dt, zero(float_type), rho)
end

# Softly (but discretely, with projection on nearest neighbour) Tuned Couplings
mutable struct DiscreteNRXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    sigma::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function DiscreteNRXY(params)
    @unpack T, sigma, symmetry, dt, float_type, rho = params
    T, sigma, dt, rho = convert.(float_type, (T, sigma, dt, rho))

    return DiscreteNRXY{float_type}(T, sigma, symmetry, dt, zero(float_type), rho)
end

# Continous Couplings, without projection on nearest neighbour
mutable struct NRXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    sigma::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function NRXY(params)
    @unpack T, sigma, symmetry, dt, float_type, rho = params
    T, sigma, dt, rho = convert.(float_type, (T, sigma, dt, rho))
    @assert symmetry == "polar" "NRXY only works for polar symmetry."
    if dt < 5E-3
        println("WARNING ! dt = $dt is too small for NRXY model. dt = 1E-1 is used instead.")
        dt = 1E-1
    end
    return NRXY{float_type}(T, sigma, symmetry, dt, zero(float_type), rho)
end

mutable struct NRXY_ExtraTerm{AbstractFloat} <: AbstractModel{AbstractFloat}
    #= Same model as NRXY, but with the extra term in the dynamics 
    stemming from the energy defined in Loos et al. (2023)
    =#
    T::AbstractFloat
    sigma::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function NRXY_ExtraTerm(params)
    @unpack T, sigma, symmetry, dt, float_type, rho = params
    T, sigma, dt, rho = convert.(float_type, (T, sigma, dt, rho))
    @assert symmetry == "polar" "NRXY_ExtraTerm only works for polar symmetry."
    if dt < 5E-3
        println("WARNING ! dt = $dt is too small for NRXY_ExtraTerm model. dt = 1E-1 is used instead.")
        dt = 1E-1
    end
    return NRXY_ExtraTerm{float_type}(T, sigma, symmetry, dt, zero(float_type), rho)
end

mutable struct ContinuousNRNematicXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    sigma::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function ContinuousNRNematicXY(params)
    @unpack T, sigma, symmetry, dt, float_type, rho = params
    T, sigma, dt, rho = convert.(float_type, (T, sigma, dt, rho))
    @assert symmetry == "nematic" "ContinuousNRNematicXY only works for nematic symmetry."

    return ContinuousNRNematicXY{float_type}(T, sigma, symmetry, dt, zero(float_type), rho)
end

mutable struct HydrodynNRXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    sigma::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function HydrodynNRXY(params)
    @unpack T, sigma, symmetry, dt, float_type, rho, alpha = params
    T, sigma, dt, rho, alpha = convert.(float_type, (T, sigma, dt, rho, alpha))

    return HydrodynNRXY{float_type}(T, sigma, symmetry, dt, zero(float_type), rho, alpha)
end

## ------------------- Propagation Models ------------------- ##
## ------------------- Propagation Models ------------------- ##
## ------------------- Propagation Models ------------------- ##
## ------------------- Propagation Models ------------------- ##

mutable struct PropagationForcedXY{AbstractFloat} <: AbstractPropagationModel{AbstractFloat}
    T::AbstractFloat
    Var::AbstractFloat
    symmetry::String
    omegas::AbstractArray{AbstractFloat}
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function PropagationForcedXY(params)
    @unpack L, T, Var, symmetry, dt, float_type, rho = params
    T, Var, dt, rho = convert.(float_type, (T, Var, dt, rho))

    omegas = zeros(float_type, L) # dummy
    # will be instanciated later during the init, when one knows the type of lattice (1D, 2D)

    rho1 = one(float_type) # for now, we don't want to investigate holes in this model
    return PropagationForcedXY{float_type}(T, Var, symmetry, omegas, dt, zero(float_type), rho1)
end


## ------------------- Hydrodynamic Vector Fields Models ------------------- ##
mutable struct HydroVector{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    Var::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
    alpha::AbstractFloat
end
function HydroVector(params)
    @unpack L, T, Var, symmetry, dt, float_type, rho, alpha = params
    T, Var, dt, rho, alpha = convert.(float_type, (T, Var, dt, rho, alpha))

    rho1 = one(float_type) # for now, we don't want to investigate holes in this model
    return HydroVector{float_type}(T, Var, symmetry, dt, zero(float_type), rho1, alpha)
end

mutable struct HydroTest{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    Var::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
    alpha::AbstractFloat
    sigma::AbstractFloat
end
function HydroTest(params)
    @unpack L, T, Var, symmetry, dt, float_type, rho, alpha, sigma = params
    T, Var, dt, rho, alpha = convert.(float_type, (T, Var, dt, rho, alpha, sigma))

    rho1 = one(float_type) # for now, we don't want to investigate holes in this model
    return HydroTest{float_type}(T, Var, symmetry, dt, zero(float_type), rho1, alpha, sigma)
end 