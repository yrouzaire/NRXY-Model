include("lattices.jl");
include("models.jl");
include("core_methods.jl");

using Hungarian
import Random.shuffle

## ----------------- Spotting Defect : core methods ----------------- ##
## ----------------- Spotting Defect : core methods ----------------- ##
## ----------------- Spotting Defect : core methods ----------------- ##
## ----------------- Spotting Defect : core methods ----------------- ##
## ----------------- Spotting Defect : core methods ----------------- ##

"""
    arclength(theta1::T,theta2::T,symm)::T where T<:AbstractFloat

-----------------------------
This function returns the signed arclength (in radians) on the unit trigonometric circle .
    Clockwise        >> sign -
    Counterclockwise >> sign +
    WARNING : Note that the inputs thetas need to lie within [0,2π] or [0,π], depending on the symmetry `symm` (resp. polar or nematic) of the model
-----------------------------
"""
function arclength(theta1::T,theta2::T,symm)::T where T<:AbstractFloat
    #= This function returns the signed arclength (in radians) on the unit trigonometric circle .
    Clockwise        >> sign -
    Counterclockwise >> sign +
    WARNING : Note that the inputs thetas need to lie within [0,2π] or [0,π], depending on the symmetry `symm` (resp. polar or nematic) of the model =#
    dtheta = theta2 - theta1
    dtheta_abs = abs(theta2 - theta1)

    shortest_unsigned_arclength = min(symm-dtheta_abs,dtheta_abs)
    if dtheta_abs ≤ symm/2
        signe = sign(dtheta)
    else
        signe = -sign(dtheta)
    end
    return signe*shortest_unsigned_arclength
end


"""
    spot_defects(thetas::Matrix{T}, model::AbstractModel, lattice::Abstract2DLattice) where T<:AbstractFloat

Returns two arrays of tuples (x,y,q,µ) where 
    - x,y are the coordinates of the defect's core (they live on the dual lattice so they are real number)
    - q is the topological charge of the defect (± 1, ± 1/2 ) 
    - µ is the shape of the defect, i.e. the angle of the director field at the right of the defect's core. It is defined from 0 to 2π (modulo 2π).

First output array : + defects, second output array : - defects
Use : defects_plus,defects_minus = spot_defects(thetas,model,lattice).

This function first performs a preconditionning of the theta field:
    - if needed (because rho < 1 for instance), it fills the holes.
    - it relaxes the field at T = 0 for Δt = 2 (200 timesteps)

It then scan all the plaquettes of the lattice and spot the defects, cf. `xyqµ_plaquette` (@xyqµ_plaquette)
"""
function spot_defects(thetas::Matrix{T}, model::AbstractModel, lattice::Union{Abstract2DLattice,TorusLattice}) where T<:AbstractFloat
    L = lattice.L
    defects_plus  = Tuple{T,T,T,T}[] # (x,y,q,µ) where q is the charge and µ the shape of the defect
    defects_minus = Tuple{T,T,T,T}[]

    ## Preconditioning
    coeff_mod = model.symmetry == "nematic" ? pi : 2pi
    thetasmod = mod.(copy(thetas),coeff_mod)
    if model.rho < 1 preconditionning!(thetasmod,model,lattice) end
    trelax = 0 ; Trelax = 0.
    # thetasmod = relax(thetasmod,model,lattice,trelax=trelax,T=Trelax)


    ## Spotting the defects
    if lattice isa TorusLattice || lattice.periodic range_bc = 1:L else range_bc = 3:L-2 end # 3:L-2 instead of 2:L-1 because for TriangularLattices, we need to take (i+2,j)
    for i in range_bc
        for j in range_bc
            xyqµ_defects = xyqµ_plaquette(thetasmod,model,lattice,i,j,coeff_mod)
            for (x,y,q,µ) in xyqµ_defects
                if     q > + 0.1 push!(defects_plus,(x,y,q,µ))
                elseif q < - 0.1 push!(defects_minus,(x,y,q,µ))
                end
            end
        end
    end

    return defects_plus,defects_minus
end

"""
    xyqµ_plaquette(thetas,model,lattice::SquareLattice,i,j,symm_coeff)

Preambule : Since the procedure for `SquareLattice` and `TriangularLattice` differ a bit (cf. next function), 
it is clearer to define two distinct functions.

Description of the function xyqµ_plaquette for `SquareLattice`.
-----------------------------
This function returns the list of the defects (x,y,q,µ) in the plaquette (i,j) of the lattice.
The procedure reads as follows:
1. Define the plaquette. 
    - SquareLattice => square dual plaquette 
    - TriangularLattice => double triangle dual plaquette (cf. next function)
2. Compute the charge q
3. If q ≠ 0, compute the shape µ (otherwise return a dummy µ)
3*. Small corrections to µ (of the order of ± 0.05)
4. Compute the position of the defect's core (x,y)
-----------------------------
"""
function xyqµ_plaquette(thetas::Matrix{T},model,lattice::Union{SquareLattice,TorusLattice},i::Int,j::Int,symm_coeff) where T<:AbstractFloat
    # 1. Define the plaquette 
    angles_plaquette_for_q,angles_plaquette_for_µ = get_angles_plaquette(thetas,lattice,i,j)
    
    # 2. Compute the charge q
    q = 0.0
    nnq = length(angles_plaquette_for_q)
    for i in 1:nnq
        q += arclength(angles_plaquette_for_q[i],angles_plaquette_for_q[mod1(i+1,nnq)],symm_coeff)
    end
    q = round(q/2π,digits=1) # to avoid numerical errors such as 0.99999 ≠ 1
    
    # 3. If q ≠ 0, compute the shape µ (otherwise return a dummy µ)
    µ = infer_mu(angles_plaquette_for_µ,lattice,q) # corrections included
    
    # 4. Corrections x,y due to the plaquette choice
    x,y = (i+0.5,j+0.5)
    
    return [(x,y,q,µ)]
end

"""
    xyqµ_plaquette(thetas,model,lattice::TriangularLattice,i,j,symm_coeff)

Preambule : Since the procedure for `SquareLattice` and `TriangularLattice` differ a bit (cf. next function), 
it is clearer to define two distinct functions.

Description of the function xyqµ_plaquette for `TriangularLattice`.
-----------------------------

The procedure is very similar to the SquareLattice one.
However it has its own specificities. Indeed, while the dual of a 
square lattice is also square and thus allows for a straightforward
tessalation of the space, it's not the case anymore for the triangular 
lattice, which admits a hexagonal dual lattice. For each lattice point (i,j), 
I then have to define two dual plaquettes, one with the (i,j) point at its top 
(up triangle), and one with the (i,j) point at its bottom (down triangle). 

The only difference between the method for TriangularLattice and SquareLattice 
is that the present one (for TriangularLattice), returns the result x,y,q,µ for 
two plaquettes instead of only one for SquareLattice.
"""
function xyqµ_plaquette(thetas,model,lattice::TriangularLattice,i,j,symm_coeff)
    #= The procedure is very similar to the SquareLattice one.
    However it has its own specificities. Indeed, while the dual of a 
    square lattice is also square and thus allows for a straightforward
    tessalation of the space, it's not the case anymore for the triangular 
    lattice, which admits a hexagonal dual lattice. For each lattice point (i,j), 
    I then have to define two dual plaquettes, one with the (i,j) point at its top 
    (up triangle), and one with the (i,j) point at its bottom (down triangle). 
    =#
    xyqµs = []

# First plaquette, "up triangle" 🔺 
    type_dual_plaquette = "up"
    angles_plaquette_for_q,angles_plaquette_for_µ = get_angles_plaquette(thetas,lattice,i,j,type_dual_plaquette)

    q = 0.0
    nnq = length(angles_plaquette_for_q)
    for i in 1:nnq 
        q += arclength(angles_plaquette_for_q[i],angles_plaquette_for_q[mod1(i+1,nnq)],symm_coeff)
    end
    q = round(q/2π,digits=1) # to avoid numerical errors such as 0.99999 ≠ 1

    µ   = infer_mu(angles_plaquette_for_µ,lattice,q,type_dual_plaquette) # corrections included
    x,y = (i+0.5,j)
    push!(xyqµs,(x,y,q,µ))
    
# Second plaquette, "down triangle"🔻
    type_dual_plaquette = "down"
    angles_plaquette_for_q,angles_plaquette_for_µ = get_angles_plaquette(thetas,lattice,i,j,type_dual_plaquette)
    
    q = 0.0
    nnq = length(angles_plaquette_for_q)
    for i in 1:nnq 
        q += arclength(angles_plaquette_for_q[i],angles_plaquette_for_q[mod1(i+1,nnq)],symm_coeff)
    end
    q = round(q/2π,digits=1) # to avoid numerical errors such as 0.99999 ≠ 1

    µ   = infer_mu(angles_plaquette_for_µ,lattice,q,type_dual_plaquette) # corrections included
    x,y = (i+0.5*cos(π/6),j+0.5)
    push!(xyqµs,(x,y,q,µ))
    
    return xyqµs
end

function infer_mu(angles_plaquette_for_µ,lattice,q,type_dual_plaquette="")
    if abs(q) < 0.1 
        return 0 
    else
        #= Arclengths from the horizontal axis to the corner of the plaquette, 
        independently of the value of the angle at the corner =#
        if lattice isa SquareLattice || lattice isa TorusLattice        
            arclengths = [3π/4,-3π/4,-π/4,π/4]
            correction = +0.01
        elseif lattice isa TriangularLattice 
            if type_dual_plaquette == "up"
                arclengths = [π/2,5π/6,-5π/6,-π/2,-π/6,π/6] # first and second NN
                arclengths = [π/2,-5π/6,-π/6] # first NN
                correction = +0.01
            elseif type_dual_plaquette == "down"
                arclengths = [5π/6,-5π/6,-π/2,-π/6,π/6,π/2] # first and second NN
                arclengths = [5π/6,-π/2,π/6] # first NN
                correction = +0.01
            end
        else error("Lattice unknown. Abort mission and come back to Earth.")
        end

        mu_hat = angles_plaquette_for_µ - q*arclengths 
        µ = angle(sum(exp.(im*mu_hat)))
        µ = mod(µ+correction,2π)   
        return µ 
    end
end

function get_angles_plaquette(thetas::Matrix{T},lattice::Union{SquareLattice,TorusLattice},i::Int,j::Int) where T<:AbstractFloat
    #= For a square lattice, I define the following arbirary convention for the plaquette :
    1:(i,j)    4:(i,j+1)
    2:(i+1,j)  3:(i+1,j+1)
    =#
    L = lattice.L
    in_bulk = is_in_bulk(i,j,L)

    if in_bulk
        angles_plaquette_for_q = [thetas[i,j], thetas[i+1,j], thetas[i+1,j+1], thetas[i,j+1]]
    else 
        angles_plaquette_for_q = [thetas[i,j], thetas[mod1(i+1,L),j], thetas[mod1(i+1,L),mod1(j+1,L)], thetas[i,mod1(j+1,L)]]
    end

    angles_plaquette_for_µ = angles_plaquette_for_q
    return angles_plaquette_for_q,angles_plaquette_for_µ
end

function get_angles_plaquette(thetas,lattice::TriangularLattice,i,j,type_dual_plaquette::String)
    #= Recall, for a triangular lattice, one has:
        If i is even :
            1:(i,j)  
        2:(i+1,j)  3:(i+1,j+1)
            
        If i is odd :
            1:(i,j)  
        2:(i+1,j-1)  3:(i+1,j)

        ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ 

        ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ 

        For type_dual_plaquette == "up", this translates into

        If i is even :
        2:(i,j-1)    1:(i,j)   6:(i,j+1)
            3:(i+1,j)   5:(i+1,j+1)
                    4:(i+2,j) 
        If i is odd 
        2:(i,j-1)    1:(i,j)   6:(i,j+1)
            3:(i+1,j-1)   5:(i+1,j)
                    4:(i+2,j)
            
        For type_dual_plaquette == "down", this translates into
        If i (of point 1)  is even : 
                    6:(i-1,j+1)
                1:(i,j)     5:(i,j+1)     
          2:(i+1,j) 3:(i+1,j+1) 4:(i+1,j+2)
        
        If i (of point 1) is even : 
                    6:(i-1,j)
                1:(i,j)     5:(i,j+1)     
          2:(i+1,j-1) 3:(i+1,j) 4:(i+1,j+2)
    =#

    L = lattice.L
    in_bulk = at_least_at_distance_X_from_border(i,j,L,X=3)
    if type_dual_plaquette == "up"
        if in_bulk
            if iseven(i) angles_plaquette_for_µ = [thetas[i,j],thetas[i,j-1],thetas[i+1,j],thetas[i+2,j],thetas[i+1,j+1],thetas[i,j+1]]
            else         angles_plaquette_for_µ = [thetas[i,j],thetas[i,j-1],thetas[i+1,j-1],thetas[i+2,j],thetas[i+1,j],thetas[i,j+1]]
            end
        else
            if iseven(i) angles_plaquette_for_µ = [thetas[i,j],thetas[mod1(i,L),mod1(j-1,L)],thetas[mod1(i+1,L),j],thetas[mod1(i+2,L),j],thetas[mod1(i+1,L),mod1(j+1,L)],thetas[mod1(i,L),mod1(j+1,L)]]
            else         angles_plaquette_for_µ = [thetas[i,j],thetas[mod1(i,L),mod1(j-1,L)],thetas[mod1(i+1,L),mod1(j-1,L)],thetas[mod1(i+2,L),j],thetas[mod1(i+1,L),j],thetas[mod1(i,L),mod1(j+1,L)]]
            end 
        end
    elseif type_dual_plaquette == "down"
        if in_bulk
            if iseven(i) angles_plaquette_for_µ = [thetas[i,j],thetas[i+1,j],thetas[i+1,j+1],thetas[i+1,j+2],thetas[i,j+1],thetas[i-1,j+1]]
            else         angles_plaquette_for_µ = [thetas[i,j],thetas[i+1,j-1],thetas[i+1,j],thetas[i+1,j+2],thetas[i,j+1],thetas[i-1,j]]
            end        
        else
            if iseven(i) angles_plaquette_for_µ = [thetas[i,j],thetas[mod1(i+1,L),j],thetas[mod1(i+1,L),mod1(j+1,L)],thetas[mod1(i+1,L),mod1(j+2,L)],thetas[i,mod1(j+1,L)],thetas[mod1(i-1,L),mod1(j+1,L)]]
            else         angles_plaquette_for_µ = [thetas[i,j],thetas[mod1(i+1,L),mod1(j-1,L)],thetas[mod1(i+1,L),j],thetas[mod1(i+1,L),mod1(j+2,L)],thetas[i,mod1(j+1,L)],thetas[mod1(i-1,L),j]]
            end
        end
    end

    angles_plaquette_for_q = angles_plaquette_for_µ[[1,3,5]] # the angles used to compute the charge are the ones of the inner triangle
    return angles_plaquette_for_q,angles_plaquette_for_q
end


# function alone_in_window(pos,pos_all,lattice,window)::Bool
#     alone_in_window = true
#     for loc in pos_all
#         distance = dist(lattice,pos,loc)
#         if 0 < distance ≤ window+1 # take margin because non integer positions might get rounded
#             alone_in_window = false
#             break
#         end
#     end
#     return alone_in_window
# end

## ----------------- Spotting Defects : preconditionning and auxiliary methods ----------------- ##
## ----------------- Spotting Defects : preconditionning and auxiliary methods ----------------- ##
## ----------------- Spotting Defects : preconditionning and auxiliary methods ----------------- ##
## ----------------- Spotting Defects : preconditionning and auxiliary methods ----------------- ##
## ----------------- Spotting Defects : preconditionning and auxiliary methods ----------------- ##

function preconditionning!(thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::Abstract2DLattice;trelax = 0.5, Trelax=0)
    remove_isolates!(thetas,model,lattice)
    fill_holes!(thetas,model,lattice)
    thetas = relax(thetas,model,lattice,trelax=trelax,T=Trelax)
end

function remove_isolates!(thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::Abstract2DLattice)
    L = size(thetas,1)
    for j in 1:L
        for i in 1:L
            if isempty(get_neighbours(thetas,model,lattice,i,j,is_in_bulk(i,j,L)))
                thetas[i,j] = NaN
            end
        end
    end
    return thetas
end

function fill_holes!(thetas::Matrix{T},model::AbstractModel{T},lattice::Abstract2DLattice) where T<:AbstractFloat
    L = size(thetas,1)
    # model.symmetry == "polar" ? symm = 2 : symm = 1
    while count(isnan,thetas) > 0
        for j in shuffle(1:L) , i in shuffle(1:L)
            if isnan(thetas[i,j])
                angle_neighbours = invoke(get_neighbours,Tuple{Matrix{T},AbstractModel{T},typeof(lattice),Int,Int,Bool},  thetas,model,lattice,i,j,is_in_bulk(i,j,L))
                # thetas[i,j] = angle(sum(exp,symm*im*angle_neighbours, init=0))/symm # moyenne des voisins
                if !isempty(angle_neighbours) thetas[i,j] = rand(angle_neighbours) end
            end
        end
    end
    return thetas
end

## ----------------- Spotting Defects : derivative methods ----------------- ##
## ----------------- Spotting Defects : derivative methods ----------------- ##
## ----------------- Spotting Defects : derivative methods ----------------- ##
## ----------------- Spotting Defects : derivative methods ----------------- ##
## ----------------- Spotting Defects : derivative methods ----------------- ##

function number_defects(thetas,model,lattice) 
    df = spot_defects(thetas,model,lattice)
    return length(df[1]) + length(df[2])
end

## ---------------------- Tracking Defects : dedicated Data Structures  ---------------------- ##
## ---------------------- Tracking Defects : dedicated Data Structures  ---------------------- ##
## ---------------------- Tracking Defects : dedicated Data Structures  ---------------------- ##
## ---------------------- Tracking Defects : dedicated Data Structures  ---------------------- ##
## ---------------------- Tracking Defects : dedicated Data Structures  ---------------------- ##
mutable struct Defect
    id::Int
    charge::Number
    shape::Vector{Float64} # shapes might change over the simulation
    pos::Vector{Tuple{Number,Number}}
    annihilation_time::Union{Float64,Nothing}
    creation_time::Float64
    id_annihilator::Union{Int,Nothing}
end
Defect(;id,charge,loc,t,shape=NaN) = Defect(id,charge,[shape],[loc],nothing,t,nothing)

first_loc(d::Defect) = d.pos[1]
last_loc(d::Defect)  = d.pos[end]

last_shape(d::Defect)  = d.shape[end]
first_shape(d::Defect) = d.shape[1]
shape(d::Defect)       = d.shape[end]


function compute_velocity_defect(df::Defect, times;over=1)
	pos = df.pos
	vel = zeros(length(pos)-over)
	for i in 1:length(pos)-over
		tmp = (pos[i+over] .- pos[i]) ./(times[i+over] - times[i])
        # doesnt take PBC into account mais pas important pour ce que j'en fais
		vel[i] = norm(tmp)
	end
	return vel
end


function compute_acceleration_defect(df::Defect, times;over=3)
	pos = df.pos
	over2 = round(Int,(over-1)/2)
	acc = zeros(length(pos)-2over2)
	for i in over2+1:length(pos)-over2-1
		tmp = (pos[i+over2] .- 2.0 .* pos[i] .+ pos[i-over2]) ./ (times[i+over2] - times[i-over2])^2
        # doesnt take PBC into account mais pas important pour ce que j'en fais
		acc[i] = norm(tmp)
	end
	return acc
end

function min_abs(x)
    x_abs = abs.(x)
    ind = -1
    tmp = Inf
    for i in eachindex(x)
        if x_abs[i] < tmp
            ind = i
            tmp  = x_abs[i]
        end
    end
    return x[ind]
end

function acceleration(d::Defect,lattice::Abstract2DLattice;over=1)::Tuple{Number,Number}
    # the velocity is symmetrically computed over "over" timesteps
    if length(d.pos) ≥ 2over + 1
        a,b = d.pos[end-over] ; x,y = d.pos[end]  ; f,g = d.pos[end-2over]
        d2x = x+f-2a
        d2y = y+g-2b
        if lattice.periodic
            d2x = min_abs([d2x,lattice.L+d2x,d2x-lattice.L])
            d2y = min_abs([d2y,lattice.L+d2y,d2y-lattice.L])
        end
        return (d2x,d2y)
    else
        println("Not enough information, the history of position has to be of length ≥ 3 to infer a second derivative.")
        return (NaN,NaN)
    end
end

function update_position_and_shape!(d::Defect,new_loc,new_shape)
    push!(d.pos,new_loc)
    if new_shape == "unknown"  push!(d.shape,last_shape(d)) # by default, if unknown, push the last known shape
    else push!(d.shape,new_shape)
    end
end

mutable struct DefectTracker
    defectsP::Vector{Defect} # the id of a defect is its index in this vector
    defectsN::Vector{Defect} # so there is a (+)defect with id=1 AND and a (-)defect with id=1
    current_time::Float64 # latest update time (by convention, the creation time of the whole data structure = 0)

    function DefectTracker(thetas,model,lattice) # constructor
        vortices,antivortices = spot_defects(thetas,model,lattice)
        defectsP = [Defect(id=i,charge=vortices[i][3],shape=vortices[i][4],loc=vortices[i][1:2],t=model.t) for i in eachindex(vortices)]
        defectsN = [Defect(id=i,charge=antivortices[i][3],shape=antivortices[i][4],loc=antivortices[i][1:2],t=model.t) for i in eachindex(antivortices)]
        new(defectsP,defectsN,model.t)
    end
end
number_defectsP(dft::DefectTracker) = length(dft.defectsP)
number_defectsN(dft::DefectTracker) = length(dft.defectsN)
number_defects(dft::DefectTracker)  = length(dft.defectsN)  + length(dft.defectsN)

number_active_defectsP(dft::DefectTracker) = count(isnothing,[d.annihilation_time for d in dft.defectsP])
number_active_defectsN(dft::DefectTracker) = count(isnothing,[d.annihilation_time for d in dft.defectsN])
number_active_defects(dft::DefectTracker)  = number_active_defectsN(dft) + number_active_defectsP(dft)

active_defectsP(dft) = [d for d in dft.defectsP if isnothing(d.annihilation_time)]
active_defectsN(dft) = [d for d in dft.defectsN if isnothing(d.annihilation_time)]
active_defects(dft)  = vcat(active_defectsP(dft),active_defectsN(dft))

get_mus(dft::DefectTracker) = [d.shape[end] for d in active_defects(dft)]
get_musP(dft::DefectTracker) = [d.shape[end] for d in active_defectsP(dft)]
get_musN(dft::DefectTracker) = [d.shape[end] for d in active_defectsN(dft)]


function coordinates_defects(dft::DefectTracker;separate_xy = true)
    coordinate_positive_defects = Vector{Tuple{Number,Number}}(undef, number_active_defectsP(dft))
    coordinate_negative_defects = Vector{Tuple{Number,Number}}(undef, number_active_defectsN(dft))
    pos_defs = active_defectsP(dft)
    for i in each(pos_defs)
        coordinate_positive_defects[i] = pos_defs[i].pos[1]
    end
    pos_defs = active_defectsN(dft)
    for i in each(pos_defs)
        coordinate_negative_defects[i] = pos_defs[i].pos[1]
    end
    if separate_xy
        x_positive_defects = [x for (x,y) in coordinate_positive_defects]
        y_positive_defects = [y for (x,y) in coordinate_positive_defects]
        x_negative_defects = [x for (x,y) in coordinate_negative_defects]
        y_negative_defects = [y for (x,y) in coordinate_negative_defects]
        return x_positive_defects,y_positive_defects,x_negative_defects,y_negative_defects
    else
        return coordinate_positive_defects,coordinate_negative_defects
    end
end


function t_bounds(dft::DefectTracker)
    alldefects = vcat(dft.defectsP,dft.defectsN)
    creation_times     = filter(!isnothing,[d.creation_time for d in alldefects])
    annihilation_times = filter(!isnothing,[d.annihilation_time for d in alldefects])
    if isempty(creation_times) tmin = 0.0
    else tmin = minimum(creation_times)
    end
    if isempty(annihilation_times) tmax = dft.current_time
    else tmax = maximum(creation_times)
    end
    return tmin, tmax
end

function ID_active_defects(dt::DefectTracker)
    activeP = Int[]
    for i in 1:number_defectsP(dt)
        if dt.defectsP[i].annihilation_time == nothing push!(activeP,i) end
    end
    activeN = Int[]
    for i in 1:number_defectsN(dt)
        if dt.defectsN[i].annihilation_time == nothing push!(activeN,i) end
    end
    return activeP,activeN
end

function add_defect!(dt::DefectTracker;charge,loc,shape=NaN)
    if charge > 0 push!(dt.defectsP,Defect(id=1+number_defectsP(dt),charge=charge,shape=shape,loc=loc,t=dt.current_time))
    else          push!(dt.defectsN,Defect(id=1+number_defectsN(dt),charge=charge,shape=shape,loc=loc,t=dt.current_time))
    end
end

function pair_up_hungarian(dt::DefectTracker,new,old,lattice::Union{TorusLattice,Abstract2DLattice},charge::String)
    # charge can be either "+" or "-"
    distance_matrixx = distance_matrix(new,old,lattice) # m_new lignes, m_old colonnes
    proposal         = hungarian(distance_matrixx)[1] # length = length(new)
    assignment       = copy(proposal) # because it will be modified in the next for loop

    #= A few comments :
    1. The proposal is the match between indices of the vectors
    new,old while the assignment matches actual IDs of the DefectTracker.
    2. If cost_matrix is a NxM matrix (workers x jobs), the output of hungarian(cost_matrix)
    is a Nx1 vector containing the assignments of each of the N workers to the indices of the jobs (1 < indices < M).
    A "0" means the worker has not been assigned to any job, hence the assignment[i] ≠ 0 condition below.
    =#
    if charge == "+"
        for i in eachindex(assignment)
            for j in 1:number_defectsP(dt)
                if assignment[i] ≠ 0 && dt.defectsP[j].annihilation_time == nothing && last_loc(dt.defectsP[j]) == old[proposal[i]]
                    # the first condition is a protection against the creation case, where the Hungarian algo matches the newly created vortex to 0
                    # the second condition ensures that one only considers currently living vortices and not vortices now annihilated
                    assignment[i] = j
                    break # breaks innerloop only
                end
            end
        end
    else # if charge  == "-"
        for i in eachindex(assignment)
            for j in 1:number_defectsN(dt)
                if assignment[i] ≠ 0 && dt.defectsN[j].annihilation_time == nothing && last_loc(dt.defectsN[j]) == old[proposal[i]]
                    assignment[i] = j
                    break
                end
            end
        end
    end
    return assignment
end

function find_closest_before_annihilation(dt,lattice,old_loc_defect)
    distance = Inf ; ID_antidefect = -1 # dummy
    # println(length((dt.defectsN)))
    for i in eachindex(dt.defectsN)
        if dt.defectsN[i].annihilation_time == dt.current_time # it has just annihilated
            tmp = dist(lattice,old_loc_defect,last_loc(dt.defectsN[i]))
            if tmp < distance
                distance = tmp
                ID_antidefect = i
                # println("Found annihilating defect: $(ID_antidefect)")
            end
        end
    end
    if ID_antidefect == -1
        @warn "Annihilating defect not found: annihilation location will not be computed."
        println([dt.defectsN[i].annihilation_time == dt.current_time for i in 1:number_defectsN(dt)])
        return -1,(-1,-1) # dummy
    else
        return ID_antidefect,last_loc(dt.defectsN[ID_antidefect])
    end
end

function annihilate_defects(dt::DefectTracker,ids_annihilated_defects,lattice)
    for i in ids_annihilated_defects
        old_loc_vortex = last_loc(dt.defectsP[i])
        ID_antivortex,old_loc_antivortex = find_closest_before_annihilation(dt,lattice,old_loc_vortex)
        if ID_antivortex ≠ -1
            dt.defectsP[i].id_annihilator = ID_antivortex
            dt.defectsN[ID_antivortex].id_annihilator = i

            # dt.defectsP[i].annihilation_time = dt.current_time
            # dt.defectsN[ID_antivortex].annihilation_time = dt.current_time

            estimate = mean_2_positions(old_loc_vortex,old_loc_antivortex,lattice.L)
            update_position_and_shape!(dt.defectsP[i],estimate,last_shape(dt.defectsP[i]))
            update_position_and_shape!(dt.defectsN[ID_antivortex],estimate,last_shape(dt.defectsN[ID_antivortex]))
        end
    end
    return dt
end

function update_and_track!(thetas::Matrix{T},model::AbstractModel{T},lattice::Union{TorusLattice,Abstract2DLattice},dft::DefectTracker,times;verbose=false, plot=false, defects=false,stop_no_defects=false) where T<:AbstractFloat
    if times[1] == 0 && number_active_defects(dft) > 0 # the initialisation has already been done
        times = times[2:end] # we skip the first time so that the lenghts of times and dft internal arrays (positions of the defects for instance) match
    end
    for tt in each(times)
        evolve!(thetas,model,lattice,times[tt])
        if verbose println("t = ",round(times[tt],digits=1)," n = ",number_active_defects(dft)) end
        
        # try 
            update_DefectTracker!(dft,thetas,model,lattice)
            if plot
                display(plot_thetas(thetas,model,lattice,defects=defects))
            end    
            if stop_no_defects && number_active_defects(dft) == 0 
                println("No defects left, simulation stopped.")
                break
            end
        # catch e
        #     println(e)
        #     println("Uncaught error, simulation aborted.")
        #     break
        # end
    end
    return dft
end


function update_DefectTracker!(dt::DefectTracker,thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::Union{TorusLattice,Abstract2DLattice})
    dt.current_time = model.t
    # NB = lattice.periodic
    vortices_new,antivortices_new = spot_defects(thetas,model,lattice)

    # if BC == "periodic" @assert length(vortices_new) == length(antivortices_new) && length(vortices_old) == length(antivortices_old) end
    locP_old    = [last_loc(dt.defectsP[i]) for i in eachindex(dt.defectsP)]
    locN_old    = [last_loc(dt.defectsN[i]) for i in eachindex(dt.defectsN)]
    chargeP_old = [dt.defectsP[i].charge    for i in eachindex(dt.defectsP)]
    chargeN_old = [dt.defectsN[i].charge    for i in eachindex(dt.defectsN)]

    locP_new    = [vortices_new[i][1:2]     for i in eachindex(vortices_new)]
    locN_new    = [antivortices_new[i][1:2] for i in eachindex(antivortices_new)]
    chargeP_new = [vortices_new[i][3]       for i in eachindex(vortices_new)]
    chargeN_new = [antivortices_new[i][3]   for i in eachindex(antivortices_new)]
    shapeP_new   = [vortices_new[i][4]       for i in eachindex(vortices_new)]
    shapeN_new   = [antivortices_new[i][4]   for i in eachindex(antivortices_new)]

    Np_new,Np_old = length(locP_new),length(locP_old)
    Nn_new,Nn_old = length(locN_new),length(locN_old)
    N_old = Np_old + Nn_old
    N_new = Np_new + Nn_new

    # Special simple cases to deal with upstream
    if N_new == N_old == 0 # do nothing, otherwise, "reducing over empty collection blablabla"

    elseif Nn_new == Nn_old == 0 && Np_new == Np_old > 0 # there are only (+) defects and no creation/annihilation
        # println("Som-hi !")
        assignment_vortices = pair_up_hungarian(dt,locP_new,locP_old,lattice,"+")
        for i in 1:Np_new update_position_and_shape!(dt.defectsP[assignment_vortices[i]],locP_new[i],shapeP_new[i]) end

    elseif Np_new == Np_old == 0 && Nn_new == Nn_old > 0 # there are only (-) defects and no creation/annihilation
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,lattice,"-")
        # println("ici")
        # println(zoom(thetas,lattice,locN_new[1]...,WINDOW)[2])
        for i in 1:Nn_new update_position_and_shape!(dt.defectsN[assignment_antivortices[i]],locN_new[i],shapeN_new[i]) end

    elseif N_new > 0 && N_old == 0
        for i in 1:Np_new add_defect!(dt,charge=chargeP_new[i],shape=shapeP_new[i],loc=locP_new[i]) end
        for i in 1:Nn_new add_defect!(dt,charge=chargeN_new[i],shape=shapeN_new[i],loc=locN_new[i]) end

    elseif N_new == 0 && N_old > 0 # (+)(-) >> plus rien
        id_just_annihilated_defectP,id_just_annihilated_defectM = ID_active_defects(dt) # seek for not yet annihilated defects

        for i in id_just_annihilated_defectP dt.defectsP[i].annihilation_time = dt.current_time end
        for i in id_just_annihilated_defectM dt.defectsN[i].annihilation_time = dt.current_time end
        # annihilation_time is already taken care of in the annihilate_defects function
        dt = annihilate_defects(dt::DefectTracker,id_just_annihilated_defectP,lattice)

    elseif Np_new > 0 && Np_old > 0 && Nn_old > 0 && Nn_new == 0  # (+)(+)(-) >> (+) par exemple
        assignment_vortices = pair_up_hungarian(dt,locP_new,locP_old,lattice,"+")
        # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_vortices) update_position_and_shape!(dt.defectsP[assignment_vortices[i]],locP_new[i],shapeP_new[i]) end
        # Identify annihilated defects
        ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
        for i in 1:number_defectsP(dt)
            if i ∉ assignment_vortices && dt.defectsP[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_vortices,i)
            end
        end
        ID_annihilated_antivortices = ID_active_defects(dt)[2] # in this special case, there is no antivortices left in the "new" timestep

        for i in ID_annihilated_vortices     dt.defectsP[i].annihilation_time = dt.current_time end
        for i in ID_annihilated_antivortices dt.defectsN[i].annihilation_time = dt.current_time end
        dt = annihilate_defects(dt,ID_annihilated_vortices,lattice)

    elseif Nn_new > 0 && Nn_old > 0 && Np_old > 0 && Np_new == 0  # (+)(-)(-) >> (-) par exemple
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,lattice,"-")
        # Update living antivortices. NB : the annihilated antivortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_antivortices) update_position_and_shape!(dt.defectsN[assignment_antivortices[i]],locN_new[i],shapeN_new[i]) end
        # Identify annihilated defects
        ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
        for i in 1:number_defectsN(dt)
            if i ∉ assignment_antivortices && dt.defectsN[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_antivortices,i)
            end
        end
        ID_annihilated_vortices = ID_active_defects(dt)[1] # in this special case, there is no vortices left in the "new" timestep
        
        for i in ID_annihilated_vortices     dt.defectsP[i].annihilation_time = dt.current_time end
        for i in ID_annihilated_antivortices dt.defectsN[i].annihilation_time = dt.current_time end
        
        dt = annihilate_defects(dt,ID_annihilated_vortices,lattice)
    elseif Nn_new > 0 && Nn_old > 0 && Np_old == 0 && Np_new > 0  # (-) >> (+)(-)(-) par exemple
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,lattice,"-")
        
        # Take care of the newly created (-) defects
        ind_created_antivortex = findall(iszero,assignment_antivortices) # newly created vortex -> the assignment vector contains a 0
        loc_created_antivortex = antivortices_new[ind_created_antivortex]
        for j in eachindex(loc_created_antivortex) add_defect!(dt,charge=chargeN_new[j],shape=shapeN_new[j],loc=loc_created_antivortex[j][1:2]) end

        # Take care of the newly created (+) defect
        loc_created_vortex = vortices_new[1]
        add_defect!(dt,charge=chargeP_new[1],shape=shapeP_new[1],loc=loc_created_vortex[1:2])

        # Update previously living antivortex
        for i in eachindex(assignment_antivortices)
            if assignment_antivortices[i] ≠ 0 # avoid newly created defects
                update_position_and_shape!(dt.defectsN[assignment_antivortices[i]],locN_new[i],shapeN_new[i])
            end
        end
    elseif Np_new > 0 && Np_old > 0 && Nn_old == 0 && Nn_new > 0  # (+) >> (+)(+)(-) par exemple
        assignment_vortices = pair_up_hungarian(dt,locP_new,locP_old,lattice,"+")
        
        # Take care of the newly created (+) defects
        ind_created_vortex = findall(iszero,assignment_vortices) # newly created vortex -> the assignment vector contains a 0
        loc_created_vortex = vortices_new[ind_created_vortex]
        for j in eachindex(loc_created_vortex) add_defect!(dt,charge=chargeP_new[j],shape=shapeP_new[j],loc=loc_created_vortex[j][1:2]) end

        # Take care of the newly created (-) defect
        loc_created_antivortex = antivortices_new[1]
        add_defect!(dt,charge=chargeN_new[1],shape=shapeN_new[1],loc=loc_created_antivortex[1:2])

        # Update previously living vortex
        for i in eachindex(assignment_vortices)
            if assignment_vortices[i] ≠ 0 # avoid newly created defects
                update_position_and_shape!(dt.defectsN[assignment_vortices[i]],locP_new[i],shapeP_new[i])
            end
        end
    else # end of special cases

    # GENERAL TREATMENT
        assignment_vortices     = pair_up_hungarian(dt,locP_new,locP_old,lattice,"+")
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,lattice,"-")

        # CASE 1 : no creation, no annihilation : simply update the data structure
        if N_new == N_old
            for i in 1:Np_new update_position_and_shape!(dt.defectsP[assignment_vortices[i]],locP_new[i],shapeP_new[i]) end
            for i in 1:Nn_new update_position_and_shape!(dt.defectsN[assignment_antivortices[i]],locN_new[i],shapeN_new[i]) end

        # CASE 2 : creation !
        elseif N_new > N_old
            # Take care of the newly created defects
            ind_created_vortex = findall(iszero,assignment_vortices) # newly created vortex -> the assignment vector contains a 0
            loc_created_vortex = vortices_new[ind_created_vortex]
            for j in eachindex(loc_created_vortex) add_defect!(dt,charge=chargeP_new[j],shape=shapeP_new[j],loc=loc_created_vortex[j][1:2]) end

            ind_created_antivortex = findall(iszero,assignment_antivortices)
            loc_created_antivortex = antivortices_new[ind_created_antivortex]
            for j in eachindex(loc_created_antivortex) add_defect!(dt,charge=chargeN_new[j],shape=shapeN_new[j],loc=loc_created_antivortex[j][1:2]) end

            # Update the ancient defects' positions
            for i in eachindex(assignment_vortices)
                if assignment_vortices[i] ≠ 0 # avoid newly created defects
                    update_position_and_shape!(dt.defectsP[assignment_vortices[i]],locP_new[i],shapeP_new[i])
                end
            end
            for i in eachindex(assignment_antivortices)
                if assignment_antivortices[i] ≠ 0 # avoid newly created defects
                    update_position_and_shape!(dt.defectsN[assignment_antivortices[i]],locN_new[i],shapeN_new[i])
                end
            end

        # CASE 3 : annihilation !
        elseif N_new < N_old
             # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
             for i in eachindex(assignment_vortices)     update_position_and_shape!(dt.defectsP[assignment_vortices[i]],locP_new[i],shapeP_new[i]) end
             for i in eachindex(assignment_antivortices) update_position_and_shape!(dt.defectsN[assignment_antivortices[i]],locN_new[i],shapeN_new[i]) end

            # Identify annihilated defects
            ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
            for i in 1:number_defectsP(dt)
                if i ∉ assignment_vortices && dt.defectsP[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                    dt.defectsP[i].annihilation_time = dt.current_time # from now on, defects that have just annihilated have annihilation_time == t
                    push!(ID_annihilated_vortices,i)
                end
            end
            for i in 1:number_defectsN(dt)
                if i ∉ assignment_antivortices && dt.defectsN[i].annihilation_time == nothing
                    dt.defectsN[i].annihilation_time = dt.current_time
                    push!(ID_annihilated_antivortices,i)
                end
            end
            dt = annihilate_defects(dt,ID_annihilated_vortices,lattice)
        end # end of general treatment
    end # end of special cases & general treatment
    return dt
end

function MSD(dfts::Union{VT1,VT2},lattice::Abstract2DLattice) where {VT1<:Vector{Union{Missing,DefectTracker}},VT2<:Vector{DefectTracker}}
    indices = [] # indices of dft defined (is simulation not finished, dfts[i] == missing)
    for i in 1:length(dfts)
        if !ismissing(dfts[i]) push!(indices,i) end
    end
    maxlength = maximum([maximum([length(d.pos) for d in vcat(dft.defectsP,dft.defectsN)]) for dft in dfts[indices]])
    MSD_P   = NaN*zeros(length(indices),maxlength)
    MSD_N   = NaN*zeros(length(indices),maxlength)
    MSD_all = NaN*zeros(length(indices),maxlength)
    for i in 1:length(indices)
        msd_all, msd_p, msd_n = MSD(dfts[indices[i]],lattice)
        MSD_P[i,1:length(msd_p)] = msd_p
        MSD_N[i,1:length(msd_n)] = msd_n
        MSD_all[i,1:length(msd_all)] = msd_all
    end

    MSD_P_avg = nanmean(MSD_P,1)[1,:]
    MSD_N_avg = nanmean(MSD_N,1)[1,:]
    MSD_all_avg = nanmean(MSD_all,1)[1,:]

    return MSD_all_avg,MSD_P_avg,MSD_N_avg
end

function MSD(dft::DefectTracker,lattice::Abstract2DLattice,maxlength=nothing)
    nP = number_defectsP(dft)
    nN = number_defectsN(dft)
    # tmin,tmax = t_bounds(dft) # (tmin,tmax) = timestamps of (first defect creation , last defect annihilation)

    # hasfield(typeof(model),:dt) ? dummy_dt = model.dt : dummy_dt = 1
    if isnothing(maxlength)
        maxlength = maximum([length(d.pos) for d in vcat(dft.defectsP,dft.defectsN)])
    end
    # Compute the SD
    SD_P = NaN*zeros(nP,maxlength)
    SD_N = NaN*zeros(nN,maxlength)
    for n in 1:nP
        defect = dft.defectsP[n]
        tmp = square_displacement(defect,lattice)
        SD_P[n,1:length(tmp)] = tmp
    end
    for n in 1:nN
        defect = dft.defectsN[n]
        tmp = square_displacement(defect,lattice)
        SD_N[n,1:length(tmp)] = tmp
    end

    # Now average to compute the MSD
    MSD_P = nanmean(SD_P,1)[1,:]
    MSD_N = nanmean(SD_N,1)[1,:]
    MSD_all = nanmean(hcat(MSD_P,MSD_N),2)[:]

    return MSD_all, MSD_P, MSD_N
end

function square_displacement(d::Defect,lattice::Abstract2DLattice)
    loc_t0 = first_loc(d)
    return [dist(lattice,loc,loc_t0) for loc in d.pos] .^ 2
end

function interdefect_distance(defect1::Defect,defect2::Defect,lattice::Union{Abstract2DLattice,TorusLattice})
    # TODO take care of case with creation and/or annihilation time different.
    # So far, this care is left to the user...
    # @assert defect1.creation_time == defect2.creation_time
    # @assert defect1.annihilation_time == defect2.annihilation_time
    tmax = min(length(defect1.pos),length(defect2.pos))
    R = [dist(lattice,defect1.pos[t],defect2.pos[t]) for t in 1:tmax]
    return R
end
interdefect_distance(dft::DefectTracker, lattice::Union{Abstract2DLattice,TorusLattice}) = interdefect_distance(dft.defectsP[1],dft.defectsN[1],lattice)

function interdefect_distance(dfts::Vector{DefectTracker},lattice::Union{Abstract2DLattice,TorusLattice})
    R = [interdefect_distance(dfts[i],lattice) for i in 1:length(dfts)]
    R_matrix = vector_of_vector2matrix(R)
    return nanmean(R_matrix,2)[:,1]
end

function mean_distance_to_annihilator(dfts::Vector{Union{Missing,DefectTracker}},lattice::Abstract2DLattice)
    indices = [] # indices of dft defined (is simulation not finished, dfts[i] == missing)
    for i in 1:length(dfts)
        if !ismissing(dfts[i]) push!(indices,i) end
    end
    Rs = [mean_distance_to_annihilator(dfts[indices[n]],lattice) for n in 1:length(indices)]
    Rs_matrix = vector_of_vector2matrix(Rs)
    return nanmean(Rs_matrix,2)[:,1]
end

function mean_distance_to_annihilator(dft::DefectTracker,lattice::Abstract2DLattice)
    nP = number_defectsP(dft)
    Rs = [distance_to_annihilator(dft,n,lattice) for n in 1:nP]
    Rs_matrix = vector_of_vector2matrix(Rs)
    return nanmean(Rs_matrix,2)[:,1]
end

function distance_to_annihilator(dft::DefectTracker,id1::Int,lattice::Abstract2DLattice;reversed=true)
    R = interdefect_distance(dft,dft.defectsP[id1],dft.defectsN[dft.defectsP[id1].id_annihilator],lattice)
    if reversed reverse!(R) end
    return R
end

