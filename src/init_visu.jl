include("lattices.jl");
include("models.jl");
include("core_methods.jl");
include("defects_methods.jl");
import StatsBase.sample
import Plots.@animate

## ------------------------ Initializations  ------------------------
function init_thetas(model::AbstractPropagationModel{T}, lattice::Abstract1DLattice, params)::Vector{T} where {T<:AbstractFloat}
    L = lattice.L
    @unpack init = params
    if init == "lowtemp"
        thetas = zeros(L)
    elseif init == "spinwave"
        model.symmetry == "polar" ? symm = 2 : symm = 1
        thetas = [symm * pi / L * i for i in 1:L]
    end
    model.omegas = sqrt(model.Var) * randn(L)
    return thetas
end

function init_thetas(model::AbstractPropagationModel{T}, lattice::Abstract2DLattice,params=nothing)::Matrix{T} where {T<:AbstractFloat}
    L = lattice.L
    # for now, only lowtemp initialisation is supported, so no need to provide params_init
    model.omegas = sqrt(model.Var) * randn(L, L)
    thetas = zeros(L, L)
    return thetas
end

function init_thetas(model::AbstractModel{float_type}, lattice::Abstract2DLattice, params) where {float_type<:AbstractFloat}
    L = lattice.L
    @unpack init, q, r0, mu0, mu1, decay_function, mu_plus, mu_minus, phi = params
    if init in ["hightemp", "disorder"]
        thetas = 2π * rand(L, L)
    elseif init in ["lowtemp", "polar_order"]
        thetas = zeros(L, L)
    elseif init in ["lowtemp_nematic", "nematic_order"]
        thetas = rand(Bool, L, L) * π
    elseif init in ["spinwave","spinwaves"]
        thetas = create_spinwave(L,nb_horizontal_waves,nb_vertical_waves)
    elseif init in ["isolated", "single"]
        # For a single untwisted defect, theta(x,y)  = q atan(y/x) + μ0 
        x0 , y0 = location_core_defect(L, lattice)
        thetas = create_single_defect(lattice, x0=x0,y0=y0, q=q, mu=mu0) # in case of something more exotic, recall that the use is : create_single_defect(q,type,L,y0,x0) (x and y swapped)
        lattice.periodic = false
    elseif init in ["single_twist","single_twisted", "twisted", "twist"]
        # For a single twisted defect, theta(x,y)  = q atan(y/x) + μ0 + μ1 * decay_function(r)
        x0 , y0 = location_core_defect(L, lattice)
        thetas = create_single_defect(lattice, x0=x0,y0=y0, q=q, mu=mu0) # in case of something more exotic, recall that the use is : create_single_defect(q,type,L,y0,x0) (x and y swapped)
        distance_to_core = [sqrt((i - x0)^2 + (j - y0)^2) for i in 1:L, j in 1:L]
        thetas += mu1 * decay_function.(distance_to_core)
        lattice.periodic = false
    elseif init == "pair"
        thetas = create_pair_vortices(lattice, r0=r0, q=abs(q), types=[mu_plus, mu_minus, phi])
        thetas = relax(thetas, model, lattice, trelax = 2,T=0)
    elseif init in ["flag","step","band"]
        thetas = 0*ones(L, L) # background color 
        width = L/5
        thetas[:, round(Int, L/2-width/2):round(Int, L/2+width/2)] .= -π/2 # stripe color
    elseif init in ["smooth_flag","smooth_step","smooth_band"]
        thetas = zeros(L, L) # background color 
        width = 30
        for y in 1:L
            thetas[:,y] .= -π/4 * exp(-(y - round(Int, L / 2))^2 / (2(width/2)^2))
        end
    elseif init in ["pair_torsed", "pair_torsion"]
        torsion_time = 200
        torsion_sigma = 0.3
        params_torsion = Dict("init" => "pair", "q" => q, "r0" => r0,
            "phi" => phi, "mu0" => mu0, "mu_minus" => mu_minus,
            "mu_plus" => mu_plus, "L" => L, "T" => 0,
            "sigma" => torsion_sigma, "symmetry" => model.symmetry,
            "dt" => model.dt, "float_type" => float_type, "rho" => model.rho)
        model_torsion = NRXY(params_torsion)
        thetas = init_thetas(model_torsion, lattice,params_torsion)
        update!(thetas, model_torsion, lattice, torsion_time)
    else
        error("ERROR : Type of initialisation unknown. Choose among \"hightemp/order\",\"lowtemp/polar_order\",\"isolated\" , \"pair\" or \"lowtemp_nematic/nematic_order\" .")
    end
    if model.rho < 1
        make_holes!(thetas, model.rho)
    end

    thetas = float_type.(thetas)

    return thetas
end

function location_core_defect(L, lattice::Abstract2DLattice)
    if lattice isa SquareLattice
        x0 = (L + 1*iseven(L))/ 2
        y0 = (L + 1*iseven(L))/ 2
    elseif lattice isa TriangularLattice
        if iseven(L)
            if iseven(L/2)
                x0 = (L+1)/ 2 
                y0 = (L-1/sqrt(3))/ 2
            else
                x0 = (L+1)/ 2 
                y0 = (L+1/sqrt(3))/ 2
            end
        else
            if iseven((L-1)/2)
                x0 = (L+1)/ 2 
                y0 = (L+1/sqrt(3))/ 2
            else
                x0 = (L)/ 2 
                y0 = (L+1+1/sqrt(3))/ 2
            end
        end
    end
    return x0, y0
end

function create_spinwave(L,nb_horizontal_waves, nb_vertical_waves)
    # 1D spinwave
    thetas = zeros(L,L)
    horizontal_wave = nb_horizontal_waves ≠ 0 && nb_vertical_waves == 0
    vertical_wave = nb_vertical_waves ≠ 0 && nb_horizontal_waves == 0
    if horizontal_wave
        println("horizontal wave")
        x_wave = [mod(nb_horizontal_waves*j*2π/L,2*π) for j in 1:L]
        for i in 1:L
            thetas[i,:] = x_wave
        end
    elseif vertical_wave
        println("vertical wave")
        y_wave = [mod(nb_vertical_waves*i*2π/L,2*π) for i in 1:L]
        for j in 1:L
            thetas[:,j] = y_wave
        end
    else
        println("Diagonal wave")
        x_wave = [mod(nb_horizontal_waves*j*2π/L,2π) for j in 1:L]
        y_wave = [mod(nb_vertical_waves*i*2π/L,2π) for i in 1:L]
        for i in 1:L, j in 1:L
            thetas[i,j] = mod(x_wave[j] + y_wave[i], 2π)
        end
    end
    #= 2D spinwave, like in Ballistic Kuramoto project.
    Useless because even for XY dynamics, it goes back 
        to ordered state.
        Code kept for reference.
    for i in 1:L, j in 1:L
        tmp = π + π/2 * (sin(nb_vertical_waves * i* 2π / L) + sin(nb_horizontal_waves * j * 2π / L))
        thetas[i,j] = mod(tmp, 2π)
    end
    =#
    
    return thetas    
end
create_spinwaves = create_spinwave # alias for convenience 

function create_single_defect(lattice;x0, y0, q=1, mu=2π*rand())
    L = lattice.L
    thetas = zeros(L, L)
    for i in 1:L, j in 1:L
        x,y = ij2xy(i, j, L, lattice)
        dy = y - y0
        dx = x - x0
        #= Because the two folloying lines are not understood 
        by the julia version I have on the cluster, 
        I recode them.
        dy = argmin(abs,[dy, L - dy, -L + dy])
        dx = argmin(abs,[dx, L - dx, -L + dx]) =#
        arr = [dy, L - dy, -L + dy]
        min_abs_value = abs(dy)
        for i in 2:length(arr)
            if abs(arr[i]) < min_abs_value
                dy = arr[i]
                min_abs_value = abs(arr[i])
            end
        end
        arr = [dx, L - dx, -L + dx]
        min_abs_value = abs(dx)
        for i in 2:length(arr)
            if abs(arr[i]) < min_abs_value
                dx = arr[i]
                min_abs_value = abs(arr[i])
            end
        end

        thetas[i, j] = q * atan(dy,dx)
    end
    return thetas .+ mu
    #= Important Notes : 
    A. Here I use the atan2(y,x) function instead of atan(y/x) 
        because atan(y/x) gives the wrong field for x<0 
    B. Recall that the mapping from the position (i,j) in the matrix 
        and the position (x,y) when visualized is given by 
        i = L - y + 1
        j = x
        equivalent to 
        x = j
        y = L - i + 1
        This is encoded in the function ij2xy(i,j,L) =#
end

function create_pair_vortices(lattice::Abstract2DLattice; r0=Int(L / 2), q, types)
    L = lattice.L
    #= Check for meaningfulness of the defaults separation,
    otherwise the defaults will annihilate with relaxation =#
    @assert 2 ≤ r0 ≤ 0.5L "Error : r0 should be bigger 2 ≤ r0 ≤ L/2. "
    if isodd(r0)
        r0 -= 1
        # println("Warning : r0 has to be even, corrected : r0 -= 1 ")
    end

    mu_plus, mu_minus, phi = types
    @assert sum(isnothing.(types)) == 1

    #= IMPORTANT ! 
    
    The following lines arise from the following continuity equation : mu_minus - mu_plus = π + 2φ, 
    to be understood as modulo 2π. φ is defined as the angle between the horizontal axis and the line 
    joining the + defect to the - defect.
   
    Because θ = ± atan(y/x) + µ = ± φ + µ, one has
        µ+ = µ1 + µ2 - φ - π
        µ- = µ1 + µ2 + φ
    Since µ1 and µ2 only enter the equation as their sum, 
    one can set µ1 = 0 and only work with µ2, obtaining:
        µ+ = µ2 - φ - π
        µ- = µ2 + φ

    or reversed, depending on which information is available:
        µ2 = µ+ + φ + π
        µ2 = µ- - φ
    
    Leading to the (equivalent) equations : 
        µ+ - µ- = - 2φ - π
        µ- - µ+ = 2φ - π
    =#
    mu1 = 0
    if isnothing(phi) # mu_plus and mu_minus are given/imposed
        phi = (mu_minus - mu_plus - π) / 2
        mu2 = mu_minus - phi # could also be mu_plus + phi + π 
        
        #= Two side comments : 
        About the 1st line: 
        Keep in mind that phi =  (mu_minus - mu_plus - pi) / 2 + π 
        could also be a solution. It changes the field and inverts the 
        position of the defects but their rescpectives µ are unchanged. 
        For experiment reproducibility, I drop the additional constant "+ rand([0, π])". 
        Nevertheless, it is possible to obtain the configuration one would have obtained with +π. 
        One just has to change µ- → µ- + π  
        (for the negative defect only, since I use the equation from µ- right below.) 
       
        About the 2nd line:  
        It could very well have been mu2 = mu_plus + phi + π. In this case, 
        one should change µ+ → µ+ + π if one wants to recover the effect of φ → φ + π.
        =#
    
    elseif isnothing(mu_plus) # phi and mu_minus are given/imposed
        mu2 = mu_minus - phi
    elseif isnothing(mu_minus) # phi and mu_plus are given/imposed
        mu2 = mu_plus + phi + π
    end

    # Location of the defects
    xp = round(Int, L / 2 + r0 / 2 * cos(phi - π))
    yp = round(Int, L / 2 + r0 / 2 * sin(phi - π))
    xm = round(Int, L / 2 + r0 / 2 * cos(phi))
    ym = round(Int, L / 2 + r0 / 2 * sin(phi))

    thetas = create_single_defect(lattice, x0=xp, y0=yp, q=+q, mu=mu1) +
             create_single_defect(lattice, x0=xm, y0=ym, q=-q, mu=mu2)

    # return smooth_border!(thetas)
    return thetas
end

function smooth_border!(thetas)
    cst = round(Int, 0.02 * size(thetas, 1)) # has to be even

    for i in 1:L
        sample1 = thetas[i, cst]
        sample2 = thetas[i, L-cst]
        arcl = arclength(sample1, sample2, 2pi) # Oddly, pi or 2pi is not important even when symmetry is polar or nematic
        smoothed_values = [sample1 + arcl * n / (2cst + 1) for n in 1:(2cst+1)]
        for j in 1:2cst+1
            thetas[i, mod1(Int(cst - j + 1), L)] = smoothed_values[j]
        end
    end
    return thetas # in fact, this procedure is fairly insensible to relax() since already smooth
end

function create_2pairs_vortices(L; r0, q, types)
    mu_plus, mu_minus, phi = types
    @assert isinteger(r0 / 4) "R0 has to be a multiple of 4 !"
    return create_pair_vortices(L, r0=r0, q=q, phi=phi, type=type) + create_pair_vortices(L, r0=r0 / 2, q=q, phi=phi, type=type)'
end

function make_holes!(thetas, rho)
    N = length(thetas)
    holes = sample(1:N, round(Int, (1 - rho) * N), replace=false)
    thetas[holes] .= NaN
    return thetas
end

## ------------------------ Initialisation on the Torus  ------------------------ ##
## ------------------------ Initialisation on the Torus  ------------------------ ##
## ------------------------ Initialisation on the Torus  ------------------------ ##

function init_thetas(model::AbstractModel{float_type}, lattice::TorusLattice, params) where {float_type<:AbstractFloat}
    L = lattice.L
    # b = lattice.b
    @unpack init = params
    if init in ["hightemp", "disorder"]
        psis = 2π * rand(L, L)
    elseif init in ["lowtemp", "polar_order"]
        psis = zeros(L, L)
    elseif init in ["pair"]
        psis = create_pair_vortices(lattice, params)
    else 
        error("ERROR : Type of initialisation unknown. Choose among \"hightemp/order\",\"lowtemp/polar_order\" , \"single\" , \"pair\" or \"spinwave_horizontal\"  \"spinwave_vertical\" .")
    end
    psis = float_type.(psis)

    return psis
end

function init_thetas(model::AbstractModel{float_type}, lattice::TorusHydroLattice, params) where {float_type<:AbstractFloat}
    L = lattice.L
    # b = lattice.b
    @unpack init, xpu = params
    if init in ["hightemp", "disorder"]
        psis = 2π * rand(L, L)
    elseif init in ["lowtemp", "polar_order"]
        psis = zeros(L, L)
    elseif init in ["pair"]
        psis = create_pair_vortices(lattice, params)
    elseif init in ["horizontal", "spinwave_horizontal"]
        psis = zeros(L, L)
        for i in 1:L
            psis[i,:] .= 2π * (i-1) / L
        end
    elseif init in ["spinwave_vertical", "vertical"]
        psis = zeros(L, L)
        for j in 1:L
            psis[:,j] .= 2π * (j-1) / L
        end
    else 
        error("ERROR : Type of initialisation unknown. Choose among \"hightemp/order\",\"lowtemp/polar_order\" , \"single\" , \"pair\" or \"spinwave_horizontal\"  \"spinwave_vertical\" .")
    end
    psis = float_type.(psis)

    S = float_type.(1) # arbitrary initial length
    thetas = S * cos.(psis) # the normalisation √gθθ = 1 so it doesnt matter.
    phis   = S * sin.(psis) # Here the normalisation √gφφ = ρ = η + cos(θ) so we have to rescale phis
    for i in 1:L
        phis[i,:] /= lattice.rho[i]
    end
    if xpu =="cpu" return thetas, phis
    elseif xpu == "gpu" return CuArray(thetas), CuArray(phis)
    end
end

function create_pair_vortices(lattice::AbstractTorusLattice, params)
    L = lattice.L
    @unpack r0, q, mu_plus, theta0 = params # phi = ±π/2 by default, so µ_minus will be computed from µ_plus
    phi = π/2
    mu2 = mu_plus + phi + π
    
    #= Check for meaningfulness of the defaults separation,
    otherwise the defaults will annihilate with relaxation =#
    @assert 2 ≤ r0 ≤ 0.5L "Error : r0 should be bigger 2 ≤ r0 ≤ L/2. "
    if isodd(r0)
        r0 -= 1
        # println("Warning : r0 has to be even, corrected : r0 -= 1 ")
    end


    # Location of the defects. Here x <-> φ and y <-> θ. 0 < theta0 < 2π 
    xp = xm = round(Int, L / 2) 
    yp = mod1(round(Int, theta0/2π*L),L)
    ym = mod1(round(Int, theta0/2π*L + r0),L)
    if ym < yp 
        mu2 -= π # because now the negative defect is below the positive one
    end
    thetas = create_single_defect(lattice, x0=xp, y0=yp, q=+q, mu=mu1) +
             create_single_defect(lattice, x0=xm, y0=ym, q=-q, mu=mu2)

    # return smooth_border!(thetas)
    return thetas
end

## ------------------------ Visualization  ------------------------
function plot_thetas(thetas::Vector{T}, model::AbstractModel, lattice::Abstract1DLattice; cols=cgrad([:black, :blue, :green, :orange, :red, :black]),size=(400, 100)) where {T<:AbstractFloat}
    if model.symmetry == "nematic"
        modd = pi
    elseif model.symmetry == "polar"
        modd = 2pi
    end
    # thetas_fattened = hcat(thetas,thetas,thetas)
    p = heatmap(mod.(thetas', modd), c=cols, clims=(0, modd), size=size, yaxis=nothing)
    return p
end

function plot_thetas(thetas::Matrix{<:AbstractFloat}, model::AbstractModel, lattice::Union{Abstract2DLattice,TorusLattice}; defects=false, quiver=false, axis=true,  force_quiver=false, title="", colorbar=true, cols=cgrad([:black, :blue, :green, :orange, :red, :black]), size=(400 + colorbar * 85, 400))
    if model.symmetry == "nematic"
        modd = pi
    elseif model.symmetry == "polar"
        modd = 2pi
    end

    p = plot_thetas_auxx(thetas, model, lattice, cols=cols, 
        size=size, colorbar=colorbar, colorbartitle="θ", axis=axis)
    title!(title)
    
    if defects
        defects_p, defects_m = spot_defects(thetas, model, lattice)
        locP = [defects_p[i][1:2] for i in each(defects_p)]
        locN = [defects_m[i][1:2] for i in each(defects_m)]
        highlight_defects!(p, lattice.L, locP, locN)
    end

    garde_fou_quiver = 25
    L = round(Int, (sqrt(length(thetas)) - 1) / 2, RoundDown)
    if quiver
        if (L < garde_fou_quiver) || force_quiver
            display_quiver!(p, thetas, lattice, window=L)
        elseif (L > garde_fou_quiver) || force_quiver == false
            println("Quiver procedure stopped because the system size is large.\n
                Set the kwarg force_quiver to `true` to force the quiver procedure.\n
                Note that it might take a while to display the figure as
                PyPlot \"quiver\" function is awfully slow.\n")
        elseif (L > garde_fou_quiver) || force_quiver == true
            println("You wanted it ! It might take a long while to display the figure.")
        end
    end
    # xlims!((0,lattice.L))
    # ylims!((0,lattice.L))
    return p
end

function plot_thetas(thetas, phis, model, lattice::TorusHydroLattice;size=(512,512), vertical = false) 
    psis, S2 = get_psis_S2(thetas, phis, lattice)
    dummy_lattice = TorusLattice(lattice.L)
    p1 = plot_thetas(psis, model, dummy_lattice)   
    p2 = heatmap(reverse(S2,dims=1),size=size, aspect_ratio = 1, clims =(0,1))
    ss = size[1]
    if vertical 
        p = plot(p1,p2,layout=(2,1), size=(ss,2ss)) # vertical
    else
        p = plot(p1,p2,layout=(1,2), size=(2ss,ss)) # horizontal
    end
    return p
end

ij2xy(i, j, L, lattice::SquareLattice) = (j, L - i + 1)
ij2xy(i, j, L, lattice::AbstractTorusLattice) = (j, L - i + 1)
ij2xy(i, j, L, lattice::TriangularLattice) = (j + 0.5 * iseven(i), L - i + 1)
function plot_thetas_auxx(thetas::Matrix{<:AbstractFloat}, model::AbstractModel, lattice::Union{Abstract2DLattice,TorusLattice}; colorbartitle="θ", colorbar=true, cols=cgrad([:black, :blue, :green, :orange, :red, :black]), axis=true,  size=(400 + colorbar * 85, 400))
    if model.symmetry == "nematic"
        modd = π
    elseif model.symmetry == "polar"
        modd = 2π
    end

    L = Int(sqrt(length(thetas))) # cannot use size(thetas,1) because "size" overridden by kwarg 
    if isa(lattice, SquareLattice) || isa(lattice, TorusLattice)
        markerr = :square
        asp = 1 # because square
        mss = 115 / L * size[2]^1.01 / 460 # benchmarked for gr()
        # mss = 165 / L * size[2]^1.01 / 460 # benchmarked for pyplot()
    elseif isa(lattice, TriangularLattice)
        # myhex = [(1, 0), (0.5, 0.866), (-0.5, 0.866), (-1, 0), (-0.5, -0.866), (0.5, -0.866), (1, 0)]
        myhex2 = [(cos(π/6),sin(π/6)),(0,1),(-cos(π/6),sin(π/6)),(-cos(π/6),-sin(π/6)),(0,-1),(cos(π/6),-sin(π/6)),(cos(π/6),sin(π/6))]
        markerr = Shape(myhex2)
        asp = 0.866 # because of the shape of an hexagon is not isotropic
        mss = 405 / L * size[2]^1.02 / 600 # benchmarked for gr()
        # mss = 157 / L * size[2]^1.15 / 460  # benchmarked for pyplot()
    end
    
    p = plot(xlims=(1, L + 0.5), ylims=(1, L + 0.5),aspect_ratio=asp, size=size, axis=axis)
    locs = vec([(ij2xy(i, j, L, lattice)) for i in 1:L, j in 1:L])
    couleurs = vec(mod.(thetas, modd))
    scatter!(locs, marker_z=couleurs, c=cols, clims=(0,modd), marker=markerr, markersize=mss, colorbartitle=colorbartitle, colorbar=colorbar)

    return p
end

function plot_mus(dft::DefectTracker, L::Int; colorbar=true, cols=cgrad([:black, :blue, :green, :orange, :red, :black]), size=(400 + colorbar * 85, 400), vel=false, acc=false, over=1)
    factor_vel = 3
    factor_acc = 1
    p = plot(colorbar=true, clims=(0, 2pi), colorbartitle="µ", aspect_ratio=1)
    for defect in dft.defectsP
        pos = last_loc(defect)
        plot!(pos, markershape=:circle, markersize=7, marker_z=defect.type / 2π, c=cols)
        if vel
            speed = velocity(defect, lattice, over=over)
            quiver!([pos[1]], [pos[2]], quiver=([pos[1] + factor_vel * speed[1]], [pos[2] + factor_vel * speed[2]]), c=:black)

        end
        if acc
            accel = acceleration(defect, lattice, over=over)
            quiver!([pos[1]], [pos[2]], quiver=([pos[1] + factor_acc * accel[1]], [pos[2] + factor_vel * accel[2]]), c=:red)
        end
    end
    for defect in dft.defectsN
        pos = last_loc(defect)
        plot!(pos, markershape=:utriangle, markersize=7, marker_z=defect.type / 2π, c=cols)
        if vel
            speed = velocity(defect, lattice, over=over)
            quiver!([pos[1]], [pos[2]], quiver=([pos[1] + factor_vel * speed[1]], [pos[2] + factor_vel * speed[2]]), c=:black)

        end
        if acc
            accel = acceleration(defect, lattice, over=over)
            quiver!([pos[1]], [pos[2]], quiver=([pos[1] + factor_acc * accel[1]], [pos[2] + factor_vel * accel[2]]), c=:red)
        end
    end

    xlims!(0, L)
    ylims!(0, L)
    return p
end

function plot_mus(thetas::Matrix{<:AbstractFloat}, model::AbstractModel, lattice::Abstract2DLattice)
    dft = DefectTracker(thetas, model, lattice, find_type=true)
    return plot_mus(dft, lattice.L)
end

function highlight_defects!(p, L, defects_p, defects_m, symbP=:circle, symbM=:utriangle)
    dummy_lattice = SquareLattice(L)
    # dummy_lattice = TriangularLattice(L)
    for defect in defects_p
        scatter!(ij2xy(defect[1:2]..., L, dummy_lattice), m=(8, 1.0, symbP, :transparent, stroke(1.2, :grey85)))
    end
    for defect in defects_m
        scatter!(ij2xy(defect[1:2]..., L, dummy_lattice), m=(8, 1.0, symbM, :transparent, stroke(1.2, :grey85)))
    end
    # xlims!((1,L))
    # ylims!((1,L))
    return p
end

function display_quiver!(p, thetas_zoom, lattice;window)
    p
    # for j in 1:2window+1
    #     quiver!(j * ones(2window + 1), collect(1:2window+1),
    #         quiver=(-sin.(thetas_zoom[j, :]), cos.(thetas_zoom[j, :])),
    #         # quiver=(cos.(thetas_zoom[j, :]), sin.(thetas_zoom[j, :])),
    #         c=:white, lw=0.8)
    # end

    LL = size(thetas_zoom)[1]
    if isa(lattice, SquareLattice)
        for y in 1:LL
            quiver!(collect(1:2window+1), y * ones(2window + 1),
                quiver=(cos.(thetas_zoom[LL-y+1,:]), sin.(thetas_zoom[LL-y+1,:])),
                c=:white, lw=0.8)
        end
    elseif isa(lattice, TriangularLattice)
        if iseven(lattice.L)
            for y in 1:LL
                quiver!(collect(1:2window+1).+0.5*isodd(y), y * ones(2window + 1),
                    quiver=(cos.(thetas_zoom[LL-y+1,:]), sin.(thetas_zoom[LL-y+1,:])),
                    c=:white, lw=0.8)
                end
        else
            for y in 1:LL
                quiver!(collect(1:2window+1).+0.5*iseven(y), y * ones(2window + 1),
                    quiver=(cos.(thetas_zoom[LL-y+1,:]), sin.(thetas_zoom[LL-y+1,:])),
                    c=:white, lw=0.8)
                end
        end
    end
    return p
end
# Note : after returning the plots with quiver, one has to add
# xlims!(1,2window+1) ; ylims!(1,2window+1)


function zoom_quiver(thetas, model, lattice::Abstract2DLattice, i, j, window=WINDOW; defects=false, size=(400, 400))
    L = lattice.L
    no_problem_go_ahead, thetas_zoom = zoom(thetas, lattice, j,i, window) # j<->i on purpose
    if no_problem_go_ahead
        p = plot_thetas(thetas_zoom, model, lattice, defects=defects, size=size)
        display_quiver!(p, thetas_zoom, lattice, window=window)
        xlims!(1, 2window + 1)
        ylims!(1, 2window + 1)
    else
        error("Zoom not possible")
    end
    return p
end


## ------------------------ Movies  ------------------------
function movies(thetas, model, lattice; defects=false, saving_times, transients=Inf, title=titre)
    anim = @animate for time in saving_times
        println("$(round(time/saving_times[end]*100,digits=2)) %")
        evolve!(thetas, model, lattice, time)  # updates until time = t
        if time < transients
            p = plot_thetas(thetas, model, lattice, defects=false, size=(512, 512))
        else
            p = plot_thetas(thetas, model, lattice, defects=defects, size=(512, 512))
        end
        title!(titre * "\n t = $(round(time,digits=1))")
    end
    return anim
end
