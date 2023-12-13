logspace(x1, x2, n; digits=1) = unique!(round.([10.0^y for y in range(log10(x1), log10(x2), length=n)], digits=digits))
each = eachindex # alias WARNING : "each" is at least 10x slower than "eachindex" ) 


function prinz(z)
    ss = "$(round(Int,z)) seconds"
    mm = "$(round(z/60,digits=2)) minutes"
    hh = "$(round(z/3600,digits=2)) hours"
    dd = "$(round(z/86400,digits=2)) days"
    if z < 60
        println("Runtime : $ss = $mm")
    elseif z < 3600
        println("Runtime : $mm = $hh")
    else
        println("Runtime : $hh = $dd")
    end
    return z
end


function hrun(runtimes,nbins=20)
    mean_runtime = nanmean(runtimes)
    if mean_runtime < 60
        unit = "seconds"
    elseif mean_runtime < 3600
        unit = "minutes"
        runtimes = runtimes/60
    elseif mean_runtime < 86400
        unit = "hours"
        runtimes = runtimes/3600
    else 
        unit = "days"
        runtimes = runtimes/86400
    end
    h = histogram(runtimes, nbins=nbins, title=unit)
    return h
end

nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,dims=y)
nanstd(x) = std(filter(!isnan,x))
nanstd(x,y) = mapslices(nanstd,x,dims=y)
replace_nan_with_zeros(v) = map(x -> isnan(x) ? zero(x) : x, v)
function complete_with_nan(x,n)
    @assert n ≥ length(x)
    x_extended = NaN*similar(x,n)
    x_extended[1:length(x)] = x
    return x_extended
end

circular_mean(x::AbstractVector) = mod(angle(mean(exp.(im*x))),2π)
circular_mean(x::AbstractArray, y) = mapslices(circular_mean,x,dims=y)


function vector_of_vector2matrix(x::Vector{Vector{T}}) where T<:AbstractFloat
    maxlength = maximum([length(x[i]) for i in eachindex(x)])
    matrice = NaN*zeros(T,maxlength,length(x))
    for i in eachindex(x)
        matrice[1:length(x[i]),i] = x[i]
    end
    return matrice
end

all_false(x::AbstractArray{Bool}) = iszero(sum(x))
all_true(x::AbstractArray{Bool}) = isone(prod(x))

import Base: +,-
+(x::Vector{Tuple{T,T}},y::Vector{Tuple{T,T}}) where T<:Number = [x[i] .+ y[i] for i in 1:length(x)]
-(x::Vector{Tuple{T,T}},y::Vector{Tuple{T,T}}) where T<:Number = [x[i] .- y[i] for i in 1:length(x)]
+(x::Vector{Tuple{T,T}},y::Tuple{T,T}) where T<:Number = [x[i] .+ y for i in 1:length(x)]
+(y::Tuple{T,T},x::Vector{Tuple{T,T}}) where T<:Number = [x[i] .+ y for i in 1:length(x)]
-(x::Vector{Tuple{T,T}},y::Tuple{T,T}) where T<:Number = [x[i] .- y for i in 1:length(x)]
-(y::Tuple{T,T},x::Vector{Tuple{T,T}}) where T<:Number = [x[i] .- y for i in 1:length(x)]


function index_element_repeated(x)
    indices = []
    L = length(x)
    for i in 1:L
        count = 0
        for j in 1:L
            if x[j] == x[i] count += 1 end
        end
        if count > 1 push!(indices,i) end
    end
    return indices
end

function remove_negative(input)
    array = Float32.(copy(input))
    for i in 1:length(array)
        if array[i] ≤ 0 array[i] = NaN end
    end
    return array
end
remn = remove_negative
removen = remove_negative
remne = remove_negative

function exp_moving_average(x,window)
    y = zeros(length(x))
    for t in 1:length(y)
        y[t] = sum([x[t-i]*((window-1)/window)^(i) for i in 0:t-1])/window
    end
    return y
end

function smooth(X;over=3) ## for smoother plots
    smoothed = copy(X)
    coeff = [2.0 ^(-(i-1)) for i in 1:over]
    coeff = coeff./sum(coeff)
    s = length(coeff)
    for i in 1+s:length(smoothed)
        smoothed[i] = X[i-s+1:i]'*coeff
    end
    return smoothed
end

function mu2string(x)
    if isnothing(x) || isnan(x)
        return x
    else
        return round(x, digits=2)
    end
end