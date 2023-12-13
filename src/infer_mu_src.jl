## ------------- New Definition of µ ------------- ##
# It is now local and not based on the field around an ideal and isotropic defect)
# function infer_mu(thetas_neighbours,lattice)
#     if isa(lattice,TriangularLattice)
#         constant = π/3
#     elseif isa(lattice,SquareLattice)
#         constant = π/2
#     end
#     mu_hats = [thetas_neighbours[i] - i*constant for i in eachindex(thetas_neighbours)]
#     mu_hat = angle(sum([exp(im*mu) for mu in mu_hats]))
#     mu_hat = mod(mu_hat,2π)
#     return mu_hat
# end


# ## Old Functions
# function infer_mu(thetas;q,window=WINDOW,decay=true)
#     if decay && q > 0 return infer_mu_decay(thetas,q=q,window=window)
#     else              return infer_mu_no_decay(thetas,q=q,window=window)
#     end
# end
# function infer_mu_no_decay(thetas::Matrix{T};q,window=WINDOW) where T<:AbstractFloat
#     L = size(thetas,1)
#     @assert L == 2window+1
#     muss = zeros(size(thetas))
#     # tmp = Complex(0)
#     range = 2:L-1
#     for j in range, i in range
#         muss[i,j] = thetas[i,j] - abs(q)*atan( (i-window) ,(j-window))
#         # i<->j irrelevant because i,j and j,i have the same weight for "mean" operation
#     end
#     moyenne = angle(mean(exp.(im*muss[range,range])))
#     if     q == +1   correction = -pi/2
#     elseif q == -1   correction = 0.01
#     elseif q == +1/2 correction = 0.1
#     elseif q == -1/2 correction = pi/4
#     end
#     #= To be honest, I don't know where the shifts might come from,
#     In the beggining, I thought maybe from the spin at the center of the defect,
#     where theta should not be defined. But if one changes mean->sum and adds the condition
#     "if (i == j == window)", the result only becomes weirder... =#
#     # return mod.(muss,2pi)
#     return mod(moyenne .+ correction,2π)
# end
# function infer_mu_decay(thetas::Matrix{T};q,window=WINDOW) where T<:AbstractFloat
#     L = size(thetas,1)
#     @assert L == 2window+1
#     muss = zeros(size(thetas))
#     tmp = Complex(0)
#     range = 2:L-1
#     for j in range, i in range
#         muss[i,j] = thetas[i,j] - abs(q)*atan( (i-window) ,(j-window))
#         # i<->j irrelevant because i,j and j,i have the same weight for "mean" or "sum" operation
#         distance = sqrt((i-window)^2 + (j-window)^2)
#         # if distance == 0 distance = Inf end # does horribly wrong
#         tmp += exp(im*muss[i,j] - 1*distance) # 1*distance seems to be the best
#     end
#     moyenne = angle(tmp)
#     if     q == 1   correction = -pi/2 - 0.4
#     elseif q == -1  correction = 0.25
#     elseif q == 1/2 correction = 0.1
#     elseif q == -1/2 correction = pi/4
#     end
#     #= To be honest, I don't know where the shifts might come from,
#     In the beggining, I thought maybe from the spin at the center of the defect,
#     where theta should not be defined. But if one changes mean->sum and adds the condition
#     "if (i == j == window)", the result only becomes weirder... =#
#     return mod(moyenne .+ correction,2π)
# end
