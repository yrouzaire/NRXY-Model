############### XY MODEL EXAMPLE on CPU #################### 
function neighbours(thetas, i,j,L)
	return [thetas[i,mod1(j-1,L)],
			thetas[i,mod1(j+1,L)],
			thetas[mod1(i-1,L),j], 
			thetas[mod1(i+1,L),j]]
end

function update(thetas)
	L = size(thetas,1)
	thetas_old = copy(thetas)
	for i in 1:L
		for j in 1:L
			neigh = neighbours(thetas, i,j,L)
			thetas[i,j] = sum(neigh .- thetas[i,j])
		end
	end
	return thetas
end

thetas = rand(Float32, 100, 100)
update(thetas)
@btime update($thetas)

############### XY MODEL EXAMPLE on GPU #################### 
using CUDA

CUthetas = 2*pi*CUDA.rand(N, N)  # a matrix stored on the GPU filled with 1.0 (Float32)

function CUneighbours(thetas, i,j,L)
	return CUDA.@allowscalar [thetas[i,mod1(j-1,L)],
			thetas[i,mod1(j+1,L)],
			thetas[mod1(i-1,L),j], 
			thetas[mod1(i+1,L),j]]
end
CUneighbours(CUthetas, 10,10,100)


function update(thetas)
	L = size(thetas,1)
	thetas_old = copy(thetas)
	for i in 1:L
		for j in 1:L
			neigh = neighbours(thetas_old, i,j,L)
			thetas[i,j] += 0.1*sum(neigh .- thetas_old[i,j])
		end
	end
	return thetas
end


update(CUthetas)
@btime update($CUthetas)

############### BASIC EXAMPLE #################### 
# from https://cuda.juliagpu.org/stable/tutorials/introduction/#Writing-a-parallel-GPU-kernel
using CUDA
CUDA.versioninfo()
N = 2^20
x = fill(1.0f0, N)  # a vector filled with 1.0 (Float32)
y = fill(2.0f0, N)  # a vector filled with 2.0
y .+= x             # increment each element of y with the corresponding element of x

using Test
@test all(y .== 3.0f0)

## Parallelization on the CPU

function sequential_add!(y, x)
    for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end
    return nothing
end

fill!(y, 2)
sequential_add!(y, x)
@test all(y .== 3.0f0)

function parallel_add!(y, x)
    Threads.@threads for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end
    return nothing
end

fill!(y, 2)
parallel_add!(y, x)
@test all(y .== 3.0f0)

@btime sequential_add!($y, $x)
@btime parallel_add!($y, $x)
	
## Your first GPU computation
using CUDA
x_d = CUDA.fill(1.0f0, N)  # a vector stored on the GPU filled with 1.0 (Float32)
y_d = CUDA.fill(2.0f0, N)  # a vector stored on the GPU filled with 2.0

y_d .+= x_d
@test all(Array(y_d) .== 3.0f0)

function add_broadcast!(y, x)
    CUDA.@sync y .+= x
    return
end
@btime add_broadcast!($y_d, $x_d)

## Writing your first GPU kernel
function gpu_add1!(y, x)
    for i = 1:length(y)
        @inbounds y[i] += x[i]
    end
    return nothing
end

fill!(y_d, 2)
@cuda gpu_add1!(y_d, x_d)
@test all(Array(y_d) .== 3.0f0)

function bench_gpu1!(y, x)
    CUDA.@sync begin
        @cuda gpu_add1!(y, x)
    end
end
@btime bench_gpu1!($y_d, $x_d)

## Writing a parallel GPU kernel
function gpu_add2!(y, x)
    index = threadIdx().x    # this example only requires linear indexing, so just use `x`
    stride = blockDim().x
    for i = index:stride:length(y)
        @inbounds y[i] += x[i]
    end
    return nothing
end

fill!(y_d, 2)
@cuda threads=256 gpu_add2!(y_d, x_d)
@test all(Array(y_d) .== 3.0f0)

function bench_gpu2!(y, x)
    CUDA.@sync begin
        @cuda threads=256 gpu_add2!(y, x)
    end
end
@btime bench_gpu2!($y_d, $x_d)

function gpu_add3!(y, x)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for i = index:stride:length(y)
        @inbounds y[i] += x[i]
    end
    return
end

numblocks = ceil(Int, N/256)

fill!(y_d, 2)
@cuda threads=256 blocks=numblocks gpu_add3!(y_d, x_d)
@test all(Array(y_d) .== 3.0f0)

function bench_gpu3!(y, x)
    numblocks = ceil(Int, length(y)/256)
    CUDA.@sync begin
        @cuda threads=256 blocks=numblocks gpu_add3!(y, x)
    end
end
@btime bench_gpu3!($y_d, $x_d)

############## Now let's do the same but for a 2D array
N = 2^10
x = fill(1.0f0, N, N)  # a matrix filled with 1.0 (Float32)
y = fill(2.0f0, N, N)  # a matrix filled with 2.0
y .+= x             # increment each element of y with the corresponding element of x

using Test
@test all(y .== 3.0f0)

## Parallelization on the CPU
function sequential_add!(y, x)
    for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end
    return nothing
end

y = fill(2.0f0, N, N)  # a matrix filled with 2.0
sequential_add!(y, x)
@test all(y .== 3.0f0)

function parallel_add!(y, x)
    Threads.@threads for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end
    return nothing
end

y = fill(2.0f0, N, N)  # a matrix filled with 2.0
parallel_add!(y, x)
@test all(y .== 3.0f0)

@btime sequential_add!($y, $x)
@btime parallel_add!($y, $x)

## Parallelization on the GPU
using CUDA
x_d = CUDA.fill(1.0f0, N, N)  # a matrix stored on the GPU filled with 1.0 (Float32)
y_d = CUDA.fill(2.0f0, N, N)  # a matrix stored on the GPU filled with 2.0
y_d .+= x_d
@test all(Array(y_d) .== 3.0f0)

function add_broadcast!(y, x)
    CUDA.@sync y .+= x
    return
end
@btime add_broadcast!($y_d, $x_d)


function gpu_add1!(y, x)
    for i = 1:size(y,1)
		for j = 1:size(y,1)
			@inbounds y[i,j] += x[i,j]
		end
    end
    return nothing
end

y_d = CUDA.fill(2.0f0, N, N)  # a matrix stored on the GPU filled with 2.0
@cuda gpu_add1!(y_d, x_d)
@test all(Array(y_d) .== 3.0f0)
@btime bench_gpu1!($y_d, $x_d)


function gpu_add2!(y, x)
    index_i = threadIdx().x    # this example only requires linear indexing, so just use `x`
    stride_i = blockDim().x
	index_j = threadIdx().y    # this example only requires linear indexing, so just use `x`
	stride_j = blockDim().y
    for i = index_i:stride_i:size(y,1)
		for j = index_j:stride_j:size(y,1)
			@inbounds y[i,j] += x[i,j]
		end
	end
    return nothing
end

y_d = CUDA.fill(2.0f0, N, N)  # a matrix stored on the GPU filled with 2.0
@cuda threads=(16,16) gpu_add2!(y_d, x_d)
@test all(Array(y_d) .== 3.0f0)
function bench_gpu2!(y, x)
    CUDA.@sync begin
        @cuda threads=256 gpu_add2!(y, x)
    end
end
@btime bench_gpu2!($y_d, $x_d)


function gpu_add3!(y, x)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for i = index:stride:length(y)
        @inbounds y[i] += x[i]
    end
    return
end

numblocks = ceil(Int, 1/256)

fill!(y_d, 2)
@cuda threads=256 blocks=numblocks gpu_add3!(y_d, x_d)
@test all(Array(y_d) .== 3.0f0)