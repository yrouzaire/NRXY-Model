## Provide Divergence and Rotationnal
function provide_div_rot(X::Array{T,2}) where T<:AbstractFloat
    LL = size(X,1)
    latt = TriangularLattice(LL,periodic=false)
    X_divrot = zeros(T,LL,LL,3)
    X_divrot[:,:,1] = X
    divmat, rotmat = get_div_rot(X,latt)
    X_divrot[:,:,2] = divmat
    X_divrot[:,:,3] = rotmat
    return X_divrot
end

function provide_div_rot(X::Array{T,3}) where T<:AbstractFloat
    LL = size(X,1)
    M  = size(X,3)
    latt = TriangularLattice(LL,periodic=false)
    X_divrot = zeros(T,LL,LL,3,M)
    for m in 1:M
        X_divrot[:,:,1,m] = X[:,:,m]
        divmat, rotmat = get_div_rot(X[:,:,m],latt)
        X_divrot[:,:,2,m] = divmat
        X_divrot[:,:,3,m] = rotmat
    end
    return X_divrot
end
function provide_div_rot(X::Array{T,4}) where T<:AbstractFloat
    LL = size(X,1)
    M  = size(X,4)
    latt = TriangularLattice(LL,periodic=false)
    X_divrot = zeros(T,LL,LL,3,M)
    for m in 1:M
        X_divrot[:,:,1,m] = X[:,:,1,m]
        divmat, rotmat = get_div_rot(X[:,:,1,m],latt)
        X_divrot[:,:,2,m] = divmat
        X_divrot[:,:,3,m] = rotmat
    end
    return X_divrot
end

## Provide muss
function provide_muss(X::Array{T,2}) where T<:AbstractFloat
    LL = size(X,1)
    latt = TriangularLattice(LL,periodic=false)
    X_divrot = zeros(T,LL,LL,2)
    X_divrot[:,:,1] = X
    range = 2:LL-1
    window = Int((LL-1)/2)
    tmp = zeros(T,LL,LL)
    for j in range, i in range
        tmp[i,j] = X[i,j] - 1/2*atan( (i-window) ,(j-window))
    end
    X_divrot[:,:,2] = tmp
    return X_divrot
end

function provide_muss(X::Array{T,3}) where T<:AbstractFloat
    LL = size(X,1)
    M  = size(X,3)
    X_divrot = zeros(T,LL,LL,2,M)
    X_divrot[:,:,1,:] = X
    tmp = zeros(T,LL,LL)
    range = 2:LL-1
    window = Int((LL-1)/2)
    for m in 1:M
        for j in range, i in range
            tmp[i,j] = X[i,j,m] - 1/2*atan( (i-window) ,(j-window))
        end
        X_divrot[:,:,2,m] = tmp
    end
    return X_divrot
end

function provide_muss(X::Array{T,4}) where T<:AbstractFloat
    LL = size(X,1)
    M  = size(X,4)
    X_divrot = zeros(T,LL,LL,2,M)
    X_divrot[:,:,1,:] = X
    tmp = zeros(T,LL,LL)
    range = 2:LL-1
    window = Int((LL-1)/2)
    for m in 1:M
        for j in range, i in range
            tmp[i,j] = X[i,j,1,m] - 1/2*atan( (i-window) ,(j-window))
        end
        X_divrot[:,:,2,m] = tmp
    end
    return X_divrot
end

## Provide Divergence and Rotationnal
function provide_div_rot_muss(X::Array{T,2}) where T<:AbstractFloat
    LL = size(X,1)
    latt = TriangularLattice(LL,periodic=false)
    X_divrot = zeros(T,LL,LL,4)
    X_divrot[:,:,1] = X
    divmat, rotmat = get_div_rot(X,latt)
    X_divrot[:,:,2] = divmat
    X_divrot[:,:,3] = rotmat
    range = 2:LL-1
    window = Int((LL-1)/2)
    tmp = zeros(T,LL,LL)
    for j in range, i in range
        tmp[i,j] = X[i,j,1,m] - 1/2*atan( (i-window) ,(j-window))
    end
    X_divrot[:,:,4] = tmp

    return X_divrot
end

function provide_div_rot_muss(X::Array{T,3}) where T<:AbstractFloat
    LL = size(X,1)
    M  = size(X,3)
    latt = TriangularLattice(LL,periodic=false)
    X_divrot = zeros(T,LL,LL,4,M)
    range = 2:LL-1
    window = Int((LL-1)/2)
    for m in 1:M
        X_divrot[:,:,1,m] = X[:,:,m]
        divmat, rotmat = get_div_rot(X[:,:,m],latt)
        X_divrot[:,:,2,m] = divmat
        X_divrot[:,:,3,m] = rotmat
        tmp = zeros(T,LL,LL)
        for j in range, i in range
            tmp[i,j] = X[i,j,m] - 1/2*atan( (i-window) ,(j-window))
        end
        X_divrot[:,:,4,m] = tmp
    end
    return X_divrot
end
function provide_div_rot_muss(X::Array{T,4}) where T<:AbstractFloat
    LL = size(X,1)
    M  = size(X,4)
    range = 2:LL-1
    window = Int((LL-1)/2)
    latt = TriangularLattice(LL,periodic=false)
    X_divrot = zeros(T,LL,LL,4,M)
    for m in 1:M
        X_divrot[:,:,1,m] = X[:,:,1,m]
        divmat, rotmat = get_div_rot(X[:,:,1,m],latt)
        X_divrot[:,:,2,m] = divmat
        X_divrot[:,:,3,m] = rotmat
        tmp = zeros(T,LL,LL)
        for j in range, i in range
            tmp[i,j] = X[i,j,1,m] - 1/2*atan( (i-window) ,(j-window))
        end
        X_divrot[:,:,4,m] = tmp
    end
    return X_divrot
end
