
function deriv(y::AbstractArray,x,order=1) 
    func = Spline1D(x, y, k=3,bc="extrapolate") 
    # func = Spline1D(x, y; w=ones(length(x)), k=3, bc="nearest", s=0.0)
    derivative(func,x,nu=order)
end

function get_e(f,T)
    return(f(T)-T*derivative(f,T))
end

function get_c(f,T)
    return(-T*derivative(f,T,nu = 2))
end

function get_e(f::AbstractArray,T::AbstractArray)
    return( f .-T .*deriv(f,T))
end

function get_c(f::AbstractArray,T::AbstractArray)
    return( -T .* deriv(f,T,2))
end

function GetThermo(PMData;skipvals = 1,smoothen = false,smoothParam = 0.001,SplineDegree = 3)

   (;T,fint_T,Chi_TR,gamma_TxN,N,NLen,NUnique,Chi_TRnu) = PMData

   T = T[1:skipvals:end]
   fint_T = fint_T[1:skipvals:end]
   Chi_TR = Chi_TR[1:skipvals:end,:]
   Chi_TRnu = Chi_TRnu[1:skipvals:end,:,:]
   gamma_TxN = gamma_TxN[1:skipvals:end,:,:]
   
   f_T = similar(fint_T)
   e_T = similar(fint_T)
   c_T = similar(fint_T)
   s_T = similar(fint_T)
   f = similar(fint_T)

   Res = (;N, NLen, NUnique, T ,Chi_TR, Chi_TRnu, gamma_TxN , fint = fint_T,f,e=e_T,c=c_T,s=s_T)
   
   if length(T) < SplineDegree
        return Res
    end
    # smooth Data for f:
    if smoothen
        spl = fit(SmoothingSpline, T, fint_T, smoothParam) # smoothing parameter low means less smoothing
        fint_T = SmoothingSplines.predict(spl,T) # fitted vector
    end
    f = -T*log(2) +fint_T
    # f_intPol = intpol(fint_T,T)
    f_intPol = Spline1D(T, f, k=SplineDegree,bc="extrapolate") 
    # return f_intPol
    for (iT,Temp) in enumerate(T)
        e_T[iT] = get_e(f_intPol,Temp)
        c_T[iT] = get_c(f_intPol,Temp)
        s_T[iT] = (e_T[iT]-f[iT])/Temp
    end
    return Res
end

function GetThermo(Filename;kwargs...)
    selecter = (endOfFirstDim,endOfLastDim)[getLambdaDim(Filename)]
    res = ReadPMResults(Filename,selecter)
    return GetThermo(res;kwargs...)
end

function reverseTOrder(T,Arrays...)
    neworder = length(T):-1:1
    reorder(x) = x[neworder,fill(:,ndims(x)-1)...]
    return reorder.((T,Arrays...))
end

function getLambdaDim(Filename)
    f_ints = readGroupElements(Filename,"f_int")
    LambdaLs = length.(readGroupElements(Filename,"Lambda"))
    NUniques = readGroupElements(Filename,"NUnique")
    dims = size.(f_ints)
    if all(getindex.(dims,1) == LambdaLs) || all(getindex.(dims,2) == NUniques)
        return 1
    elseif all(getindex.(dims,1) == NUniques) || all(getindex.(dims,2) == LambdaLs)
        return 2
    end
    @warn "Could not read Lambda convention from file"
end

function PMResults(Filename::AbstractString;kwargs...)
    h5open(Filename) do f
        return PMResults(f;kwargs...)
    end
end

function PMResults(Filename::AbstractString,key::AbstractString;kwargs...)
    h5open(Filename) do f
        return PMResults(f[key];kwargs...)
    end
end

function PMResults(f::Union{HDF5.File,HDF5.Group};kwargs...)
    selecter = (endOfFirstDim,endOfLastDim)[getLambdaDim(f)]
    res = ReadPMResults(f,selecter)
    return GetThermo(res;kwargs...)
end