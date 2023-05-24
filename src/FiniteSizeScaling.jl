const SO3ETA = 0.035
function rescale(Chi,NLen,eta) 
    Chi ./(NLen^(2-eta))
end

function getChiIntPol(kmax::StaticVector,Chi_TR::AbstractMatrix,T::AbstractVector,NLen::Integer,Lattice;eta= SO3ETA)
    Chik = getFlow(kmax,Chi_TR, T,Lattice)
    Chik_res = rescale(Chik,NLen,eta)
    return Dict(
        :T => T,
        :NLen => NLen,
        :Chi => Chi_TR,
        :Chi_k => Chik,
        :Chi_k_Res => Chik_res,
        :intpol_Res => Spline1D(T,Chik_res),
        :kmax => kmax
    )
end

getChiIntPol(kmax::StaticVector,Res,Lattice,eta= SO3ETA) = getChiIntPol(kmax,Res.Chi_TR,Res.T,Res.NLen,Lattice;eta=eta)

function getChiIntPol(Chi_TR::AbstractMatrix,T::AbstractVector,NLen::Integer,Lattice,RegionFunc::Function;eta= SO3ETA,kwargs...)
    kmax = getkMax(Chi_TR[1,:],Lattice,RegionFunc;kwargs...)
    getChiIntPol(kmax,Chi_TR,T,NLen,Lattice;eta=eta)
end

function getChiIntPol(Res,Lattice,RegionFunc::Function;eta = SO3ETA,kwargs...)
    kmax = getkMax(Res.Chi_TR[1,:],Lattice,RegionFunc;kwargs...)
    getChiIntPol(kmax,Res.Chi_TR,Res.T,Res.NLen,Lattice;eta=eta)
end

function getChiIntPol(Res,Lattice::LatticeInfo{B,R,FT,Dim};eta = SO3ETA,kwargs...) where {B,R,FT,Dim}
    kmax = getkMax(Res.Chi_TR[1,:],Lattice;kwargs...)
    getChiIntPol(SVector{Dim,Float64}(kmax),Res.Chi_TR,Res.T,Res.NLen,Lattice;eta=eta)
end

function getIntersect(f1,f2,guess)
    find_zero(x->f1(x)-f2(x), guess)
end

function getCrossingPoint(guess,args)
    allIntersect = Float64[]
    for I in args
        for J in args
            if I != J
                push!(allIntersect,getIntersect(I,J,guess))
            end
        end
    end
    Cross = mean(allIntersect)
    maxErr = maximum(allIntersect) - minimum(allIntersect)
    stderr = std(allIntersect)
    return (;Cross,maxErr,stderr)
end

getCrossingPoint(guess,args...) = getCrossingPoint(guess,[args...])
