using Roots

function rescale(Chi,NLen,eta) 
    Chi ./(NLen^(2-eta))
end

function getChiIntPol(kmax::StaticArray,Res::PMResults,Lattice,eta= 0.04)
    NLen = Res.NLen
    Chik = getFlow(kmax,Res.Chi_TR, Res.T,Lattice)
    Chik_res = rescale(Chik,NLen,eta)
    return Dict(
        :T => Res.T,
        :NLen => NLen,
        :Chi => Res.Chi_TR,
        :Chi_k => Chik,
        :Chi_k_Res => Chik_res,
        :intpol_Res => Spline1D(Res.T,Chik_res),
        :kmax => kmax
    )
end

function getChiIntPol(Res::PMResults,Lattice,RegionFunc::Function,eta= 0.04;kwargs...)
    kmax = getkMax(Res.Chi_TR[1,:],Lattice,RegionFunc;kwargs...)
    getChiIntPol(kmax,Res,Lattice,eta)
end

function getIntersect(f1,f2,guess)
    find_zero(x->f1(x)-f2(x), guess)
end

function getCrossingPoint(guess,args...)
    allIntersect = Float64[]
    for I in args
        for J in args
            if I != J
                push!(allIntersect,getIntersect(I,J,guess))
            end
        end
    end
    meanInt = mean(allIntersect)
    deltaI = maximum(allIntersect) - minimum(allIntersect)
end
