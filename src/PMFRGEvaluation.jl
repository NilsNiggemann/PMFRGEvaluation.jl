module PMFRGEvaluation
using Dierckx,SmoothingSplines,RecursiveArrayTools,Reexport, Roots

@reexport using HDF5,SpinFRGLattices,FRGLatticeEvaluation,StaticArrays,LaTeXStrings,HDF5Helpers

export  PMResults,  cutData, cutDataAndRecompute

include("FileReading.jl")
export ReadPMResults, GetThermo,readLastGroupElements,getMaxVertexFlow,getMaxChiTR,getChiTRnu,h5keys,getCorr,getNorms,AllPMResults

include("Thermodynamics.jl")
export deriv, get_e, get_c, get_e, get_c, reverseTOrder,getHTSE,HTSE_keys

include("FiniteSizeScaling.jl")
export getChiIntPol,getCrossingPoint

include("ConsistencyCheck.jl")
export getChilocal,wardIdentityviolation

include("EqualTimeCorrelator.jl")
export equalTimeChiBeta,EnergyBeta,get_e_Chi

"""Removes T points and re-computes the derivative"""
function cutData(Results,index1,index2 = 0)
    (;T,fint,f,e,s,c,Chi_TR,N,NLen,NUnique,gamma_TxN) = Results
    fields = deepcopy((T,fint,f,e,s,c,Chi_TR,N,NLen,NUnique,gamma_TxN))
    Names = (:T,:fint,:f,:e,:s,:c,:Chi_TR,:N,:NLen,:NUnique,:gamma_TxN)
    res = (Name = field for (Name,field) in zip(Names,fields))

    slice(x) = x[begin+index1:end-index2,fill(:,ndims(x)-1)...] #slices array along first dim (Temperature)
    for (key,val) in res
        if val isa AbstractArray
            res[key] = slice(val)
        end
    end
    return res
end

"""Removes T points and re-computes the derivative"""
function cutDataAndRecompute(Results,removeinds::Vector;kwargs...)
    (;T,fint,Chi_TR,N,NLen,NUnique,gamma_TxN) = Results
    fields = (;T, fint, Chi_TR, N, NLen, NUnique, gamma_TxN)
    inds = deleteat!(collect(eachindex(T)),removeinds)
    slice(x) = x[inds,fill(:,ndims(x)-1)...] #slices array along first dim (Temperature)
    for (key,val) in fields
        if val isa AbstractArray
            fields[key] = slice(val)
        end
    end
    return GetThermo(fields;kwargs...)
end

end#module