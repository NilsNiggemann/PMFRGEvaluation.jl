module PMFRGEvaluation
using Dierckx,SmoothingSplines,DelimitedFiles,EllipsisNotation,RecursiveArrayTools,Reexport, Roots,Plots

@reexport using HDF5,SpinFRGLattices,FRGLatticePlotting,StaticArrays,LaTeXStrings,HDF5Helpers

export  PMResults, Thermoplots,  cutData, cutDataAndRecompute
Base.@kwdef struct PMResults
    T::Vector{Float64}
    Chi_TR::Array{Float64,2}
    gamma_TxN::Array{Float64,3}
    fint::Vector{Float64}
    N::Int
    NLen::Int
    NUnique::Int
    f::Vector{Float64}
    e::Vector{Float64}
    c::Vector{Float64}
    s::Vector{Float64}
end

include("FileReading.jl")
export ReadPMResults, GetThermo,readLastGroupElements,getMaxVertexFlow,getMaxChiTR,getChiTRnu,h5keys,getCorr,getNorms

include("Thermodynamics.jl")
export deriv, get_e, get_c, get_e, get_c, reverseTOrder,getHTSE,HTSE_keys

include("FiniteSizeScaling.jl")
export getChiIntPol,getCrossingPoint

include("PlotFunctions.jl")
export plotgamma_T, plotgamma,plotMaxVertexFlow,VertexRplot

include("ConsistencyCheck.jl")
export getChilocal

include("EqualTimeCorrelator.jl")
export equalTimeChiBeta,EnergyBeta,get_e_Chi

"""Removes T points and re-computes the derivative"""
function cutData(Results,index1,index2 = 0)
    (;T,fint,f,e,s,c,Chi_TR,N,NLen,NUnique,gamma_TxN) = Results
    fields = deepcopy((T,fint,f,e,s,c,Chi_TR,N,NLen,NUnique,gamma_TxN))
    Names = (:T,:fint,:f,:e,:s,:c,:Chi_TR,:N,:NLen,:NUnique,:gamma_TxN)
    init = Dict(Name => field for (Name,field) in zip(Names,fields))

    slice(x) = x[begin+index1:end-index2,fill(:,ndims(x)-1)...] #slices array along first dim (Temperature)
    for (key,val) in init
        if val isa AbstractArray
            init[key] = slice(val)
        end
    end
    return PMResults(;init...)
end

"""Removes T points and re-computes the derivative"""
function cutDataAndRecompute(Results,removeinds::Vector;kwargs...)
    (;T,fint,Chi_TR,N,NLen,NUnique,gamma_TxN) = Results
    fields = Dict(:T =>T,:fint_T => fint,:Chi_TR => Chi_TR,:N => N,:NLen => NLen,:NUnique => NUnique,:gamma_TxN => gamma_TxN)
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