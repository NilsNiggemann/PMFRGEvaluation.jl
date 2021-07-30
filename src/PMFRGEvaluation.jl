module PMFRGEvaluation
using Dierckx,SmoothingSplines,DelimitedFiles,EllipsisNotation,RecursiveArrayTools,Reexport

@reexport using Parameters,HDF5,SpinFRGLattices,FRGLatticePlotting,StaticArrays,Plots,LaTeXStrings

export deriv, get_e, get_c, get_e, get_c, getNumberFromName, ReadPMResults, GetThermo, reverseTOrder, PMResults, Thermoplots, plotgamma_T, plotgamma, getHTSE, cutData, cutDataAndRecompute, HTSE_keys,readGroupElements,readLastGroupElements,plotMaxVertexFlow,stringLatex,VertexRplot,getChiTRnu,h5keys,getCorr

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

@with_kw struct PMResults
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

function getNumberFromName(Name,subName)
    res_string = split(Name,subName*"=")[end]
    for i in length(res_string):-1:1
        N = tryparse(Int,res_string[1:i])
        if N !== nothing
            return N
        end
    end
    error("Could not get ", subName, "from string ",Name)
end

function ReadPMResults_old(Filename)

    T = h5read(Filename,"Trange")
    # sort values to ascending order in T
    keylist = sortperm(T)
    T = T[keylist]

    fint_T = mean(h5read(Filename,"fint_Tx")[keylist,:],dims=2)[:,1]
    Chi_TR = h5read(Filename,"Chi_TR")[keylist,:]
    gamma_TxN = h5read(Filename,"gamma_TxN")[keylist,:,:]
    N = h5read(Filename,"N")
    NLen= 0
    try
        NLen = h5read(Filename,"NLen")
    catch
        NLen = getNumberFromName(Filename,"NLen")
    end
    NUnique = h5read(Filename,"NUnique")
    # skip values 
   
    return (Dict(:T => T ,:fint_T => fint_T ,:Chi_TR => Chi_TR ,:gamma_TxN => gamma_TxN ,:N => N ,:NLen => NLen ,:NUnique => NUnique))
end

"""Fetches key from file for each group and appends results to a list"""
function readGroupElements(File,key)
    h5open(File,"r") do f
        return [Array(f[string(Group,"/",key)]) for Group in keys(f)]
    end
end

function readLastGroupElements(File,key)
    h5open(File,"r") do f
        return VectorOfArray([ Array(f[string(Group,"/",key)])[end,..] for Group in keys(f)]) #using EllipsisNotation to get index in first dimension
    end
end

function h5keys(Filename::String)
    h5open(Filename,"r") do f
        return keys(f)
    end
end

function plotMaxVertexFlow(Filename,index,xlims = (0.,2.5);kwargs...)
    h5open(Filename,"r") do f
        key = keys(f)[index]
        Lambda = Array(f[key*"/Lambda"])
        T = Array(f[key*"/T"])
        MaxVa = maximum( Array(f[key*"/MaxVa"]), dims = 2)
        MaxVb = maximum( Array(f[key*"/MaxVb"]), dims = 2)
        MaxVc = maximum( Array(f[key*"/MaxVc"]), dims = 2)
        pl = plot(Lambda,MaxVa,xlims = xlims,label = L"max(\Gamma_a)",title = "\$T = $T \$",xlabel = L"\Lambda")
        plot!(Lambda,MaxVb,label = L"max(\Gamma_b)")
        plot!(Lambda,MaxVc,label = L"max(\Gamma_c)")
        plot!(;kwargs...)
        return pl
    end
end

function getChiTRnu(Filename)
    h5open(Filename,"r") do File
        T = readGroupElements(Filename,"T")
        keylist = sortperm(T)
        Tlen = length(T)
        k1 = first(keys(File))
        Chi_nuDims = size( File[k1*"/Chi_nu"][:,:])

        Chi_TRnu = zeros(Tlen,Chi_nuDims...)
        #write to arrays
        for (i,key) in enumerate(keys(File))
            Chi_TRnu[i,:,:] .= File[key*"/Chi_nu"][:,:]
        end
        return Chi_TRnu[keylist,:,:]
    end
end

function getCorr(key,Filename,index,Lattice)
    groups = h5keys(Filename)
    @unpack Basis,PairList,PairTypes = Lattice
    @unpack refSites = Basis
    # R1 = refSites[1]
    # norm(R) = dist(R1,R,Basis)
    norm(i) = dist(refSites[PairTypes[i].xi],PairList[i],Basis)
    
    norms = norm.(eachindex(PairList))
    corr = abs.(h5read(Filename,string(groups[index],"/",key))[end,:])
    return norms,corr
end

function VertexRplot(Filename,index,Lattice;kwargs...)
    T = readGroupElements(Filename,"T")[index]
    @unpack Basis,PairList,PairTypes = Lattice
    @unpack refSites = Basis
    # R1 = refSites[1]
    # norm(R) = dist(R1,R,Basis)
    norm(i) = dist(refSites[PairTypes[i].xi],PairList[i],Basis)
    
    norms = norm.(eachindex(PairList))
    MaxVa = readLastGroupElements(Filename,"MaxVa")
    MaxVb = readLastGroupElements(Filename,"MaxVb")
    MaxVc = readLastGroupElements(Filename,"MaxVc")
    scatter(norms,MaxVa[:,index],label = stringLatex(L"max_{s,t,u}(\Gamma_a)"),xlabel = L"R/a",title = "\$T = $T \$")
    scatter!(norms,MaxVb[:,index],label = stringLatex(L"max_{s,t,u}(\Gamma_b)"))
    scatter!(norms,MaxVc[:,index],label = stringLatex(L"max_{s,t,u}(\Gamma_c)"))
    plot!(;kwargs...)
end
function ReadPMResults(Filename)
    h5open(Filename,"r") do File
        T = readGroupElements(Filename,"T")
        Tlen = length(T)
        # sort values to ascending order in T
        
        N = only(unique(readGroupElements(Filename,"N")))
        NLen= 0
        try
            NLen = only(unique(readGroupElements(Filename,"NLen")))
        catch
            NLen = getNumberFromName(Filename,"NLen")
        end
        NUnique = only(unique(readGroupElements(Filename,"NUnique")))
        #allocate correct memory
        k1 = first(keys(File))
        Chidims = size( File[k1*"/Chi"][end,:])
        gammadims = size( File[k1*"/gamma"][end,:,:])

        fint_T = zeros(Tlen)
        Chi_TR = zeros(Tlen,Chidims...)
        gamma_TxN = zeros(Tlen,gammadims...)
        #write to arrays
        for (i,key) in enumerate(keys(File))
            fint_T[i] =  mean(File[key*"/f_int"][end,:]) # read fint for last value of Lambda
            Chi_TR[i,:] .= File[key*"/Chi"][end,:]
            gamma_TxN[i,:,:] .= File[key*"/gamma"][end,:,:]
        end
        #sort according to T
        keylist = sortperm(T)
        T = T[keylist]
        
        fint_T = fint_T[keylist]
        Chi_TR = Chi_TR[keylist,:]
        gamma_TxN = gamma_TxN[keylist,:,:]
        return Dict(:T => T ,:fint_T => fint_T ,:Chi_TR => Array(Chi_TR) ,:gamma_TxN => Array(gamma_TxN) ,:N => N ,:NLen => NLen ,:NUnique => NUnique)
    end
end

function GetThermo(PMData::Dict;skipvals = 1,smoothen = false,smoothParam = 0.001)

    T,fint_T,Chi_TR,gamma_TxN,N,NLen,NUnique = getindex.(Ref(PMData),(:T,:fint_T,:Chi_TR,:gamma_TxN,:N,:NLen,:NUnique))

    T = T[1:skipvals:end]
    fint_T = fint_T[1:skipvals:end]
    Chi_TR = Chi_TR[1:skipvals:end,:]
    gamma_TxN = gamma_TxN[1:skipvals:end,:,:]

    f_T = similar(fint_T)
    e_T = similar(fint_T)
    c_T = similar(fint_T)
    s_T = similar(fint_T)
    # smooth Data for f:
    if smoothen
        spl = fit(SmoothingSpline, T, fint_T, smoothParam) # smoothing parameter low means less smoothing
        fint_T = SmoothingSplines.predict(spl,T) # fitted vector
    end
    f = -T*log(2) +fint_T
    # f_intPol = intpol(fint_T,T)
    f_intPol = Spline1D(T, f, k=3,bc="extrapolate") 
    # return f_intPol
    for (iT,Temp) in enumerate(T)
        e_T[iT] = get_e(f_intPol,Temp)
        c_T[iT] = get_c(f_intPol,Temp)
        s_T[iT] = (e_T[iT]-f[iT])/Temp
    end
    return PMResults(T=T,N = N, NLen = NLen, NUnique = NUnique, Chi_TR=Chi_TR,gamma_TxN=gamma_TxN,fint=fint_T,f=f,e=e_T,c=c_T,s=s_T)
end

function reverseTOrder(T,Arrays...)
    neworder = length(T):-1:1
    reorder(x) = x[neworder,fill(:,ndims(x)-1)...]
    return reorder.((T,Arrays...))
end

function PMResults(Filename;kwargs...)
    res = ReadPMResults(Filename)
    return GetThermo(res;kwargs...)

end

function Thermoplots(Results,pl =plot(layout = (4,1));xAxis = "T",method = plot!,shape = :circle,kwargs...)
    @unpack T,f,e,s,c = Results
    ThermQuantities = (f,e,s,c)
    # linestyles = (:solid,:dash,:dot,:dashdot)
    # shapes = (:circle,:rect,:diamond,:cross)
    # colors = ("blue","red","black","cyan","magenta","green","pink")
    Labels = [L"f",L"e",L"s",L"c"]
    if :label in keys(kwargs)
        legendLabel = (kwargs[:label],"","","")
    else
        legendLabel = ("","","","")
    end
    x = T
    if xAxis in ("Beta", "beta","1/T")
        x = 1 ./T
        xAxis = "\\beta"
    end

    for (i,(obs,lab)) in enumerate(zip(ThermQuantities,Labels))
        method(pl[i],x,obs,ylabel = lab,xlabel = "",shape=shape,xformatter=_->"",top_margin = -20*Plots.px,xlims =  [0,maximum(x)];kwargs...,label = legendLabel[i])
    end
    # plot!(pl[1],,legend = true;kwargs...)
    plot!(pl[end],[],[],label = "",xlabel = latexstring(xAxis),xformatter=x->x,size = (500,700),left_margin=20*Plots.px)
    return pl
end

function Thermoplots_Makie(Results,fig =Figure(resolution = (500,700));xAxis = "T",shape = :circle,kwargs...)
    @unpack T,f,e,s,c = Results
    ThermQuantities = (f,e,s,c)
    # linestyles = (:solid,:dash,:dot,:dashdot)
    # shapes = (:circle,:rect,:diamond,:cross)
    # colors = ("blue","red","black","cyan","magenta","green","pink")
    Labels = [L"f",L"e",L"s",L"c"]
    x = T
    if xAxis in ("Beta", "beta","1/T")
        x = 1 ./T
        xAxis = "\\beta"
    end
    
    for (i,(obs,lab)) in enumerate(zip(ThermQuantities,Labels))
        ax = Axis(fig[i,1], xlabel = "", ylabel = lab,ylabelsize = 20)
        line = lines!(x,obs,kwargs...)
        points = scatter!(x,obs,lw = 2,shape=shape,xformatter=_->"";kwargs...)
    end
    if :label in keys(kwargs)
        leg = Legend(fig, [[line, points]], [kwargs[:label]])
        fig[4,1] = leg
    end
    # fontsize_theme = Theme(fontsize = 10)
    # set_theme!(fontsize_theme)
    trim!(fig.layout)
    # axislegend()
    return fig
end
function plotgamma_T(Results,iT,pl = plot())
    @unpack N,gamma_TxN = Results
    for x in 1:NUnique
        scatter!(pl,1:N,gamma_TxN[iT,x,:],ylabel = L"\gamma",xlabel = L"N",label = nothing)
    end
end

function plotgamma(Results,x=1,Nmax=size(Results.gamma_TxN,3))
    @unpack T,gamma_TxN = Results
    surface(1:Nmax,T,gamma_TxN[:,x,1:Nmax],zlabel = L"\gamma",ylabel=L"T",xlabel = L"N",label = nothing,c= colorscheme)
end


HTSE_keys = (
    "T",
    "Chi_10",
    "U_9",
    "C/k_10",
    "Chi_Pade_4,6",
    "Chi_Pade_5,5",
    "Chi_Pade_6,4",
    "U_Pade_4,5",
    "U_Pade_5,4",
    "C/k_Pade_4,6",
    "C/k_Pade_5,5",
    "C/k_Pade_6,4",
    )
function getHTSE(FileName)
    HTSE_Data = readdlm(FileName,skipstart = 155)
    HTSE = Dict( key => HTSE_Data[:,i] for (i,key) in enumerate(HTSE_keys))
end

"""Removes T points and re-computes the derivative"""
function cutData(Results,index1,index2 = 0)
    @unpack T,fint,f,e,s,c,Chi_TR,N,NLen,NUnique,gamma_TxN = Results
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
function cutDataAndRecompute(Results,index1,index2 = 0;kwargs...)
    @unpack T,fint,Chi_TR,N,NLen,NUnique,gamma_TxN = Results
    fields = Dict(:T =>T,:fint_Tx => fint,:Chi_TR => Chi_TR,:N => N,:NLen => NLen,:NUnique => NUnique,:gamma_TxN => gamma_TxN)

    slice(x) = x[begin+index1:end-index2,fill(:,ndims(x)-1)...] #slices array along first dim (Temperature)
    for (key,val) in fields
        if val isa AbstractArray
            fields[key] = slice(val)
        end
    end
    return GetThermo(fields;kwargs...)
end
function stringLatex(args...)
    res = ""
    for arg in args
        if typeof(arg) == LaTeXString
            res = string(res,"\\ ",arg.s[2:end-1])
        elseif typeof(arg) == String
            res = string(res,"\\ ","\\textrm{$arg}")
        end
    end
    return string("\$",res,"\$")
end

end#module