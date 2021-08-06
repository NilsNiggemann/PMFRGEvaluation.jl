
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
