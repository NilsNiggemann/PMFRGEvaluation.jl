function plotMaxVertexFlow(Filename,index,pl = plot(),RDim = 1,kwargs...)
    T,Lambda,MaxVa,MaxVb,MaxVc = getMaxVertexFlow(Filename,index,RDim)
    Lambda 
    plot!(Lambda,MaxVa,xlims = xlims,label = L"max(\Gamma_a)",title = "\$T = $T \$",xlabel = L"\Lambda")
    plot!(Lambda,MaxVb,label = L"max(\Gamma_b)")
    plot!(Lambda,MaxVc,label = L"max(\Gamma_c)")
    plot!(;kwargs...)
    return pl
end
function VertexRplot!(pl,Vertex::AbstractVector,Lattice;kwargs...)
    (;Basis,PairList,PairTypes) = Lattice
    (;refSites) = Basis
    norm(i) = dist(refSites[PairTypes[i].xi],PairList[i],Basis)
    norms = norm.(eachindex(PairList))
    scatter!(norms,Vertex;kwargs...)
end

VertexRplot!(Vertex::AbstractVector,Lattice;kwargs...) = VertexRplot!(current(),Vertex,Lattice;kwargs...)

VertexRplot(Vertex::AbstractVector,Lattice;kwargs...) = VertexRplot!(plot(),Vertex,Lattice;kwargs...)

function VertexRplot(Filename::String,index,Lattice;RDim = 2,kwargs...)
    T = readGroupElements(Filename,"T")[index]
    key = h5keys(Filename)[index]
    selecter = (endOfFirstDim,endOfLastDim)[RDim]
    MaxVa = h5read(Filename,key*"/MaxVa") |>selecter |> vec
    MaxVb = h5read(Filename,key*"/MaxVb") |>selecter |> vec
    MaxVc = h5read(Filename,key*"/MaxVc") |>selecter |> vec
    # println(MaxVc)
    VertexRplot(MaxVa,Lattice,label = stringLatex(L"max_{s,t,u}(\Gamma_a)"),xlabel = L"R/a",title = "\$T = $T \$")
    VertexRplot!(MaxVb,Lattice,label = stringLatex(L"max_{s,t,u}(\Gamma_b)"))
    VertexRplot!(MaxVc,Lattice,label = stringLatex(L"max_{s,t,u}(\Gamma_c)"))
    plot!(;kwargs...)
end
function plotgamma_T(Results,iT,pl = plot())
    (;N,gamma_TxN) = Results
    for x in 1:NUnique
        scatter!(pl,1:N,gamma_TxN[iT,x,:],ylabel = L"\gamma",xlabel = L"N",label = nothing)
    end
end

function plotgamma(Results,x=1,Nmax=size(Results.gamma_TxN,3))
    (;T,gamma_TxN) = Results
    surface(1:Nmax,T,gamma_TxN[:,x,1:Nmax],zlabel = L"\gamma",ylabel=L"T",xlabel = L"N",label = nothing,c= colorscheme)
end

