function getMaxVertexFlow(Filename,index,RDim = 1)
    h5open(Filename,"r") do f
        key = keys(f)[index]
        Lambda = Array(f[key*"/Lambda"])
        T = Array(f[key*"/T"])
        MaxVa = maximum( Array(f[key*"/MaxVa"]), dims = RDim) |> vec
        MaxVb = maximum( Array(f[key*"/MaxVb"]), dims = RDim) |> vec
        MaxVc = maximum( Array(f[key*"/MaxVc"]), dims = RDim) |> vec
        return (T=T,Lambda = Lambda,MaxVa=MaxVa,MaxVb=MaxVb,MaxVc=MaxVc )
    end
end

"""returns susceptibility at the value of Lambda where chi(R,Lam) is largest. Assumes that this point will be the position of a peak of the full susceptibility"""
function getMaxChiTR(Filename)
    Chi_RLam_T = readGroupElements(Filename,"Chi")
    TLen = length(Chi_RLam_T)
    NPairs = size(first(Chi_RLam_T),1)
    Chi_TR = Matrix{Float64}(undef,TLen,NPairs)
    for (iT,chiRL) in enumerate(Chi_RLam_T)
        maxR,MaxLam = Tuple(argmax(chiRL)) #gets position of global(!) maximum
        # println(length.((Chi_TR[iT,:], chiRL[:,MaxLam])))
        Chi_TR[iT,:] .= @view chiRL[:,MaxLam]
    end
    return Chi_TR
end

function getChiTRnu(File::HDF5.H5DataStore)
    T = readGroupElements(File,"T")
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

function getChiTRnu(Filename,key = "/")
    h5open(Filename,"r") do File1
        File = File1[key]
        T = readGroupElements(File,"T")
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

areParallel(v1::AbstractVector,v2::AbstractVector) =  isapprox(norm(v1)*norm(v2),abs(v1' * v2),atol = 1E-14)

function getNorms(Lattice)
    (;Basis,PairList,PairTypes) = Lattice
    (;refSites) = Basis
    norm(i) = dist(refSites[PairTypes[i].xi],PairList[i],Basis)
    
    norms = norm.(eachindex(PairList))
    return norms
end
function getCorr(key,Filename,index,Lattice)
    groups = h5keys(Filename)
    (;Basis,PairList,PairTypes) = Lattice
    (;refSites) = Basis
    # R1 = refSites[1]
    # norm(R) = dist(R1,R,Basis)
    norm(i) = dist(refSites[PairTypes[i].xi],PairList[i],Basis)
    
    norms = norm.(eachindex(PairList))
    corr = abs.(h5read(Filename,string(groups[index],"/",key))[:,end])
    return norms,corr
end
function getCorr(Direction,key,Filename,index,Lattice)
    (;Basis,SiteList,PairList,PairTypes,pairToInequiv) = Lattice
    cartDirection = Direction' *inv(Basis.T) |> vec
    inDirection(x) = areParallel(cartDirection,x)
    # println(cartDirection)
    # println(getCartesian.(SiteList,Ref(Basis)))
    SiteInds = findall(inDirection,getCartesian.(SiteList,Ref(Basis)))
    Rvecs = SiteList[SiteInds]
    Ineq_Rvecs = unique(last.(pairToInequiv.(Ref(Basis.refSites[1]),Rvecs)))
    println(getCartesian.(Ineq_Rvecs,Ref(Basis)))
    norms,corr = getCorr(key,Filename,index,Lattice)
    indices = findall(x-> x in Ineq_Rvecs,PairList)
    # println(Rvecs,Ineq_Rvecs,indices)
    return norms[indices],corr[indices]
end

function getCorr_t0(Filename,index,Lattice)
    groups = h5keys(Filename)
    (;Basis,PairList,PairTypes) = Lattice
    (;refSites) = Basis
    norm(i) = dist(refSites[PairTypes[i].xi],PairList[i],Basis)
    
    norms = norm.(eachindex(PairList))
    chiRnu = h5read(Filename,string(groups[index],"/","Chi_nu"))
    Chit0 = equalTimeChi(chiRnu)
    corr = abs.(Chit0)
    return norms,corr
end

function endOfFirstDim(Arr::AbstractArray)
    dims = size(Arr)
    selectdim(Arr,1,dims[begin])
end
endOfFirstDim(A::HDF5.Dataset) = endOfFirstDim(Array(A))

function endOfLastDim(Arr::AbstractArray)
    dims = size(Arr)
    d = length(dims)
    selectdim(Arr,d,dims[end])
end
endOfLastDim(A::HDF5.Dataset) = endOfLastDim(Array(A))

getOnly(Filename,key) =only(unique(readGroupElements(Filename,key)))

function ReadPMResults(File,selecter=endOfLastDim)
    T = readGroupElements(File,"T")
    Tlen = length(T)
    # sort values to ascending order in T
    
    N = getOnly(File,"N")
    NLen = getOnly(File,"NLen")
    NUnique = getOnly(File,"NUnique")
    #allocate correct memory
    k1 = first(keys(File))
    Chidims = size( selecter(File[k1*"/Chi"]))
    gammadims = size( selecter(File[k1*"/gamma"]))

    fint_T = zeros(Tlen)
    Chi_TR = zeros(Tlen,Chidims...)
    gamma_TxN = zeros(Tlen,gammadims...)
    #write to arrays
    for (i,key) in enumerate(keys(File))
        fint_T[i] =  mean(selecter(File[key*"/f_int"])) # read fint for last value of Lambda
        Chi_TR[i,:] .= selecter(File[key*"/Chi"])
        gamma_TxN[i,:,:] .= selecter(File[key*"/gamma"])
    end
    #sort according to T
    keylist = sortperm(T)
    T = T[keylist]
    
    fint_T = fint_T[keylist]
    Chi_TR = Chi_TR[keylist,:]
    gamma_TxN = gamma_TxN[keylist,:,:]
    Chi_TRnu = getChiTRnu(File)[keylist,:,:]
    return (;N, NLen, NUnique, T ,fint_T ,Chi_TR, Chi_TRnu, gamma_TxN)
end

function ReadPMResults(Filename::String,selecter=endOfLastDim)
    h5open(Filename,"r") do File
        ReadPMResults(File,selecter)
    end
end

function getNumberFromName(Name,subName)
    res_string = split(Name,subName)[end]
    for i in length(res_string):-1:1
        N = tryparse(Int,res_string[1:i])
        if N !== nothing
            return N
        end
    end
    error("Could not get ", subName, "from string ",Name)
end


"""Read all results from file and wrap them in a named tuple for convenience. Also creates interpolations for thermodynamic observables and fourier transforms of the susceptibility."""
function AllPMResults(filename,Lattice::AbstractLattice)
    allkeys = h5keys(filename,1)
    T = readGroupElements(filename,"T")
    order = sortperm(T)
    T .= T[order]
    function read(key)
        value = readGroupElements(filename,key)
        allequal(value) && return only(unique!(Array(value)))
        return value[order]
    end
    allRes = Dict([Symbol(k)=>read(k) for k in allkeys])

    fint = mean.(last.(allRes[:f_int]))
    System = Lattice.System
    M = 
    if haskey(allRes,:Name) && occursin("2S_",Name) 
        Name = allRes[:Name]
        getNumberFromName(Name,"2S_")
        else
            1
    end
    Spin = M/2
    
    Thermos = getThermoIntPol(allRes[:T],fint,M)

    Chi_TR_12 = Array.(endOfLastDim.(allRes[:Chi]) )
    Chi_TR = convertSusceptibilityToSpinS.(Chi_TR_12,Ref(Spin),Ref(System.OnsitePairs))

    @assert length(Chi_TR) == length(allRes[:T]) "Could not find a susceptibility for every temperature"
    u = only(unique!(length.(Chi_TR)))
    @assert (u == System.Npairs) "number of pairs does not match: $u vs. $(System.Npairs) This means that you are likely trying to use the wrong geometry for your data!"
    @assert allRes[:NUnique] == System.NUnique "NUnique incompatible with geometry"

    Chi_Tq = getFourier.(Chi_TR,Ref(Lattice))
    ChiT_Rnu = allRes[:Chi_nu]

    Sij_TR = equalTimeChiBeta.(ChiT_Rnu) .*T
    Sij_TR = convertSusceptibilityToSpinS!.(Sij_TR,Ref(Spin),Ref(System.OnsitePairs))
    S_Tq = getFourier.(Sij_TR,Ref(Lattice))
    
    gamma = endOfLastDim.(allRes[:gamma])
    @inline Δ(x=1,M_sum = 1000) = [wardIdentityviolation(Ti,gammai[x,:],Chi_R[System.OnsitePairs[x]],M_sum) for (Ti,gammai,Chi_R) in zip(T,gamma,Chi_TR_12)]
    delete!(allRes,(:Chi,:Npairs,:NUnique,))
    return (;allRes..., Thermos...,Spin,Lattice,Chi_TR,Chi_Tq,Δ,Sij_TR,S_Tq)
end

function AllPMResults(filename::AbstractString,LatticeGenerator::Function)
    NLen = getOnly(filename,"NLen")
    Lattice = LatticeGenerator(NLen)
    AllPMResults(filename,Lattice)
end

function AllPMResults(filenames::AbstractArray,LatticeGenerator::Function)
    NLens = [getOnly(f,"NLen") for f in filenames]
    Lattices = Dict([NLen => LatticeGenerator(NLen) for NLen in unique(NLens)])
    [AllPMResults(f,Lattices[NLen]) for (NLen,f) in zip(NLens,filenames)]
end

function AllPMResults(filenames::AbstractArray,geometryGenerator::Function,Module::Function)
    LatticeGenerator(NLen) = LatticeInfo(geometryGenerator(NLen),Module)
    AllPMResults(filenames,LatticeGenerator)
end