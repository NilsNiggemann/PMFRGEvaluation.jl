##
allOccurIn(name,args...) = all((occursin(arg,name) for arg in args))
findNames(names,args...) = findall(x->allOccurIn(x,args...),names)
NameFilter(names,args...) = filter(x->allOccurIn(x,args...),names)
OnlyIndex(names,args...) = only(findNames(names,args...))

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

function ArrayReadGroupElements(File,key)
    Data = readGroupElements(File,key)
    d = size(first(Data)) |> length
    return cat(Data...,dims = d+1)
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

areParallel(v1::AbstractVector,v2::AbstractVector) =  isapprox(norm(v1)*norm(v2),abs(v1' * v2),atol = 1E-14)

function getNorms(Lattice)
    @unpack Basis,PairList,PairTypes = Lattice
    @unpack refSites = Basis
    norm(i) = dist(refSites[PairTypes[i].xi],PairList[i],Basis)
    
    norms = norm.(eachindex(PairList))
    return norms
end
function getCorr(key,Filename,index,Lattice)
    groups = h5keys(Filename)
    @unpack Basis,PairList,PairTypes = Lattice
    @unpack refSites = Basis
    # R1 = refSites[1]
    # norm(R) = dist(R1,R,Basis)
    norm(i) = dist(refSites[PairTypes[i].xi],PairList[i],Basis)
    
    norms = norm.(eachindex(PairList))
    corr = abs.(h5read(Filename,string(groups[index],"/",key))[:,end])
    return norms,corr
end
function getCorr(Direction,key,Filename,index,Lattice)
    @unpack Basis,SiteList,PairList,PairTypes,pairToInequiv = Lattice
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
    @unpack Basis,PairList,PairTypes = Lattice
    @unpack refSites = Basis
    # R1 = refSites[1]
    # norm(R) = dist(R1,R,Basis)
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

function ReadPMResults(Filename,selecter=endOfLastDim)
    h5open(Filename,"r") do File
        T = readGroupElements(Filename,"T")
        Tlen = length(T)
        # sort values to ascending order in T
        
        N = getOnly(Filename,"N")
        NLen= 0
        try
            NLen = getOnly(Filename,"NLen")
        catch
            NLen = getNumberFromName(Filename,"NLen")
        end
        NUnique = getOnly(Filename,"NUnique")
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
        return Dict(:T => T ,:fint_T => fint_T ,:Chi_TR => Array(Chi_TR) ,:gamma_TxN => Array(gamma_TxN) ,:N => N ,:NLen => NLen ,:NUnique => NUnique)
    end
end

