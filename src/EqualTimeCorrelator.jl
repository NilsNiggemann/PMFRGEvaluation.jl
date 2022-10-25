
"""Computes χ_ij(τ=0)/T"""
function equalTimeChiBeta(Chi_RNu::AbstractArray,N_nu = size(Chi_RNu,2))
    # Npairs,N_nu = size(Chi_RNu)
    Chi_Tau0 = Chi_RNu[:,begin] # add static nu=0 component, appears only once in sum
    for R in eachindex(Chi_Tau0)
        sum = 0.
        for n in 2:N_nu # dynamic components
            sum += Chi_RNu[R,n]
        end
        Chi_Tau0[R] += 2*sum# Chi(nu=0)
    end
    return Chi_Tau0
end

"""Compute energy divided by temperature from spin correlations"""
function EnergyBeta(Chi_RNu, Lattice,Nnu)
    (;PairList,SiteList,PairTypes,Basis,UnitCell,pairToInequiv,Npairs )= Lattice
    J_ij = Lattice.System.couplings
    E = 0.
    Chi_Tau0 = equalTimeChiBeta(Chi_RNu,Nnu)

    for i_site in UnitCell
        # println(Ri)
        for j_site in SiteList # site summation
            R_Ref,ij = pairToInequiv(i_site,j_site) #Map j to correct pair so that we may use Chi_0,j'
            xi = getSiteType(R_Ref,Basis)
            pair = MapToPair(xi,ij,PairList,PairTypes)
            if pair !== 0
                E += 3/(2*Basis.NCell) *J_ij[pair] * Chi_Tau0[pair]
            end
            # println(j_site,Chi_R[pair])
        end
    end
    return E
end

function get_e_Chi(Chi_TRnu,Trange,Lattice,Nnu)
    e_Chi = similar(Trange)
    for (iT,T) in enumerate(Trange)
        e_Chi[iT] = @views T*EnergyBeta(Chi_TRnu[iT,:,:],Lattice,Nnu)
    end
    return e_Chi
end

function equalTimeChiBeta(f::Union{HDF5.Group,HDF5.File},key ="",args...)
    chiRnu = f[key*"/Chi_nu"] |> read

    equalTimeChiBeta(chiRnu,args...)
end

function equalTimeChiBeta(Filename::AbstractString,key = "",args...)
    h5open(Filename) do f
        equalTimeChiBeta(f,key,args...)
    end
end
