w(n,T) =pi*T*(2*n+1.)

function gamma_(gamma::AbstractArray, nw::Integer)
    Ngamma = size(gamma,1)
    s = 1
    if nw<0
        nw = -nw -1
        s = -1
    end
    iw = get_sign_iw(nw,Ngamma)
    return s*gamma[iw]
end

function get_sign_iw(nw::Integer,N::Integer)
    s = sign(nw)
    nw_bounds = min( abs(nw), N-1)
    return s*nw_bounds+1
end

function g(n::Integer,T,gamma::AbstractVector)
    wn= w(n,T)
    return 1/(wn+gamma_(gamma,n))
end

function getChilocal(T,gamma::AbstractVector,N = length(gamma))
    res = 0. # where does the factor 2 come from?
    for n in -N:N-1
        wn = w(n,T)
        res += T/wn * g(n,T,gamma)
    end
    return res
end
