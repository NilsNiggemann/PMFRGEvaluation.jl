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
    return sum(n -> T/w(n,T) * g(n,T,gamma),-N:N-1)
end

function wardIdentityviolation(T::Real,gamma::AbstractVector,Chi00::Real,args...)
    chigamma = getChilocal(T,gamma,args...) 
    return 1 - Chi00/chigamma
end