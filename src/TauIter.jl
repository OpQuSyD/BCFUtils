export TauExpIter, get_TauExpIter_for_Ohmic

abstract type AbstractTauIter{T <: Real} end

"""
Distribution of time arguments used to calculate the difference
between a given BCF and an exponential approximation.

The distribution is given by the expression

    τ(u) = τ_0 (exp(a u) - 1)

where the values `u` are sampled uniformly from the closed interval [0, 1].

To construct such a range use the function `get_TauExpIter(t_max, s, a_tilde, u_tilde, len)`
which calculates the raw parameters `τ_0` and `a` from more meaning full values.

`t_max`: specifies the time interval `τ ∈ [0, t_max]` of interest
`s`: power of the Ohmic spectral density (use `BCF.red_ohmic_bcf` as reference bcf, i.e., (1+ i τ)^(-s-1).
`a_tilde`, `u_tilde`: 
    Determines the distribution such that a_tilde = |α(τ(u_tilde))| holds (approximately).
    The default values a_tilde=0.1 and u_tilde=0.5 mean that half of the non-linear tau range
    make up the the region where the bcf (abs. value) drops from 1 to 0.1.
"""
struct TauExpIter{T<:AbstractFloat, L<:Integer} <: AbstractTauIter{T}
    a::T
    tau0::T
    len::L
end

function calc_tau_0_a(t_max, s, a_tilde, u_tilde)
    t_tilde = sqrt(a_tilde^(-2/(s+1)) - 1)
    c = (1/t_tilde - u_tilde/t_max) / (u_tilde - 1)
    b = (u_tilde * log(t_max) - log(t_tilde)) / (u_tilde - 1)
    tau_0 = exp(b)/(1 + c*exp(b))
    a = log((t_max + tau_0) / tau_0)
    return tau_0, a
end


"""
convenient function to construct TauExpIter from the meaningful parameters

ToDo: explain parameters
"""
function get_TauExpIter_for_Ohmic(t_max, s, len, a_tilde=0.1, u_tilde=0.5)
    0 < a_tilde < 1 || thorw(ArgumentError("0 < a_tilde < 1 needs to hold"))
    0 < u_tilde < 1 || thorw(ArgumentError("0 < u_tilde < 1 needs to hold"))
    
    tau_0, a = promote(calc_tau_0_a(t_max, s, a_tilde, u_tilde)...)
    TauExpIter(a, tau_0, len)
end


function Base.iterate(::TauExpIter{T, L}) where {T, L}
    return zero(T), one(L)
end

function Base.iterate(t::TauExpIter{T, L}, state::L) where {T, L}
    state < t.len || return nothing
    return t.tau0*(exp(t.a*state/(t.len-1))-1), state+1
end

function Base.length(t::TauExpIter)
    return t.len
end

