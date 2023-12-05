# Ohmic bath
export OhmicExpCO
# multi exp bath 
export MultiExpBCF, GW_2_u, u_2_GW
# utility functions
export sd, bcf

using SpecialFunctions

Γ = gamma


"""
A structured bath is defined via its spectral density J(ω).
From that one can deduce 

* the (zero temperature) bath correlation function (BCF)
    α(t) = 1/π ∫ dω J(ω) exp(-i ω t)

* the thermal contribution to the BCF
    α_β(t) = 1/π ∫ dω n(β, ω) J(ω) exp(-i ω t)
    with the Bose-Einstein distribution n(β, ω) = 1/(1 - exp(β ω))

* the thermal BCF
    α(β, t) = α(t) + α_β(t) =  1/π ∫ dω J(ω) (2 n(β, ω) cos(ω t) - i sin(ω t))

* as well as the half-sided Fourier-transform of the BCF
    F(ω) = ∫_0^∞ d t α(β, t) exp(i ω t) = J(ω) + i S(ω)
"""
abstract type AbstractBath{T<:AbstractFloat} end



#################################################
##              Ohmic class
#################################################

"""
OhmicExpCO defines the class of (sub-/super-) Ohmic baths with exponential cutoff of the spectral density.

    J(ω) = η ω^s exp(-ω / wc)
"""
struct OhmicExpCO{T} <: AbstractBath{T}
    s::T
    η::T
    wc::T
    c1::T    # used to cache prefactor of the BCF
    "if 'c1' (prefctor of the BCF) is given, η will be set accordingly"
    function OhmicExpCO(s, η=1, wc=1; c1=nothing)
        if c1 === nothing
            c1 = η/π*Γ(s+1)
        else
            η = c1*π/Γ(s+1)
        end
        s, η, wc, c1 = promote(s, η, wc, c1)
        new{typeof(s)}(s, η, wc, c1)
    end     
end



sd(w, b::OhmicExpCO) = b.η * w^b.s * exp(-w/b.wc)
bcf(t, b::OhmicExpCO) = b.c1 * (b.wc/(1.0 + im*b.wc*t) )^(b.s+1)



#################################################
##              multi exponential class
#################################################


"convert G,W -> u, based on G e^-Wt = exp(-u12 - u34t)"
function GW_2_u(G, W)
    T = promote_type(typeof(real(G[1])), typeof(real(W[1])), typeof(log(2)))
    u = Vector{T}(undef, 4*length(G))
    u[1:4:end] .= -real.(log.(G))
    u[2:4:end] .= -imag.(log.(G))
    u[3:4:end] .= real.(W)
    u[4:4:end] .= imag.(W)
    return u
end


"convert u -> G,W, based on G e^-Wt = exp(-u12 - u34t)"
function u_2_GW(u::Vector{T}) where {T <: Real}
    G = exp.( -complex.(u[1:4:end], u[2:4:end]))
    W = complex.(u[3:4:end], u[4:4:end])
    return G, W
end

"""
Defines the class of BCF that take a multi-exponential form, i.e.,
    α(t) = ∑ exp( -a_i - b_i |t|) ≡ ∑ G_i exp(-W_i |t|)
with arbitrary complex parameters G_i = exp(-a_i) and W_i = b_i.
For the sub-class of real positive parameters G_i, the corresponding spectral density
amounts to a sum of Lorentzians. 
Otherwise the spectral in not guaranteed to be real and positive.
Nonetheless, as mathematical vehicle, such a BCF is perfectly sound.
"""
struct MultiExpBCF{T} <: AbstractBath{T}
    u::Vector{T}
    n::Integer
    MultiExpBCF(u::Vector{T}) where T<:Real = length(u) % 4 == 0 ? new{T}(u, length(u) ÷ 4) : throw("length of u must be multiple of 4") 
end
MultiExpBCF(G, W) = MultiExpBCF(GW_2_u(G, W))


"direct acces to the 'u-form' used for fitting"
_bcf_sum_exp(t, u) = sum(@. exp( -complex(u[1:4:end], u[2:4:end]) - complex(u[3:4:end], u[4:4:end])*t ))

function sd(w, b::MultiExpBCF)
    G, W = u_2_GW(b.u)
    return real.(sum.(G / (W - im*w)))
end
bcf(t, b::MultiExpBCF) = _bcf_sum_exp(t, b.u)
