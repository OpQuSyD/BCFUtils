export fit_ohmic_exp

using SciMLBase
using Optimization
using OptimizationOptimJL
using Sobol
using Distances: minkowski
using NLsolve
using Logging
using Printf
using Serialization

struct OhmicFitCfg{T<:AbstractFloat}
    t_max::T
    s::T 
    num_exp_terms::Int64
    p::T
    u_init_min::NTuple{4, T}
    u_init_max::NTuple{4, T}
    diff_kind::Symbol
    num_tau::Int64
    a_tilde::T
    u_tilde::T
    maxiters::Int64

    function OhmicFitCfg{T}(
        t_max::T, s::T, num_exp_terms::Int64, p::T, u_init_min::NTuple{4, T}, u_init_max::NTuple{4, T}, 
        diff_kind::Symbol, num_tau::Int64, a_tilde::T, u_tilde::T, maxiters::Int64
        ) where T<:AbstractFloat
        0 < t_max|| throw("0 < t_max required")
        0 < s || throw("0 < s required")
        0 < num_exp_terms || throw("0 < num_exp_terms required")
        0 < p || throw("0 < p required")
        0 < p || throw("0 < p required")
        1 < num_tau || throw("1 < num_tau required")
        0 < a_tilde || throw("0 < a_tilde required")
        0 < u_tilde < 1 || throw("0 < u_tilde < 1 required")
        0 < maxiters || throw("0 < maxiters required")
    
        new(t_max, s, num_exp_terms, p, u_init_min, u_init_max, diff_kind, num_tau, a_tilde, u_tilde, maxiters)
    end
end

"convenient init with some default values (convert all floats to T)"
function OhmicFitCfg{T}(
    t_max, s, num_exp_terms, p, u_init_min, u_init_max, 
    diff_kind; num_tau=150, a_tilde=0.1, u_tilde=0.5, maxiters=1000
) where T <: AbstractFloat

    u_init_min = NTuple{4, T}(u_init_min)
    
    OhmicFitCfg{T}(
        convert(T, t_max), 
        convert(T, s), 
        convert(Int64, num_exp_terms),
        convert(T, p), 
        NTuple{4, T}(u_init_min),
        NTuple{4, T}(u_init_max), 
        convert(Symbol, diff_kind), 
        convert(Int64, num_tau),
        convert(T, a_tilde), 
        convert(T, u_tilde), 
        convert(Int64, maxiters)
    )
end

"convenient init with some default values (auto determine float type)"
function OhmicFitCfg(
    t_max, s, num_exp_terms, p, u_init_min, u_init_max, 
    diff_kind; num_tau=150, a_tilde=0.1, u_tilde=0.5, maxiters=1000
)
    T = promote_type(
        typeof(t_max), typeof(s), typeof(p), 
        typeof(u_init_min[1]), typeof(u_init_min[2]), typeof(u_init_min[3]), typeof(u_init_min[4]), 
        typeof(u_init_max[1]), typeof(u_init_max[2]), typeof(u_init_max[3]), typeof(u_init_max[4]), 
        typeof(a_tilde), typeof(u_tilde))

    OhmicFitCfg{T}(
        t_max, s, num_exp_terms, p, u_init_min, u_init_max, diff_kind; 
        num_tau=num_tau, a_tilde=a_tilde, u_tilde=u_tilde, maxiters=maxiters
    )
end



struct OhmicFitSol{T<:AbstractFloat}
    idx::Int64
    u::Vector{T}
    objective::T
    returncode::SciMLBase.ReturnCode.T
end

mutable struct OhmicFitState{T<:AbstractFloat}
    solutions::Vector{OhmicFitSol{T}}
    sobol_idx::Int64
    best_objective::T
end

function OhmicFitState{T}() where {T<:AbstractFloat}
    OhmicFitState(Vector{OhmicFitSol{T}}(undef, 0), 0, convert(T, Inf))
end

function OhmicFitState(::OhmicFitCfg{T}) where {T<:AbstractFloat}
    OhmicFitState(Vector{OhmicFitSol{T}}(undef, 0), 0, convert(T, Inf))
end

# supposed to be optimized (use is since it is included in NLsolve anyway)
"absolute p-norm difference"
diff(::Type{Val{:abs_p_diff}}, f_t::AbstractVector, params::Tuple{AbstractVector, Real}) = minkowski(f_t, params...)

"relative p-norm difference"
function diff(::Type{Val{:rel_p_diff}}, f_t::AbstractVector, params::Tuple{AbstractVector, AbstractVector, Real})
    f_ref, abs2_f_ref, p = params
    sum([(abs2(f_t[i] - f_ref[i]) / abs2_f_ref[i])^(p/2) for i in 1:length(f_ref)])^(1/p)
end

# dispatch diff_kind::Symbol to specific function
diff(diff_kind::Symbol, f_t::AbstractVector, params::Tuple) = diff(Val{diff_kind}, f_t, params)


"calculate the relative p-norm difference"
function diff_from_exp(u, param)
    t_range, params_for_diff, diff_kind = param
    f_t = [exp_bcf_loop(ti, u) for ti in t_range]
    return diff(diff_kind, f_t, params_for_diff)
end


"""for benchmarking only"""
function _diff_ohmic_exp_full_loop(u, param)
    t_range, params_for_diff, diff_kind = param

    @assert diff_kind == :abs_p_diff
    f1_t, p = params_for_diff

    n = Integer(length(u) / 4)

    diff = zero(u[1])

    for (i, ti) in enumerate(t_range)

        a = Complex(u[1], u[2])
        b = Complex(u[3], u[4])
        exp_bcf_ti = exp(-(a + b*ti))
        for i in 1:n-1
            a = Complex(u[4i+1], u[4i+2])
            b = Complex(u[4i+3], u[4i+4])
            exp_bcf_ti += exp(-(a + b*ti))
        end

        diff += abs2(f1_t[i] - exp_bcf_ti)^(p/2)
    end

    return diff^(1/p)
end

function fit_ohmic_exp(tau_test, bcf_tau, p, u0, method, maxiters, diff_kind::Symbol)
    if diff_kind == :abs_p_diff
        params_for_diff = bcf_tau, p
    elseif diff_kind == :rel_p_diff
        params_for_diff = bcf_tau, abs2.(bcf_tau), p
    else
        throw(ArgumentError("unknown diff_kind $diff_kind"))
    end
    f_opt = OptimizationFunction(diff_from_exp, AutoForwardDiff())
    param = (tau_test, params_for_diff, diff_kind)
    prob = OptimizationProblem(f_opt, u0, param)
    sol = solve(prob, method, maxiters=maxiters)
    return sol
end

function get_TauExpRange(ofc::OhmicFitCfg)
    BCF.get_TauExpRange(ofc.t_max, ofc.s, ofc.num_tau, ofc.a_tilde, ofc.u_tilde)
end

"""
# fit exp bcf to Ohmic-kind bcf for many Sobol generated initial conditions

Based on the state given by `ofs::OhmicFitState` continue to process new Sobol samples 
until `num_samples` have been processed.
Only results with lower objective value are added to `ofs::OhmicFitState`.

An InterruptException will stop the routine early and return gracefully with the current state.

# Arguments
- `ofc::OhmicFitCfg{T}`: fit configuration (see `struct OhmicFitCfg` for details)
- `ofs::OhmicFitState{T}`: this mutable struct keeps track of the processed samples 
    and stores sequentially improving results
- `num_samples` the maximum number of samples to process
"""
function sobol_scan_fit!(ofc::OhmicFitCfg{T}, ofs::OhmicFitState{T}, num_samples; verbose=true) where {T<:AbstractFloat}
    # we need the parametric type {T} here, as its needed to instanciate the new solution object

    tau_range = BCF.get_TauExpRange(ofc)
    tau_test = collect(tau_range)
    bcf_tau_test = [BCF.red_ohmic_bcf(ti, ofc.s) for ti in tau_test]
    
    method_NM = NelderMead()
    method_BFGS = BFGS()
       
    # init  low-discrepancy sobol sequence
    sobol_seq = Sobol.SobolSeq(4*ofc.num_exp_terms)
    # skip already checked samples
    Base.skip(sobol_seq, ofs.sobol_idx, exact=true)

    i = ofs.sobol_idx+1

    cnt_new = 0
    was_interrupted = false

    try
        while i <= num_samples
            verbose && print("sample $i, ")
        
            # draw initial condition using sobol samples
            u0 = next!(sobol_seq)
            for i in 1:4
                u0[i:4:end] = ofc.u_init_min[i] .+ (ofc.u_init_max[i]-ofc.u_init_min[i])*u0[i:4:end]
            end 
        
            # run NelderMead followed by BFGS Optimization
            sol_pre = fit_ohmic_exp(tau_test, bcf_tau_test, ofc.p, u0, method_NM, Integer(ceil(ofc.maxiters/10)), ofc.diff_kind)
            verbose && @printf("NelderMead pre solution: fmin %.4e %s, ",sol_pre.objective, sol_pre.retcode)
            sol = fit_ohmic_exp(tau_test, bcf_tau_test, ofc.p, sol_pre.u, method_BFGS, ofc.maxiters, ofc.diff_kind)
            verbose && @printf("BFGS solution: fmin %.4e %s\n", sol.objective, sol.retcode)

            if sol.objective < ofs.best_objective
                cnt_new += 1
                # store new solution which improves previous one
                new_sol = OhmicFitSol{T}(i, sol.u, sol.objective, sol.retcode)
                push!(ofs.solutions, new_sol)
                ofs.best_objective = sol.objective

                verbose && printstyled(@sprintf("new improoved result: idx %i fmin %e %s\n", i, sol.objective, sol.retcode), color=:red)
            end
            i += 1
        end
    catch e
        println()
        @info "caught exception $e"
        isa(e, InterruptException) || throw(e)
        was_interrupted = true
    end
    ofs.sobol_idx = i-1

    return cnt_new, was_interrupted
end


"""
hash the string repr of `ofc::OhmicFitCfg` and return its hex string

### ToDO 

come up with something more robust that `repr(ofc)`
"""
function ofc_identifyer(ofc::OhmicFitCfg)
    r = repr(ofc)
    string(hash(r), base=16)
end

"save fit state, contruct unique file name based on `ofc::OhmicFitCfg`"
function save(ofc::OhmicFitCfg, ofs::OhmicFitState, path=".fit")
    save(ofc_identifyer(ofc), ofc, ofs, path)
end

"save fit state"
function save(fname, ofc::OhmicFitCfg, ofs::OhmicFitState, path=".fit")
    ispath(path) || mkdir(path)
    serialize(joinpath(path, fname), (ofc, ofs))
    @info "fit state saved"
end


"load fit data from an existing file"
function load(full_name)
    deserialize(full_name)
end

"load fit for a given data OhmicFitCfg, return (sobol_state, data)"
function load(ofc::OhmicFitCfg, path=".fit")
    fname = ofc_identifyer(ofc)
    full_name = joinpath(path, fname)

    if ispath(full_name)
        # load from file
        ofc_ff, ofs = load(full_name)
        # sanity check
        @assert ofc_ff == ofc
    else
        # create new state based on config
        ofs = OhmicFitState(ofc)
    end

    return ofs, ofc
end


"load state, run fits, save state and return numbes of processed sobol samples"
function fit(ofc::OhmicFitCfg, num_samples, path=".fit"; verbose=true)
    ofs, _ = load(ofc, path)
    i_0 = ofs.sobol_idx
    _, was_interrupted = sobol_scan_fit!(ofc, ofs, num_samples, verbose=verbose)
    save(ofc, ofs, path)
    was_interrupted && throw(InterruptException())
    return ofs.sobol_idx - i_0
end


function get_fit(full_name::AbstractString, idx::Integer)
    _, ofs = load(full_name)
    idx <= 0 && return ofs.solutions[end]
    return ofs.solutions[idx]
end


function get_fit(ofc::OhmicFitCfg, idx=0, path=".fit")
    full_name = joinpath(path, ofc_identifyer(ofc))
    get_fit(full_name, idx)
end


function u_limits(u)
    u1 = u[1:4:end]
    u2 = u[2:4:end]
    u3 = u[3:4:end]
    u4 = u[4:4:end]

    u_mins = (minimum(u1), minimum(u2), minimum(u3), minimum(u4))
    u_maxs = (maximum(u1), maximum(u2), maximum(u3), maximum(u4))
    return u_mins, u_maxs
end



"""
calculate the parameters u of a single exponential term such that the 
zeroth and the first derivative agrees with an (sub-) Ohmic
bath correlation function at a specific time t

     I)      α(t) = exp(-(c0 + c1 t))
    II)  ddt α(t) = -c1 exp(-(c0 + c1 t)) = -c1 α(t)

    ⇒c1 = - (ddt α(t)) / α(t)
    ⇒c0 = -(log(α(t)) + c1 t)
"""
function smooth_single_exp_to_ohmic_bcf_at_t(t, s)
    al_tau = red_ohmic_bcf(t, s)
    grad_al_tau = BCF.ddt_red_ohmic_bcf(t, s)
    c1 = -grad_al_tau/al_tau
    c0 = -(log(al_tau) + c1*t)
    return [real(c0), imag(c0), real(c1), imag(c1)]
end


"""
calculate the parameters u of a single exponential term such that the 
zeroth and the first derivative agrees with an (sub-) Ohmic
bath correlation function minus an exponential contribution 
(parameterized by u_given) at a specific time t

     I)         α(t) - exp(u_given...) = exp(-(c0 + c1 t))
    II)  ddt [α(t)  - exp(u_given...)] = -c1 exp(-(c0 + c1 t)) = -c1 α(t)

    ⇒c1 = - (ddt α(t)) / α(t)
    ⇒c0 = -(log(α(t)) + c1 t)
"""
function smooth_single_exp_to_ohmic_bcf_minus_exp_at_t(t, s, u_given)
    al_tau = red_ohmic_bcf(t, s) - exp_bcf(t, u_given)
    grad_al_tau = BCF.ddt_red_ohmic_bcf(t, s) - ddt_exp_bcf(t, u_given)
    c1 = -grad_al_tau/al_tau
    c0 = -(log(al_tau) + c1*t)
    return [real(c0), imag(c0), real(c1), imag(c1)]
end


"""
the residual of the nonlinear multi exp system
"""
function _residual_solve_multi_exp_for_ohmic_bcf!(F, u, ts, s)
    n = length(u)
    len_ts = Integer(n/4)
    for i in 1:len_ts
        r = red_ohmic_bcf(ts[i], s) - exp_bcf(ts[i], u)
        F[4*(i-1) + 1] = real(r)
        F[4*(i-1) + 2] = imag(r)
        r_ddt = ddt_red_ohmic_bcf(ts[i], s) - ddt_exp_bcf(ts[i], u)
        F[4*(i-1) + 3] = real(r_ddt)
        F[4*(i-1) + 4] = imag(r_ddt)
    end
end


function solve_multi_exp_for_ohmic_bcf(ts, u0, s)
    residual! = (F, x) -> _residual_solve_multi_exp_for_ohmic_bcf!(F, x, ts, s)
    r = nlsolve(residual!, u0, autodiff = :forward)
    return r
end

function get_n_time_exp_repr_for_ahmic_bcf(ts, s)

    println("add t1 = $(ts[1])")
    u_exact = BCF.smooth_single_exp_to_ohmic_bcf_at_t(ts[1], s)

    for i in 2:length(ts)
        print("add t$i = $(ts[i]) ")        
        # most trivial guess, completely independent of the parameters u that have already been determined
        # ui = BCF.smooth_single_exp_to_ohmic_bcf_at_t(ts[i], s)
    
        # here we get the exp tha is smooth in α(ti, s) - ∑ exp(ti ...)
        ui = BCF.smooth_single_exp_to_ohmic_bcf_minus_exp_at_t(ts[i], s, u_exact)
    
        u_guess = [u_exact ; ui]
        r = BCF.solve_multi_exp_for_ohmic_bcf(ts[1:i], u_guess, s)

        if ! converged(r)
            println("FAILED 🗲 (return latest good solution which has $(Integer(length(u_exact)/4)) exp terms")
            return u_exact
        end

        println("✓")
        u_exact = r.zero
    end
    return u_exact
end