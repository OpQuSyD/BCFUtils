# BCFUtils

A collection of utilities to manage bath correlation function and related quantities, including fitting of multi-exponential representations in time domain.

### basics

The package provides the functions:

* **spectral density**, `sd(w, b::AbstractBath)`
* (zero temperature) **bath correlation function**, `bcf(t, b::AbstractBath)`

for some bath `b`.

> **ToDo:** Add at some point
>   * effective thermal spectral density
>   * thermal bath correlation function 
>   * thermal contribution of the BCF

Available baths:

* (sub- / super-) **Ohmic with exponential cutoff**, `OhmicExpCO(s, η, wc)`
* **multi exponential BCF**, `MultiExpBCF(G, W)` or `MultiExpBCF(u)` (see Examples)

> **ToDo:** Use numeric integration to handle more complicated baths

### multi exponential approximation (fitting in time)

A multi exponential representation $α(τ) = ∑_i exp(-u^{(i)}_{12} - u^{(i)}_{34} τ)$ of the BCF has turned out very useful for numerical methods solving open quantum system dynamics.
The subscript "12" and "34" indicate that $u$ is complex valued and, thus, has two independent parameters ($u_{12} = u_1 + i u_2$).
Minimizing the difference to a give given reference BCF amount to a high  dimensional optimization problem.
The search for the global minimum is done by deterministically sampling initial condition for the minimization using Sobol sequences.
The `BFCUtil` package allows to conveniently save the sequences of improved minima as well as the latest Sobol state.
In that way restarting the minimization with the same parameters will continue a previous search.

Parameters for the minimization are stored in 
`struct FitCfg` which has members

* `b::AbstractBath`, the bath `b` defines the reference BCF
* `tau_range::AbstractTauIter`, an iterator which yields the times `t_i` at which the difference with the reference is calculated. The minimizer will then minimize the p-norm of the difference vector `d_i = |bcf_ref(t_i) - bcf_apprx(t_i)|`
* `num_exp_terms`, the number of exponential terms used for the multi-exp-representation
* `p` of the p-norm 
* `u_init_min` minimum of $u$ when sampling initial condition as a four-tuple `[u_1,u_2,_3,u_4]` (same limits for all $u^{(i)}$). Note that $u_1$ is the scale and $u_2$ the phase of the pre-factor. $u_3$ is usually positive, otherwise the there is exponential growth in time. The angular velocity $u_4$ is arbitrary, for Ohmic baths its usually positive.
* `u_init_max` maximum of $u$ when sampling initial condition as a four-tuple `[u_1,u_2,u_3,u_4]`. Since $u_{12}$ comes with a minus sign, large u_1 mean pre-factors with small magnitude.
* `diff_kind`, so far absolute `:abs_p_diff` and relative `:rel_p_diff` has been implemented
* `maxiters`, maximum number of allowed iterations for the minimizer

To trigger the minimization call `fit(ofc::FitCfg, num_samples, path=".fit"; verbose=true)` with arguments
* `ofc`, an instance of the above explained `FitCfg`
* `num_samples`, the maximum number of Sobol samples (initial conditions) to consider. 
* `path=".fit"`, the path where results are saved
* `verbose=false`, if `false` show only a message when a new minima outreaches previous ones, otherwise give detailed information for each sample.

### convenience `FitCfg` for Ohmic baths

> **ToDo:** Explain `FitCfg_for_Ohmic_bath` 

# examples

### Ohmic spectral density with exponential cutoff 

$$J(ω) = η ω^s e^{-ω / ω_c}$$

Instantiate a bath `b` with s=0.3, η=0.1 and wc = 7.
```jldoctest
julia> using BCFUtils
julia> b = OhmicExpCO(0.3, 0.1, 7)
OhmicExpCO{Float64}(0.3, 0.1, 7.0, 0.028567379519453845)
```

Evaluate the spectral density and bath correlation function for the bath `b`.
```jldoctest
julia> sd(2.4, b)
0.09229164933483118
julia> bcf(3.6, b)
-0.0021995074140912564 - 0.004929562274049971im
```

### multi exponential BCF

$$α(τ) = ∑_i G_i e^{-W_i τ} ≡ ∑_i e^{-(u^{(i)}_0 + u^{(i)}_1 τ)}$$

The parameters $G_i$, $W_i$ and, likewise, $u^{(i)}_0$, and $u^{(i)}_1$ can be complex valued.
As an example, a 2-exp-term BCF with`G=[1,2]` and `W=[3,4]` can be created as

```jldoctest
b = MultiExpBCF([1,2im], [3im,4])
MultiExpBCF{Float64}([-0.0, -0.0, 0.0, 3.0, -0.6931471805599453, -1.5707963267948966, 4.0, 0.0], 2)
```

Note that the parameters are stored in $u$-representation with separated real and imaginary part.
So for each exponential term four parameters are stored.

Alternatively the $u$-representation can be passed directly.



# MIT License

Copyright © 2024 Richard Hartmann

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

