# BCFUtils

A collection of utilities to manage bath correlation function and related quantities, including fitting of multi-exponential representations in time domain.

# Usage



# Examples

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

$$α(t) = ∑_i G_i e^{-W_i t} ≡ ∑_i e^{-(u^{(i)}_0 + u^{(i)}_1 t)}$$

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

