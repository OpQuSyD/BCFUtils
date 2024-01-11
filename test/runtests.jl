using Revise

using BCFUtils
using Test

x = NTuple{2, Float64}((2, 2))
convert(NTuple{2, Float16}, x)

###########################################################
# simple tests of the Ohmic sd with exponential cutoff

b = OhmicExpCO(0.2)
@test BCFUtils.typeof_bath(b) == typeof(2.0)

b = OhmicExpCO(1)
@test BCFUtils.typeof_bath(b) == typeof(2.0)

b = OhmicExpCO(1, c1=1)
@test BCFUtils.typeof_bath(b) == typeof(2.0)

# explicit type given for bath
b = OhmicExpCO{Float32}(1)
@test BCFUtils.typeof_bath(b) == Float32

setprecision(BigFloat, 400)
b = OhmicExpCO(BigFloat(1))
@test BCFUtils.typeof_bath(b) == BigFloat

# check specific values
eta = 2.3
b = OhmicExpCO(0.5, eta, 10000)
@test abs(sd(1, b) - eta) < 0.001

b = OhmicExpCO(1, c1=1)
@test bcf(0, b) == 1
@test abs(bcf(1, b)) == 1/2

# type consistency for bcf call
b = OhmicExpCO{Float16}(1, c1=1)
t = Float16(1.0)
w = Float16(2.0)
@test typeof(bcf(t, b)) == ComplexF16
@test typeof(sd(w, b)) == Float16

p = 400
setprecision(BigFloat, p)
b = OhmicExpCO{BigFloat}(1, c1=1)
@test precision(b.s) == p
@test precision(b.c1) == p
t = BigFloat("2.3")
c = bcf(t, b)
@test typeof(c) == Complex{BigFloat}
@test precision(real(c)) == p

w = BigFloat("3.3")
c = sd(w, b)
@test typeof(c) == BigFloat
@test precision(real(c)) == p


###########################################################
# simple tests of the exp bcf

# test instanciation -- auto convert int to float
G = [2 + 0im]
W = [0 + 1im]
b = MultiExpBCF(G, W)
@test eltype(b.u) == typeof(log(3))

# test instanciation -- auto promote BigFloat
G = [BigFloat(2, precision=20) + 0im]
W = [0 + 1im]
b = MultiExpBCF(G, W)
@test eltype(b.u) == BigFloat

# auto deduce Float32 as type from G and W
G = ComplexF32[2 + 0im]
W = ComplexF32[0 + 1im]
b = MultiExpBCF(G, W)
@test eltype(b.u) == Float32

# call with explicit Type

G = [2 + 0im]
W = [0 + 1im]
b = MultiExpBCF{Float64}(G, W)
@test eltype(b.u) == Float64

# test specific value
@test bcf(0, b) == 2
@test bcf(1, b) ≈ 2 * exp(-1im)
@test bcf(pi, b) ≈ -2 atol=1e-12


###########################################################
# test TauExpRange

n = 16001
t_max = 1500
al_tilde = 0.1
s = 0.5

# test tmax value
for u_tilde in [0.4, 0.5, 0.8]
    tau0, a = BCFUtils.calc_tau_0_a(t_max, s, al_tilde, u_tilde)
    # final value (u=1)
    local u = 1
    @test tau0*(exp(a*u) - 1) ≈ t_max atol=1e-8
end

# check that α(t(u_tilde)) = α_tilde
for s in [0.3, 0.8, 1.6]
    local b = OhmicExpCO(s, c1=1)
    for u_tilde in [0.2, 0.5, 0.8]
        tr = get_TauExpIter_for_Ohmic(t_max, s, n, al_tilde, u_tilde)
        @test length(tr) == n
        tr = collect(tr)
        @test length(tr) == n
        @test tr[1] == 0                  # initial value (u=0)
        @test tr[end] ≈ t_max atol=1e-8   # final value (u=1)
        tau_u = tr[Integer(floor(u_tilde*n))]
        @test abs(bcf(tau_u, b)) ≈ al_tilde atol=al_tilde/10
    end
end


###########################################################
# minimize difference between exp bcf and ohmic bcf

@test iszero(BCFUtils.diff(:abs_p_diff, [1,2], ([1,2], 2)))
@test iszero(BCFUtils.diff(:abs_p_diff, [1,2], ([1,2], 5.67)))


d1 = 2
d2 = 3.4
# abs diff
@test BCFUtils.diff(:abs_p_diff, [2,5im], ([2+d1,(5+d2)*1im], 3)) ≈ (d1^3 + d2^3)^(1/3) atol=1e-12
# relative diff, actually weighted diff, here with weights 1.0 -> so same result as abs diff
@test BCFUtils.diff(:rel_p_diff, [2,5im], ([2+d1,(5+d2)*1im], [1.0, 1.0], 3)) ≈ (d1^3 + d2^3)^(1/3) atol=1e-12

@test_throws MethodError BCFUtils.diff(:new_diff, [2,5im], ([2+d1,(5+d2)*1im], [1.0, 1.0], 3))


n = 150
t_max = 100
s = 0.5
b = OhmicExpCO(s, c1=1)
tr = get_TauExpIter_for_Ohmic(t_max, s, n)
tr = collect(tr)
bcf_tau = [bcf(ti, b) for ti in tr]
p = 5
u0 = [1,1,1,1, 0.5, 0.5, 0.5, 0.5]
method = BCFUtils.NelderMead()
maxiters = 1000

sol = BCFUtils.fit_bcf_to_exp(tr, bcf_tau, p, u0, method, maxiters, :abs_p_diff)
@test Integer(sol.retcode) == 1
@test sol.objective ≈ 0.05585703941513527 atol=1e-10

sol = BCFUtils.fit_bcf_to_exp(tr, bcf_tau, p, u0, method, maxiters, :rel_p_diff)
@test Integer(sol.retcode) == 1
@test sol.objective ≈ 0.979426791250466 atol=1e-10



t_max = 100
s = 0.5
num_exp_terms = 4
p = 5
u_init_min = (0, 0, 0, 0)
u_init_max = (10, 2pi, 10, 2pi)
diff_kind = :abs_p_diff

ofc = BCFUtils.FitCfg_for_Ohmic_bath(
        t_max, s, num_exp_terms, diff_kind, 
        p=p, 
        u_init_min=u_init_min, 
        u_init_max=u_init_max,
    )
r = repr(ofc)
@test hash(r) == 0x2dabfc9e34870160
ofs = BCFUtils.FitState(ofc)

ofc2 = BCFUtils.FitCfg_for_Ohmic_bath(
    t_max, s, num_exp_terms, diff_kind, 
    p=p, 
    u_init_min=u_init_min, 
    u_init_max=u_init_max
)

@test ofc == ofc2

setprecision(20)
ofcBF = BCFUtils.FitCfg_for_Ohmic_bath(
    BigFloat(t_max), s, num_exp_terms, diff_kind, 
    p=p, 
    u_init_min=u_init_min, 
    u_init_max=u_init_max
)
r = repr(ofcBF)
@test hash(r) == 0xb0e5279e46b7e4e5
ofcBF = BCFUtils.FitState(ofcBF)
@test ofc != ofcBF


ofcBF = BCFUtils.FitCfg_for_Ohmic_bath(
    t_max, BigFloat(s), num_exp_terms, diff_kind, 
    p=p, 
    u_init_min=u_init_min, 
    u_init_max=u_init_max
)
r = repr(ofcBF)
@test hash(r) == 0xb0e5279e46b7e4e5
ofcBF2 = BCFUtils.FitState(ofcBF)
@test ofcBF2 != ofcBF


# crunch the first 5 samples
ofs = BCFUtils.FitState(ofc)
cnt_new, _ = BCFUtils.sobol_scan_fit!(ofc, ofs, 5)
@test ofs.sobol_idx == 5
@test cnt_new == length(ofs.solutions)

# continue until we have 10 samples
cnt_new2, _ = BCFUtils.sobol_scan_fit!(ofc, ofs, 10)
@test cnt_new + cnt_new2 == length(ofs.solutions)
@test ofs.sobol_idx == 10

path = "test_fits"
fname = BCFUtils.fc_identifyer(ofc)
fullname = joinpath(path, fname)
ispath(fullname) && rm(fullname)
c = BCFUtils.fit(ofc, 10, path, verbose=false)
@test c == 10

c = BCFUtils.fit(ofc, 10, path, verbose=false)
@test c == 0


