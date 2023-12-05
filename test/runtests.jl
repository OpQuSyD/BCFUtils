using Revise

using BCFUtils
using Test


###########################################################
# simple tests of the exp bcf

# test instanciation -- convert int to float
G = [2 + 0im]
W = [0 + 1im]
b = MultiExpBCF(G, W)
@test eltype(b.u) == typeof(log(3))

# test instanciation -- promote BigFloat
G = [BigFloat(2, precision=20) + 0im]
W = [0 + 1im]
b = MultiExpBCF(G, W)
@test eltype(b.u) == BigFloat

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
@test hash(r) == 0x95edec8626492c31
ofs = BCFUtils.FitState(ofc)

ofc2 = BCFUtils.FitCfg_for_Ohmic_bath(
    t_max, s, num_exp_terms, diff_kind, 
    p=p, 
    u_init_min=u_init_min, 
    u_init_max=u_init_max
)

@test ofc == ofc2

ofcBF = BCFUtils.FitCfg_for_Ohmic_bath(
    BigFloat(t_max, precision=20), s, num_exp_terms, diff_kind, 
    p=p, 
    u_init_min=u_init_min, 
    u_init_max=u_init_max
)
r = repr(ofcBF)
@test hash(r) == 0x00c525b5d8bc586c
ofcBF = BCFUtils.FitState(ofcBF)

@test ofc != ofcBF



# cunrch the first 5 samples
ofs = BCFUtils.FitState(ofc)
cnt_new, _ = BCFUtils.sobol_scan_fit!(ofc, ofs, 5)
@test ofs.solutions[1].idx == 1
@test ofs.solutions[end].idx == 4
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

println("all tests passed")

