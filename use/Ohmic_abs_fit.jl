using BCFUtils

s = 0.5
t_max = 2500
len_tau = 250
N = 5
p = 5
u_init_min = ( 0,  0,  0, 0)
u_init_max = (10, 2π, 10, 5)
diff_kind = :abs_p_diff

b = OhmicExpCO(s, c1=1)
tau_iter = get_TauExpIter_for_Ohmic(t_max, s, len_tau)


maxiters = 500
fc = FitCfg(b, tau_iter, N, p, u_init_min, u_init_max, diff_kind, maxiters)

fit(fc, 15, verbose=true)
sol = get_fit(fc)
u_min, u_max = u_limits(sol.u)
println("u_limits best fit")
println("  u_min $u_min")
println("  u_max $u_max")

