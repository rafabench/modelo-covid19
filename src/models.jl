using DifferentialEquations, ParameterizedFunctions

g = @ode_def Jia begin
    dS = - β*S*(I + θ*A) - p*S + λ*Q
    dQ = p*S - λ*Q
    dE = β*S*(I + θ*A) - σ*E
    dA = σ*(1-ρ)*E - ϵ_a*A - γ_a*A
    dI = σ*ρ*E - ϵ_i*I - γ_i*I - d_i*I
    dDA = ϵ_a*A - γ_a*DA
    dDI = ϵ_i*I - γ_d*DI - d_d*DI
    dR = γ_a*A + γ_i*I + γ_a*DA + γ_d*DI
    # Extra dimensions for fitting: cumulative cases and deaths
    dC = ϵ_a*A + ϵ_i*I
    dM = d_i*I + d_d*DI
end β θ p λ σ ρ ϵ_a γ_a ϵ_i γ_i d_i γ_d d_d;
#   1 2 3 4 5 6  7   8   9   10  11  12  13

# Idx       :         1,    2,   3,      4,     5,    6,     7,      8,    9,     10,     11,     12,      13
# Parâmetros:      beta, theta,   p, lambda, sigma,  rho,  epsA, gammaA, epsI, gammaI, deathI, gammaD, deathD
params_rio_base = [1e-7,  1/10, 2/9,   2/11,   1/7, 0.55, 1/200, 0.0996, 1/50, 0.0766, 0.0026, 0.0843, 0.0017];

popRio = 6718903.0
sucetiveis = 0.4*popRio - 15.0*498.0
# Idx  :        1,          2,        3,         4,          5,             6,          7,       8,     9,     10
# CI   :  sucetíveis; quarentena;  expostos ;  assintom; sintomáticos; diagnost A; diagnost I; recup; tot_d; mortos
CI_Rio = [sucetiveis; 0.6*popRio; 7.0*498.0 ; 8.0*498.0;    8.0*498.0;  0.1*498.0;  0.9*498.0;   0.0; 498.0;   15.0];

function model(params, CI, ts)
    u0 = CI
    p = params
    prob = ODEProblem(g, u0, [ts[1], ts[end]], p)
    tmp_prob = remake(prob, u0 = u0, tspan = [ts[1], ts[end]], p = p)
    sol = solve(tmp_prob,Feagin14(),abstol=1e-14,reltol=1e-14,saveat=ts);
end

function piecewise_model(params, CI, ts, cut_at, new_params)
    u0 = CI
    p = params
    ts_new = cut_at[1]+1.0:1.0:ts[end]
    ts = ts[1]:1.0:cut_at[1]
    prob = ODEProblem(g, u0, [ts[1], ts[end]], p)
    sol = solve(prob,Feagin14(),abstol=1e-14,reltol=1e-14,saveat=ts);
    tmp_prob = remake(prob, u0 = sol[end], tspan = [ts_new[1], ts_new[end]], p=new_params)
    sol_1 = solve(tmp_prob,Feagin14(),abstol=1e-14,reltol=1e-14,saveat=ts_new);
    return hcat(hcat(sol.u[1:end-1]...),hcat(sol_1.u...))
end