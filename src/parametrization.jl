using Random
using GlobalSensitivityAnalysis
using Distributions
using DataStructures

Random.seed!(1);

function build_p_ci(p_ci)
    # Parameters
    lambda = (1/p_ci[3] - 1)/4
    p = [p_ci[1], p_ci[2], 1/4, lambda, 1/7, p_ci[4], p_ci[5]*p_ci[6], 0.0996, p_ci[5], 0.0766, 0.0026, 0.0843, 0.0017];

    # Initial condition
    quar = p_ci[3] * 6718903.0
    a_plus_i = 498.0 * (1/taxa_diagn(p) - 1)
    a = a_plus_i * (1 - p_ci[4])
    i = a_plus_i * p_ci[4]
    e = 498.0 * p_ci[7]
    ci_base = [quar, e, a, i]
    # Calculando a população fora da quarentena
    s = 498.0 + 16.
    for i = 1:4
        s += ci_base[i]
    end
    pop_out = 6718903.0 - s
    full_ci = [pop_out, ci_base[1:4]..., 0.1*498.0, 0.9*498.0, 0.0, 498.0+16.0, 16.0]
    
    full_ci, p
end

fator_assint(p) = 1 - p[6]

function taxa_diagn(p)
    p[9]/(p[9] + p[10] + p[11])*p[6] + p[7]/(p[7] + p[8])*(1 - p[6])
end

function taxa_mort(p)
    p[11]/(p[9] + p[10] + p[11]) + p[13]/(p[12] + p[13]) * p[9]/(p[9] + p[10] + p[11])
end

function taxa_quarentena(p)
    p[3]/(p[3] + p[4])
end

function R0(CI,p)
    p[1]*(p[2]*(1-p[6])/(p[7]+p[8])+p[6]/(p[9]+p[10]+p[11]))*CI[1]
end

function Rt(sol,p)
    p[1]*(p[2]*(1-p[6])/(p[7]+p[8])+p[6]/(p[9]+p[10]+p[11]))*sol[1,:]
end

function dRt(sol,p)
    p[1]*(p[2]*(1-p[6])/(p[7]+p[8])+p[6]/(p[9]+p[10]+p[11]))*((-p[1]*sol[1,:]).*(sol[5,:] .+ p[2] .* sol[4,:]) .- p[3] .* sol[1,:] .+ p[4] .* sol[2,:])
end

function pp_p_ci(p_ci)
    full_ci, p = build_p_ci(p_ci)
    print("""População em quarentena = $(full_ci[2]) = $(100*full_ci[2]/6718903) %
              Exposta       = $(full_ci[3])
              Assintom      = $(full_ci[4])
              Sintomática   = $(full_ci[5])
        
    Taxa de contágio    = $(R0(full_ci,p))
    Fator cont assintom = $(p[2])
    Em quarentena       = $(100*p[3]/(p[3] + p[4])) % 
    Prop sintomática    = $(p[6])
    Taxa testagem A     = $(p[7]) ~ $(1/p[7]) dias
    Taxa testagem I     = $(p[9]) ~ $(1/p[9]) dias
    """)
end

p_start = [1e-7, 0.1, 0.6, 0.88, 1/30, 0.3, 7.0];

function loss(sol_)
    v = sol_[9,:]
    new_diagn = v[2:n_pts+1] - v[1:n_pts]
    v = sol_[10,:]
    new_deaths = v[2:n_pts+1] - v[1:n_pts]
    return sum(abs2,diagn_Rio .- new_diagn) + 16^2 * sum(abs2,obitos_Rio .- new_deaths)
end

function loss_obj(p_ci)
    # Initial condition and Parameters
    full_ci, p = build_p_ci(p_ci)
    
    prob = ODEProblem(g, full_ci, (0.0,float(n_pts+1)), p)
    # build_loss_objective(prob, Tsit5(), loss, maxiters=10000, verbose=false)
    sol = Array(concrete_solve(prob,Feagin14(), full_ci, p, saveat=0.0:1.0:float(n_pts+1)))
    loss(sol)
end

function build_params_variations(res)
    min_arrs = [[res.minimizer[1] * taxa for taxa in collect(0.7:0.05:1.3)], # Beta
            [res.minimizer[2]], # Fator assintomatico
            [res.minimizer[3] + x for x in collect(-0.2:0.02:0.2)], # Quarentena
            [res.minimizer[4]], # Prop sintomatica
            [1/(1/res.minimizer[5] + x) for x in collect(-20:2:20)], 
            [res.minimizer[6]],
            [res.minimizer[7]*taxa for taxa in collect(0.5:0.5:2.0)],
    ];
end

function build_possibilities(min_arrs)
    possibilities = Dict(
    "1" => min_arrs[1],
    "2" => min_arrs[2],
    "3" => min_arrs[3],
    "4" => min_arrs[4],
    "5" => min_arrs[5],
    "6" => min_arrs[6],
    "7" => min_arrs[7]
    )
    hr_poss = Dict(
        "taxa_contagio" => min_arrs[1],
        "fator_assint" => min_arrs[2],
        "quarentena" => [100*(1/4)/(1/4 + ((1/p - 1)/4)) for p in min_arrs[3]],
        "prop_sint" => min_arrs[4],
        "taxa_teste_a" => sort(vcat([1 / (p5 * p6) for (p5,p6) in Iterators.product(min_arrs[5],min_arrs[6])]...)),
        "taxa_teste_i" => 1 ./ min_arrs[5]
    )
    return possibilities, hr_poss
end

function num_series(possibilities)
    N_series = 1
    poss_keys = collect(keys(possibilities))
    for k in poss_keys
        arr = possibilities[k]
        b = length(arr)
        N_series = N_series*b
    end
    N_series
end

function make_iter(possibilities)
    a = ()
    poss_keys = collect(keys(possibilities))
    for k in poss_keys
        arr = possibilities[k]
        b = length(arr)
        a = tuple(a...,b...)
    end
    return CartesianIndices(a)
end

function param_range(possibilities)
    iter = make_iter(possibilities)
    arr_poss = []
    N = length(possibilities)
    poss_keys = collect(keys(possibilities))
    for idxs in iter
        p_ci = []
        for k = 1:N
            push!(p_ci,possibilities[string(k)][idxs[findfirst(x -> x == string(k),poss_keys)]])
        end
        full_ci, p = build_p_ci(p_ci)
        poss_dict = Dict(
            "taxa_contagio" => p[1],
            "fator_assint" => p[2],
            "quarentena" => (100*p[3]/(p[3] + p[4])),
            "prop_sint" => p[6],
            "taxa_teste_a" => 1/p[7],
            "taxa_teste_i" => 1/p[9]
        )
        push!(arr_poss, [full_ci,p,poss_dict])
    end
    return arr_poss
end

function param_range_uq(samples)
    arr_poss = []
    for idxs in eachrow(samples)
        full_ci, p = build_p_ci(idxs)
        poss_dict = Dict(
            "taxa_contagio" => p[1],
            "fator_assint" => p[2],
            "quarentena" => (100*p[3]/(p[3] + p[4])),
            "prop_sint" => p[6],
            "taxa_teste_a" => 1/p[7],
            "taxa_teste_i" => 1/p[9]
        )
        push!(arr_poss, [full_ci,p,poss_dict])
    end
    return arr_poss
end

likelihood(res; σ=1.0) = exp(-loss(res)/(2*σ))

function prob_adjust(model, prior_param_prob, likelihood, param_range; ts = 1.0:1.0:200.0, tol = 0)
    all_results = []
    for poss in param_range
        CI, p, poss_dict = poss[1], poss[2], poss[3]
        res = model(p, CI, ts)
        bayes_rule = prior_param_prob(p, CI) * likelihood(res)
        if bayes_rule < tol
            continue
        end
        push!(all_results, (Dict("res" => res, "p" => p, "CI" => CI, "poss" => poss_dict), bayes_rule))
    end
    all_results
end

function build_random_samples(possibilities; n_samples = 2000)
    N = num_series(possibilities)
    Random.seed!(1);
    Random.rand(1:N, n_samples);
end

function build_series_subnot_infec(series, idxs, ts)
    xs = [[series[j][1]["res"](t)[5]/series[j][1]["res"](t)[7] for j in idxs] for t in ts]
    ws = [[series[j][2] for j in idxs] for t in ts]
    return xs, ws
end

function build_series_subnot_assint(series, idxs, ts)
    xs = [[series[j][1]["res"](t)[4]/series[j][1]["res"](t)[6] for j in idxs] for t in ts]
    ws = [[series[j][2] for j in idxs] for t in ts]
    return xs, ws
end

function build_series(series, idxs, ts; index = 10)
    xs = [[series[j][1]["res"](t)[index] for j in idxs] for t in ts]
    ws = [[series[j][2] for j in idxs] for t in ts]
    return xs, ws
end

function build_series_daily(series, idxs, ts; index = 10)
    xs = [[series[j][1]["res"](t)[index]-series[j][1]["res"](t-1)[index] for j in idxs] for t in ts]
    ws = [[series[j][2] for j in idxs] for t in ts]
    return xs, ws
end

function build_series_parameter(series, parameter, possibilities, hr_poss; day = 40, samples = 100, serie = 10)
    N_series = 1
    poss_keys = collect(keys(possibilities))
    for k in poss_keys
        arr = possibilities[k]
        b = length(arr)
        N_series = N_series*b
    end
    N = length(hr_poss[parameter])
    xs = [Float64[] for i = 1:N]
    ws = [Float64[] for i = 1:N]
    Random.seed!(2)
    ℓ = samples
    if samples == "all"
        idxs = collect(1:1:N_series)
    else
        @assert N_series >= samples
        idxs = Random.rand(1:N_series, ℓ);
    end
    for j in idxs
        if floor(parameter != "taxa_contagio")
            for i = 1:N
                if Int(floor(series[j][1]["poss"][parameter])) == Int(floor(hr_poss[parameter][i]))
                    push!(xs[i], series[j][1]["res"](day)[serie])
                    push!(ws[i], series[j][2])
                end
            end
        else
            for i = 1:N
                if series[j][1]["poss"][parameter] == hr_poss[parameter][i]
                    push!(xs[i], series[j][1]["res"](day)[serie])
                    push!(ws[i], series[j][2])
                end
            end    
        end
    end
    return xs,ws
end

