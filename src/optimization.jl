using DiffEqFlux,Flux, Optim

cb = function (p,l,pred) #callback function to observe training
    display(l)
    return false # Tell it to not halt the optimization. If return true, then optimization stops
end

# TODO: use "new_diagn, new_death" as monitors?
function loss_adjoint(p_ci)
    # Initial condition and Parameters
    full_ci, p = build_p_ci(p_ci)
    
    prob = ODEProblem(g, full_ci, [0.0, float(n_pts)], p)
    sol = Array(concrete_solve(prob,Feagin14(), full_ci, p, saveat=0.0:1.0:float(n_pts)))
    loss(sol), sol
end