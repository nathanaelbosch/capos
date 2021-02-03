using ProbNumODE
using ProbNumODE: remake_prob_with_jac, multitime, make_savable, compute_errors
using BSON: @save
using DifferentialEquations
using Plots
using ProgressMeter

using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_lotkavoltera, prob_ode_rigidbody
lotkavolterra = remake(prob_ode_lotkavoltera, tspan=(0.0, 10.0))


failed_runs = []


problist = [
    (logistic_equation(), "logistic"),
    (lotkavolterra, "lotkavolterra"),
    (fitzhugh_nagumo_iip(), "fitzhughnagumo"),
    (prob_ode_rigidbody, "rigidbody"),
    (brusselator(), "brusselator"),
]

for (prob, probname) in problist
    @info "New Problem: $probname"
    prob = remake_prob_with_jac(prob)
    appxsol = DiffEqBase.has_analytic(prob.f) ? nothing :
        solve(remake(prob, u0=big.(prob.u0)), abstol=1e-18, reltol=1e-18)


    result_list = []


    algs = [EKF0(), EKF1()]
    qs = 1:5
    sigmas = [:fixedMLE, :fixedMLEMV, :dynamic, :dynamicMV,]
    steprules = [:standard]
    smooths = [false, true]

    abstols, reltols = Dict(
        "logistic" => (4:13, 1:10),
        "lotkavolterra" => (5:13, 2:10),
        "fitzhughnagumo" => (4:13, 1:10),
        "rigidbody" => (4:13, 1:10),
        "brusselator" => (4:13, 1:10),
    )[probname]
    abstols, reltols = 10.0 .^ -abstols, 10.0 .^ -reltols

    n_combinations = length(algs) * length(qs) * length(steprules) * length(sigmas) * length(smooths) * length(abstols)
    p = Progress(n_combinations, 1, "Working on $probname: ")

    for alg in algs,
        q in qs,
        sigmarule in sigmas,
        steprule in steprules,
        smooth in smooths,
        (abstol, reltol) in zip(abstols, reltols)

        next!(p, showvalues=[(:alg, alg), (:q, q), (:sigmarule, sigmarule),
                             (:smooth, smooth),
                             (:abstol, abstol), (:reltol, reltol)])

        if alg isa EKF1 && sigmarule in (:dynamicMV, :fixedMLEMV) continue end

        config = (
            alg=string(nameof(typeof(alg))),
            q=q,
            steprule=steprule,
            sigmarule=sigmarule,
            abstol=abstol,
            reltol=reltol,
            smooth=smooth,
        )

        @debug "Run started:" probname config

        try

            solve_call = () -> solve(prob, alg; config..., maxiters=1e6)
            sol = solve_call()
            if sol.retcode != :Success
                if sol.retcode == :MaxIters
                    # @info "MaxIters" sol.retcode probname config
                else
                    @warn "Run not successful!" sol.retcode probname config
                    push!(failed_runs, (sol.retcode, probname, config))
                end
                continue
            end

            runtime = multitime(solve_call; numruns=10)
            summaries = Dict(
                :runtime=>runtime,
                :nf => sol.destats.nf,
                :num_evals => sol.destats.nf + sol.destats.njacs,
                :naccept => sol.destats.naccept,
                :nreject => sol.destats.nreject,
            )
            savable_sol = make_savable(sol; appxsol=appxsol)
            final_vals = (
                x=savable_sol.x[end],
                t=savable_sol.t[end],
                pu=savable_sol.pu[end],
                u_ref=savable_sol.u_ref[end],
            )

            errors = compute_errors(savable_sol)
            run = (
                problem=probname,
                config=config,
                # solution=savable_sol,
                final_vals=final_vals,
                summaries=merge(summaries, errors),
                destats=sol.destats,
            )

            push!(result_list, run)
            @debug "Run finished:" run
        catch err
            if err isa InterruptException
                throw(err)
            end
            @warn "Encountered error!" err probname config
            push!(failed_runs, (err, probname, config))
            # throw(err)
        end
    end

    varname = "results"
    @save "data/$(probname)_fixedMLEMV.bson" result_list
end

@info "All unsuccessful runs again:" failed_runs
