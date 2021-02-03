using ProbNumODE
using ProbNumODE: remake_prob_with_jac, multitime
using DifferentialEquations
using Statistics
using LinearAlgebra
using StatsPlots
using LaTeXStrings

FILEDIR = dirname(@__FILE__)

for (prob, name) in [
    (remake_prob_with_jac(lotka_volterra()), "lotkavolterra"),
    (remake_prob_with_jac(fitzhugh_nagumo_iip()), "fitzhughnagumo"),
    (remake_prob_with_jac(logistic_equation()), "logistic"),
    (remake_prob_with_jac(brusselator()), "brusselator"),
]
    appxsol = solve(remake(prob, u0=big.(prob.u0)), abstol=1e-30, reltol=1e-30)

    densetimes = collect(range(prob.tspan[1], stop=prob.tspan[end], length=1000))
    dense_ref = appxsol.(densetimes)


    function get_res(solve_call)
        sol = solve_call()
        runtime = multitime(solve_call; numruns=10)
        errs = sol.(densetimes) - dense_ref
        L2 = norm(sqrt.(mean.([float.(d) .^ 2 for d in errs])))
        nevals = sol.destats.nf + sol.destats.njacs
        return L2, nevals, runtime
    end


    ekf0_results = []
    ekf0_smooth_results = []
    ekf1_results = []
    ekf1_smooth_results = []
    dp5_results = []
    for (i, (_a, _r)) in enumerate(zip(7:17, 2:14))
        @info "Iteration" name _a _r
        solve_call = () -> solve(prob, EKF0(), sigmarule=:dynamicMV, q=5, abstol=10.0^-_a, reltol=10.0^-_r, smooth=false)
        push!(ekf0_results, get_res(solve_call))

        solve_call = () -> solve(prob, EKF0(), sigmarule=:dynamicMV, q=5, abstol=10.0^-_a, reltol=10.0^-_r, smooth=true)
        push!(ekf0_smooth_results, get_res(solve_call))

        solve_call = () -> solve(prob, EKF1(), sigmarule=:dynamic, q=5, abstol=10.0^-_a, reltol=10.0^-_r, smooth=false)
        push!(ekf1_results, get_res(solve_call))

        solve_call = () -> solve(prob, EKF1(), sigmarule=:dynamic, q=5, abstol=10.0^-_a, reltol=10.0^-_r, smooth=true)
        push!(ekf1_smooth_results, get_res(solve_call))

        solve_call = () -> solve(prob, DP5(), abstol=10.0^-_a, reltol=10.0^-_r)
        push!(dp5_results, get_res(solve_call))
    end


    p = plot([r[2] for r in ekf0_results],
             [r[1] for r in ekf0_results],
             xscale=:log10,
             yscale=:log10,
             marker=:auto,
             label="EKF0",
             );
    plot!(p, [r[2] for r in ekf0_smooth_results],
          [r[1] for r in ekf0_smooth_results],
          marker=:auto,
          label="EKS0",
          );
    plot!(p, [r[2] for r in ekf1_results],
          [r[1] for r in ekf1_results],
          marker=:auto,
          label="EKF1",
          );
    plot!(p, [r[2] for r in ekf1_smooth_results],
          [r[1] for r in ekf1_smooth_results],
          marker=:auto,
          label="EKS1",
          );
    plot!(p, [r[2] for r in dp5_results],
          [r[1] for r in dp5_results],
          marker=:auto,
          label="DP5",
          );

    x0_dp5 = dp5_results[1][2]
    y0_dp5 = dp5_results[1][1]
    xmax_dp5 = dp5_results[end][2]
    xrange = 10 .^ range(log10(x0_dp5), log10(xmax_dp5), length=1000)
    plot!(p, xrange, x -> x^(-6)*(x0_dp5^6 * y0_dp5), color=:black, linestyle=:dot, label=L"$h^6$ conv.")

    x0 = ekf1_results[1][2]
    y0 = ekf1_results[1][1]
    xrange = 10 .^ range(log10(x0), log10(xmax_dp5), length=1000)
    plot!(p, xrange, x -> x^(-6)*(x0^6 * y0), color=:black, linestyle=:dot, label="")

    x0 = ekf1_smooth_results[1][2]
    y0 = ekf1_smooth_results[1][1]
    xrange = 10 .^ range(log10(x0), log10(xmax_dp5), length=1000)
    plot!(p, xrange, x -> x^(-6)*(x0^6 * y0), color=:black, linestyle=:dot, label="")

    plot!(p, xlabel=L"\#evaluations", ylabel="L2 error")
    savefig(p, joinpath(FILEDIR, "dp5_comparison_$(name).png"))
end
