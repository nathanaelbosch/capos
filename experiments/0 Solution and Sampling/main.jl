using ProbNumODE
using ProbNumODE: remake_prob_with_jac
using DifferentialEquations
using LinearAlgebra
using Distributions
using StatsPlots
pyplot()

FILEDIR = dirname(@__FILE__)

# Problem and reference solution
prob = remake_prob_with_jac(lotka_volterra())
appxsol = solve(remake(prob, u0=big.(prob.u0)), abstol=1e-30, reltol=1e-30)

# Probabilistic solve and errors
sol = solve(prob, EKF1(), sigmarule=:fixedMLE, q=5,
            abstol=1e-7, reltol=1e-4, smooth=true)
errors = sol.u .- appxsol.(sol.t)


# Sampling and plotting
sp = ProbNumODE.sample(sol, 10)
sp_errs = sp .- ProbNumODE.stack(sol.u)
p = plot(
    sol.t, ProbNumODE.stack(errors), xlabel="\$t\$", legend=:bottomleft,
    label=["\$(\\hat{y}(t)-y(t))_1\$" "\$(\\hat{y}(t)-y(t))_2\$"],
)
for i in 1:size(sp)[3]
    plot!(p, sol.t, sp_errs[:, :, i], color=[1 2], label="", linewidth=0.5)
end
plot!(sol.t, ProbNumODE.stack(errors), color=:black, label="")
plot!(p, size=(400, 250))
savefig(joinpath(FILEDIR, "sampled_errors_ekf1.png"))

println("Plot saved to $FILEDIR")
