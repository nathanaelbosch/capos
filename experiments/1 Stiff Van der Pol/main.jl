using ProbNumODE
using DifferentialEquations
using Plots
using ProgressMeter
using LinearAlgebra

FILEDIR = dirname(@__FILE__)

prob = ProbNumODE.remake_prob_with_jac(van_der_pol(p=[1e6]))



# The EKF1-IWP3 is able to solve this very stiff problem
println("Timed solve of the VdP problem with EKS1-IBM3:")
@time sol = solve(prob, EKF1(), q=3, smooth=true); sol.destats


# Visualizations
# 1. Solution
p1 = plot(sol.t,
          max.(min.(ProbNumODE.stack(sol.u), 7.0), -7.0),
          xlabel="t", ylabel="y(t)", label="",
          );
savefig(p1, joinpath(FILEDIR, "vdp_solution.png"))


# 2. Step sizes
p2 = plot(sol.t[1:end-1], sol.t[2:end]-sol.t[1:end-1],
          yscale=:log10, xlabel="t", ylabel="dt", label="");
savefig(p2, joinpath(FILEDIR, "vdp_stepsizes.png"))


# 3. Errors
println("Timed very-low-tolerance-solve of the VdP problem with RadauIIA5:")
@time appxsol = solve(
    prob, RadauIIA5(), abstol=1e-14, reltol=1e-14, maxiters=1e7); appxsol.destats
errs = sol.u - appxsol.(sol.t)
ProbNumODE.stack(errs) ./ ProbNumODE.stack(appxsol.(sol.t))

p3 = plot(sol.t, replace(abs.(ProbNumODE.stack(errs)), 0=>missing),
          color=[1 2], label="", yscale=:log10, xlabel="t", ylabel="Error");
plot!(p3, sol.t, replace(sqrt.(ProbNumODE.stack(diag.(sol.pu.Σ))), 0=>missing),
      color=[1 2], linestyle=:dash, label="");
savefig(p3, joinpath(FILEDIR, "vdp_errors.png"))


# Additional plot: Phase plot
u = ProbNumODE.stack(sol.u)
p4 = plot(u[:, 1], u[:, 2], xlabel="y_1(t)", ylabel="y_2(t)", label="");
savefig(p4, joinpath(FILEDIR, "vdp_solution_phase.png"))


println("Plots saved to $FILEDIR")
