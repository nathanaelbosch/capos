using ProbNumODE
using DifferentialEquations
using Plots
using LinearAlgebra
using LaTeXStrings


function makeplot(sol, n=10)
    errs = sol.u - appxsol.(sol.t)
    p = plot(
        xlabel=L"$t$", ylabel=L"$(\hat{y}(t) - y^*(t))$",
    );
    stds = sqrt.(ProbNumODE.stack(diag.(sol.pu.Σ)))
    plot!(sol.t, [zero(sol.t) zero(sol.t)], ribbon=3stds,
          color=[1 2], linealpha=0, fillalpha=0.2, label="");

    if n > 0
        sp = ProbNumODE.sample(sol, n)
        sp_errs = ProbNumODE.stack(sol.u) .- sp
        for i in 1:size(sp)[3]
            plot!(p, sol.t, sp_errs[:, :, i], color=[1 2], label="", linewidth=0.5)
        end
        plot!(p, sol.t, ProbNumODE.stack(errs), color=:black, label="");
    else
        plot!(p, sol.t, ProbNumODE.stack(errs), color=[1 2], label="");
    end

    return p
end



prob = ProbNumODE.remake_prob_with_jac(fitzhugh_nagumo_iip())
appxsol = solve(remake(prob, u0=big.(prob.u0)), abstol=1e-30, reltol=1e-30);


sol = solve(prob, EKF1(), q=3, sigmarule=:dynamic, smooth=true, abstol=1e-10, reltol=1e-7)
p = makeplot(sol, 0);
savefig("sigmas-tv-eks1_fitzhughnagumo.png")


sol = solve(prob, EKF0(), q=3, sigmarule=:dynamicMV, smooth=true, abstol=1e-10, reltol=1e-7)
p = makeplot(sol, 0);
savefig("sigmas-tvmv-eks0_fitzhughnagumo.png")


sol = solve(prob, EKF0(), q=3, sigmarule=:dynamic, smooth=true, abstol=1e-10, reltol=1e-7)
p = makeplot(sol, 0);
savefig("sigmas-tv-eks0_fitzhughnagumo_nosamples.png")
