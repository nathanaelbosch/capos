########################################################################################
# Post-Processing: Smoothing and uncertainty calibration
########################################################################################
function calibrate!(integ::ODEFilterIntegrator)

    @unpack state_estimates, sigmas = integ
    for s in state_estimates
        s.Σ .*= sigmas[end]
    end
end




function smooth_all!(integ::ODEFilterIntegrator)

    @unpack state_estimates, times, sigmas = integ
    @unpack A!, Q!, Precond, InvPrecond = integ.constants
    @unpack Ah, Qh = integ.cache
    # x_pred is just used as a cache here

    for i in length(state_estimates)-1:-1:2
        dt = times[i+1] - times[i]

        P = Precond(dt)
        PI = InvPrecond(dt)

        A!(Ah, dt)
        Q!(Qh, dt)
        Qh .*= sigmas[i]

        state_estimates[i] = P * state_estimates[i]
        smooth!(state_estimates[i], P*state_estimates[i+1], Ah, Qh, integ, PI)
        any(isnan.(state_estimates[i].μ)) && error("NaN mean after smoothing")
        any(isnan.(state_estimates[i].Σ)) && error("NaN cov after smoothing")
        state_estimates[i] = PI * state_estimates[i]
    end
end


function smooth!(x_curr, x_next, Ah, Qh, integ, PI=I)

    @unpack d, q = integ.constants
    @unpack x_tmp = integ.cache

    if all((Qh) .< eps(eltype(Qh)))
        @warn "smooth: Qh is really small! The system is basically deterministic, so we just \"predict backwards\"."
        return inv(Ah) * x_next
    end


    # Prediction: t -> t+1
    mul!(x_tmp.μ, Ah, x_curr.μ)
    x_tmp.Σ .= Ah * x_curr.Σ * Ah' .+ Qh


    # Smoothing
    cov_before = copy(x_curr.Σ)
    cov_pred = copy(x_tmp.Σ)
    P_p = Symmetric(cov_pred)
    P_p_inv = inv(P_p)
    G = x_curr.Σ * Ah' * P_p_inv
    x_curr.μ .+= G * (x_next.μ .- x_tmp.μ)

    # Vanilla:
    cov_diff = x_next.Σ .- x_tmp.Σ
    GDG = G * cov_diff * G'
    x_tmp.Σ .= x_curr.Σ .+ GDG

    assert_nonnegative_diagonal(x_curr.Σ)

    return nothing
end
