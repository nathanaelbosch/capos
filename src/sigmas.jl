abstract type AbstractSigmaRule end
abstract type AbstractStaticSigmaRule <: AbstractSigmaRule end
abstract type AbstractDynamicSigmaRule <: AbstractSigmaRule end
isstatic(sigmarule::AbstractStaticSigmaRule) = true
isdynamic(sigmarule::AbstractStaticSigmaRule) = false
isstatic(sigmarule::AbstractDynamicSigmaRule) = false
isdynamic(sigmarule::AbstractDynamicSigmaRule) = true
initial_sigma(sigmarule::AbstractSigmaRule, d, q) = 1.0

struct MLESigma <: AbstractStaticSigmaRule end
function static_sigma_estimation(rule::MLESigma, integ)
    @unpack d = integ.constants
    @unpack measurement = integ.cache

    v, S = measurement.μ, measurement.Σ

    if iszero(v)
        return zero(integ.cache.σ_sq)
    end
    if iszero(S)
        return Inf
    end

    sigma_t = v' * inv(S) * v / d

    if integ.success_iter == 0
        @assert length(integ.sigmas) == 0
        return sigma_t
    else
        @assert length(integ.sigmas) == integ.success_iter
        sigma_prev = integ.sigmas[end]
        sigma = sigma_prev + (sigma_t - sigma_prev) / integ.success_iter
        return sigma
    end
end


struct DynamicSigma <: AbstractDynamicSigmaRule end
function dynamic_sigma_estimation(kind::DynamicSigma, integ)
    @unpack d, R = integ.constants
    @unpack h, H, Qh = integ.cache
    σ² = h' * inv(H*Qh*H') * h / d
    return σ²
end


struct MVDynamicSigma <: AbstractDynamicSigmaRule end
initial_sigma(sigmarule::MVDynamicSigma, d, q) = kron(ones(q+1, q+1), diagm(0 => ones(d)))
function dynamic_sigma_estimation(kind::MVDynamicSigma, integ)
    @unpack dt = integ
    @unpack d, q, R, InvPrecond, E1 = integ.constants
    @unpack h, H, Qh = integ.cache

    # Assert EKF0
    PI = InvPrecond(dt)
    @assert all(H .== E1 * PI)

    # More safety checks
    @assert isdiag(H*Qh*H')
    @assert length(unique(diag(H*Qh*H'))) == 1
    Q0_11 = diag(H*Qh*H')[1]

    Σ_ii = h .^ 2 ./ Q0_11
    Σ_ii .= max.(Σ_ii, eps(eltype(Σ_ii)))
    Σ = Diagonal(Σ_ii)

    Σ_out = kron(ones(q+1, q+1), Σ)
    return Σ_out
end


struct MVMLESigma <: AbstractStaticSigmaRule end
initial_sigma(sigmarule::MVMLESigma, d, q) = kron(ones(q+1, q+1), diagm(0 => ones(d)))
function static_sigma_estimation(kind::MVMLESigma, integ)
    @unpack dt = integ
    @unpack d, q, R, InvPrecond, E1 = integ.constants
    @unpack measurement, H = integ.cache

    # Assert EKF0
    PI = InvPrecond(dt)
    @assert all(H .== E1 * PI)

    @unpack measurement = integ.cache
    v, S = measurement.μ, measurement.Σ

    # More safety checks
    @assert isdiag(S)
    @assert length(unique(diag(S))) == 1
    S_11 = diag(S)[1]

    Σ_ii = v .^ 2 ./ S_11
    Σ = Diagonal(Σ_ii)
    Σ_out = kron(ones(q+1, q+1), Σ)
    # @info "MV-MLE-Sigma" v S Σ Σ_out

    if integ.success_iter == 0
        @assert length(integ.sigmas) == 0
        return Σ_out
    else
        @assert length(integ.sigmas) == integ.success_iter
        sigma_prev = integ.sigmas[end]
        sigma = sigma_prev + (Σ_out - sigma_prev) / integ.success_iter
        return sigma
    end
end
