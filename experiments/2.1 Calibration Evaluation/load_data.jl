using ProbNumODE
using BSON
using DataFrames
using GaussianDistributions
using StructArrays

files = [
    "logistic.bson",
    "lotkavolterra.bson",
    "fitzhughnagumo.bson",
    "rigidbody.bson",
    "brusselator.bson",
    "logistic_fixedMLEMV.bson",
    "lotkavolterra_fixedMLEMV.bson",
    "fitzhughnagumo_fixedMLEMV.bson",
    "rigidbody_fixedMLEMV.bson",
    "brusselator_fixedMLEMV.bson",
]
files = "./data/" .* files



# r = results[collect(keys(results))[1]][1]
# config_keys = keys(r.config)
df = DataFrame(
    probname=String[],
    alg=String[],
    q=Int[],
    steprule=Symbol[],
    sigmarule=Symbol[],
    abstol=Real[],
    reltol=Real[],
    dt=Real[],
    smooth=Bool[],
    num_evals=Int[],
    nf=Int[],
    runtime=Real[],
    l∞=Real[],
    l2=Real[],
    L∞=Real[],
    L2=Real[],
    final=Real[],
    χ²=Real[],
    Χ²=Real[],
    final_χ²=Real[],
    final_est=Real[],
    l2_est=Real[],
    L2_est=Real[],
)



for file in files
    results = BSON.load(file)
    for (key, result_list) in results  # Iterate over the individual problems
        # probname = split(splitdir(file)[2], "_")[1]
        probname = splitext(basename(file))[1]
        if '_' in probname
            probname = split(probname, "_")[1]
        end
        for r in result_list
            config = r.config
            summaries = r.summaries
            # ssol = r.solution
            # final_vals = r.final_vals
            # errors = ProbNumODE.compute_errors(ssol)
            errors = summaries

            push!(df, (
                probname,
                config.alg,
                config.q,
                config.steprule,
                config.sigmarule,
                config.abstol,
                config.reltol,
                get(config, :dt, 0.0),
                config.smooth,
                summaries[:num_evals],
                summaries[:nf],
                summaries[:runtime],
                errors[:l∞],
                errors[:l2],
                errors[:L∞],
                errors[:L2],
                errors[:final],
                errors[:χ²],
                errors[:Χ²],
                errors[:final_χ²],
                errors[:final_est],
                errors[:l2_est],
                errors[:L2_est],
            ))
        end
    end
end

categorical!(df, :alg)
categorical!(df, :steprule)
categorical!(df, :sigmarule)

df.sigma_str = "sigma=" .* string.(Symbol.(df.sigmarule))
df.q_str = "IBM" .* string.(df.q)
df.step_str = "step=" .* string.(df.steprule)
df.smooth_str = "smooth=" .* string.(df.smooth)


@show df
