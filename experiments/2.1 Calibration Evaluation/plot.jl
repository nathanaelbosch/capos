using StatsPlots
using Query
using DataFrames
using Distributions
using LaTeXStrings
include("load_data.jl")

using Plots: PX_PER_INCH
const LINEWIDTH = 6.75*PX_PER_INCH
const COLWIDTH = 3.25*PX_PER_INCH
height(width, ratio=16/9) = width / ratio

dpl = 6.75 * 100
dph = dpl / (16/9)


plot_kwargs = (
    size=(1000, 400),
    markersize=8,
    layout=(2, 4),
    tickfontsize=8,
    guidefontsize=9,
    titlefontsize=10,
    legendfontsize=8,
    legendtitlefontsize=9,
)


####################################################################################################
# Performance: Alg, Order, Smoothing, Sigma
####################################################################################################
for probname in unique(df.probname)
    _df = df |>
        @filter((_.steprule .== :standard)) |>
        @filter((_.probname .== probname)) |>
        @filter((_.L2 .< 1e2)) |>
        DataFrame

    xlims, ylims = Dict(
        "logistic" => ((5e0, 1e6), (1e-13, 1e1),),
        "lotkavolterra" => ((1e2, 3e6), (1e-12, 1e2)),
        "fitzhughnagumo" => ((1e2, 3e6), (1e-13, 1e2)),
        "rigidbody" => ((4e1, 1.5e6), (5e-14, 1e2)),
        "brusselator" => ((4e1, 2.5e6), (1e-12, 2e2),),
    )[probname]

    legendplot = plot(
        repeat([missing], 1, 10),
        label=["EKF-IWP1" "EKF-IWP2" "EKF-IWP3" "EKF-IWP4" "EKF-IWP5" "EKS-IWP1" "EKS-IWP2" "EKS-IWP3" "EKS-IWP4" "EKS-IWP5"],
        color=[1 2 3 4 5 1 2 3 4 5],
        linestyle=hcat(repeat([:solid], 5)..., repeat([:dash], 5)...),
        marker=hcat(repeat([:o], 5)..., repeat([:star], 5)...),
        grid=false, xlims=(20,3), showaxis=false,
        legend=:left,
    );

    # plots = []
    i = 0
    for alg in ["EKF0", "EKF1"],
        sigma in [:fixedMLE, :dynamic, :fixedMLEMV,]
        i+=1

        if sigma in (:fixedMLEMV,) && alg == "EKF1"
            sigma = :dynamicMV
            alg = "EKF0"
        end

        sigmaname = Dict(
            :fixedMLE => "fixed",
            :dynamic => "TV",
            :fixedMLEMV => "fixed-MV",
            :dynamicMV => "TV-MV",
        )[sigma]

        p = _df |>
            @filter((_.sigmarule .== sigma)) |>
            @filter((_.alg .== alg)) |>
            @df StatsPlots.plot(
                :num_evals, :L2,
                group=(:q_str, :smooth_str),
                color=:q,
                linestyle=[:solid, :dash, :dot, :dashdot][:smooth.+1],
                markershape=[:o, :utriangle][:smooth.+1],
                xscale=:log10, yscale=:log10,
                xlabel = (alg == "EKF1" || sigma in (:dynamicMV,)) ? L"\#Evaluations" : " ",
                ylabel = sigma == :fixedMLE ? "L2 Error" : "",
                title="EKF$(alg[end])/EKS$(alg[end]), Γ=$sigmaname",
                legend=false,
            );
        plot!(p, xlims=xlims, ylims=ylims,)
        savefig(p, "/home/nath/PhD-Projects/ode-filter-uq-paper/figures/$(probname)_performance_$i.tex")

    end
    l = @layout [
        [a b c
         d e f] g{0.13w}
    ]
    p = plot(plots..., legendplot;
         xlims=xlims, ylims=ylims,
         plot_kwargs...,
         layout=l,
         );
    pgfsave("figures/$(probname)_performance.pdf")
    savefig("figures/$(probname)_performance.pdf")
end


####################################################################################################
# Calibration
####################################################################################################
for probname in unique(df.probname)
    _df = df |>
        @filter((_.steprule .== :standard)) |>
        @filter((_.probname .== probname)) |>
        @filter((_.L2 .< 1e2)) |>
        DataFrame

    xlims, ylims = Dict(
        "logistic" => ((1e-13, 1e1), (1e-7, 1e3),),
        "lotkavolterra" => ((1e-12, 1e2), (5e-9, 1.1e3)),
        "fitzhughnagumo" => ((1e-13, 1e2), (1e-8, 1e4)),
        "rigidbody" => ((5e-14, 1e2), (8e-7, 2e4)),
        "brusselator" => ((1e-12, 2e2), (5e-9, 1e3),),
    )[probname]

    legendplot = plot(
        repeat([missing], 1, 12),
        label=["EKF-IWP1" "EKF-IWP2" "EKF-IWP3" "EKF-IWP4" "EKF-IWP5" "EKS-IWP1" "EKS-IWP2" "EKS-IWP3" "EKS-IWP4" "EKS-IWP5" "E[χ²(2)]" "99% interval"],
        color=[1 2 3 4 5 1 2 3 4 5 :black :black],
        linestyle=hcat(repeat([:solid], 5)..., repeat([:dash], 5)..., :dash, :dot),
        marker=hcat(repeat([:o], 5)..., repeat([:star], 5)..., :none, :none),
        grid=false, xlims=(20,3), showaxis=false,
        legend=:left,
    )

    plots = []
    for alg in ["EKF0", "EKF1"],
        sigma in [:fixedMLE, :dynamic, :fixedMLEMV]

        if sigma in (:fixedMLEMV,) && alg == "EKF1"
            sigma = :dynamicMV
            alg = "EKF0"
        end

        sigmaname = Dict(
            :fixedMLE => "fixed",
            :dynamic => "TV",
            :fixedMLEMV => "fixed-MV",
            :dynamicMV => "TV-MV",
        )[sigma]

        p = plot()
        _p = 1-0.99
        dist = Chisq(2)
        m, q_low, q_high = mean(dist), quantile(dist, 1-(_p/2)), quantile(dist, _p/2)
        hline!(p, [m], color=:black, linestyle=:dash, label="")
        hline!(p, [q_low, q_high], color=:black, linestyle=:dot, label="")

        _df |>
            @filter((_.sigmarule .== sigma)) |>
            @filter((_.alg .== alg)) |>
            @df StatsPlots.plot!(p,
                # :num_evals, :Χ²,
                :L2, :Χ²,
                group=(:q_str, :smooth_str),
                color=:q,
                linestyle=[:solid, :dash, :dot, :dashdot][:smooth.+1],
                markershape=[:o, :star][:smooth.+1],
                xscale=:log10, yscale=:log10,
                xlabel = (alg == "EKF1" || sigma == :dynamicMV) ? "L2 Error" : " ",
                ylabel = sigma == :fixedMLE ? "χ²-statistic" : "",
                title="EKF$(alg[end])/EKS$(alg[end]), Γ=$sigmaname",
                legend=false,
            )
        push!(plots, p)
    end
    l = @layout [
    [a b c
     d e f] g{0.13w}
    ]
    plot(plots..., legendplot;
         xlims=xlims, ylims=ylims,
         plot_kwargs...,
         layout=l,
         )
    savefig("figures/$(probname)_calibration.pdf")
end
