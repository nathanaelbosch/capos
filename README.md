# Calibrated Adaptive Probabilistic ODE Solvers - Code

<p>
  <img src="figures/solution.svg" width="400">
  <img src="figures/samples.svg" width="400">
</p>

This repo contains the code which was used to compute the results of the paper "Calibrated Adaptive Probabilistic ODE Solvers", currently published in
[arXiv](https://arxiv.org/abs/2012.08202).


---

__To solve differential equations in Julia with probabilistic numerical solvers, please use
[ProbNumDiffEq.jl](https://github.com/nathanaelbosch/ProbNumDiffEq.jl)!__<br />
The code in this repository is not meant to be used as generic ODE solvers, whereas
[ProbNumDiffEq.jl](https://github.com/nathanaelbosch/ProbNumDiffEq.jl)
is a Julia package under active development.
It is more stable and documented, its solvers are more efficent, and it contains more features.
The DE solvers it provides are compatible with the
[DifferentialEquations.jl](https://docs.sciml.ai/stable/)
ecosystem.

---

A __Python__ implementation of these solvers, as well as of additional probabilistic numerical methods, is maintained in [ProbNum](https://github.com/probabilistic-numerics/probnum).


## Usage
The experiments can be found in `./experiments`, and the actual solvers are implemented in `./src`.
To start, open a Julia console with `julia --project=.`.
You can then run the experiments with, e.g., `include(experiments/1 Stiff Van der Pol/main.jl)`.


## Reference
```
@InProceedings{pmlr-v130-bosch21a,
  title = 	     { Calibrated Adaptive Probabilistic ODE Solvers },
  author =       {Bosch, Nathanael and Hennig, Philipp and Tronarp, Filip},
  booktitle = 	 {Proceedings of The 24th International Conference on Artificial Intelligence and Statistics},
  pages = 	     {3466--3474},
  year = 	     {2021},
  editor = 	     {Banerjee, Arindam and Fukumizu, Kenji},
  volume = 	     {130},
  series = 	     {Proceedings of Machine Learning Research},
  month = 	     {13--15 Apr},
  publisher =    {PMLR},
  pdf = 	     {http://proceedings.mlr.press/v130/bosch21a/bosch21a.pdf},
  url = 	     {http://proceedings.mlr.press/v130/bosch21a.html},
}
```
