#!/usr/bin/env python3
"""Exercises 5 and 6 -- reference solution.

Fits a Gaussian process to the designs the campaign has already run, reads the
ARD length scales as a sensitivity analysis, and asks the acquisition function
where design 13 should go.

NO CFD IS RUN.  This is the whole point: the surrogate half of Bayesian
optimization costs about a second on one core, which is why it can live on a
login node while the simulations go to the queue.

    apptainer exec --bind "$PWD/.." ../apptainer/padm.sif \
        python3 docs/exercises/surrogate.py
"""

import csv
import pathlib

import torch
from botorch.acquisition import UpperConfidenceBound
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.models.transforms.input import Normalize
from botorch.models.transforms.outcome import Standardize
from botorch.optim import optimize_acqf
from gpytorch.kernels import MaternKernel, ScaleKernel
from gpytorch.mlls import ExactMarginalLogLikelihood

torch.set_default_dtype(torch.float64)

STUDY = pathlib.Path(__file__).resolve().parents[2]
SAMPLES = STUDY / "results" / "corrected_boundary_v3" / "all_samples.csv"
NAMES = ["a_weak", "a_strong_ratio", "t_s", "t_m_ratio", "L_c", "L_s_ratio"]


def main() -> None:
    rows = [r for r in csv.DictReader(open(SAMPLES))
            if r.get("metric_flux_weighted_mixing_index", "")]

    # The BO coordinates are already the unit box: no scaling needed here.
    X = torch.tensor([[float(r["bo_" + n]) for n in NAMES] for r in rows])
    # Maximise the mixing index, so the GP is fitted to +MI directly.
    Y = torch.tensor([[float(r["metric_flux_weighted_mixing_index"])] for r in rows])
    print(f"{len(rows)} completed designs, {X.shape[1]} parameters\n")

    # The same kernel the campaign uses: ARD Matern-5/2, one length scale per
    # parameter.  ARD is what turns the fit into a sensitivity analysis.
    model = SingleTaskGP(
        X, Y,
        covar_module=ScaleKernel(MaternKernel(nu=2.5, ard_num_dims=X.shape[-1])),
        input_transform=Normalize(d=X.shape[-1]),
        outcome_transform=Standardize(m=1),
    )
    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model))

    # EXERCISE 5 -----------------------------------------------------------
    # A SHORT length scale means the objective changes quickly along that
    # parameter, i.e. it matters.  A length scale far larger than the unit box
    # means the model saw no dependence on it at all.
    lengthscales = model.covar_module.base_kernel.lengthscale.squeeze().detach()
    print("ARD length scales  (small = sensitive, huge = inert):")
    for name, value in sorted(zip(NAMES, lengthscales.tolist()), key=lambda t: t[1]):
        verdict = "drives it" if value < 5 else ("some effect" if value < 100 else "inert")
        print(f"  {name:16s} {value:9.2f}   {verdict}")

    # EXERCISE 6 -----------------------------------------------------------
    # UpperConfidenceBound takes beta = kappa**2.
    bounds = torch.stack([torch.zeros(len(NAMES)), torch.ones(len(NAMES))])
    print("\nWhere should design 13 go?")
    for kappa in (0.0, 2.0, 10.0):
        candidate, _ = optimize_acqf(
            UpperConfidenceBound(model, beta=kappa ** 2),
            bounds=bounds, q=1, num_restarts=16, raw_samples=256,
        )
        values = candidate.squeeze().tolist()
        joined = "  ".join(f"{n}={v:.2f}" for n, v in zip(NAMES, values))
        print(f"  kappa={kappa:5.1f}  ->  {joined}")


if __name__ == "__main__":
    main()
