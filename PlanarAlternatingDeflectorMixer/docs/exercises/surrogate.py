#!/usr/bin/env python3
"""Exercises 5 and 6 -- reference solution.

Fits a Gaussian process to the designs the campaign has already run, reads the
ARD length scales as a sensitivity analysis, and asks the acquisition function
where design 13 should go.

NO CFD IS RUN.  This is the whole point: the surrogate half of Bayesian
optimization costs about a second on one core, which is why it can live on a
login node while the simulations go to the queue.

    apptainer exec --bind "$PWD/.." "$PADM_SIF" \
        python3 docs/exercises/surrogate.py

The participants' version, with the steps left as TODOs, is
surrogate_skeleton.py next to this file.
"""

import csv
import pathlib

import torch
import yaml
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
SAMPLES = STUDY / "reference" / "corrected_boundary_v3" / "all_samples.csv"
# The archived designs came from the research campaign, so its box applies.
CAMPAIGN = STUDY / "bayes_optimize_sequential.yaml"
NAMES = ["a_weak", "a_strong_ratio", "t_s", "t_m_ratio", "L_c", "L_s_ratio"]


def main() -> None:
    box = yaml.safe_load(CAMPAIGN.read_text())["bo_parameters"]
    lower = torch.tensor([float(box[n]["lower"]) for n in NAMES])
    upper = torch.tensor([float(box[n]["upper"]) for n in NAMES])

    rows = [r for r in csv.DictReader(open(SAMPLES))
            if r.get("metric_flux_weighted_mixing_index", "")]

    # The bo_* columns are PHYSICAL values: a_weak in [0, 0.12] H, t_s in
    # [0.04, 0.15] H, L_c in [0.4, 2.4] H, and three ratios in [0, 1].  Map them
    # to the unit box, exactly as the campaign's own driver does, so that every
    # length scale below is measured in units of that parameter's full range.
    X_phys = torch.tensor([[float(r["bo_" + n]) for n in NAMES] for r in rows])
    X = (X_phys - lower) / (upper - lower)
    # Maximize the mixing index, so the GP is fitted to +MI directly.
    Y = torch.tensor([[float(r["metric_flux_weighted_mixing_index"])] for r in rows])
    print(f"{len(rows)} completed designs, {X.shape[1]} parameters")
    print("(12 points for 9 hyperparameters: read what follows as indicative.)\n")

    # The same kernel the campaign uses: ARD Matern-5/2, one length scale per
    # parameter.  ARD is what turns the fit into a sensitivity analysis.
    unit_box = torch.stack([torch.zeros(len(NAMES)), torch.ones(len(NAMES))])
    model = SingleTaskGP(
        X, Y,
        covar_module=ScaleKernel(MaternKernel(nu=2.5, ard_num_dims=X.shape[-1])),
        input_transform=Normalize(d=X.shape[-1], bounds=unit_box),
        outcome_transform=Standardize(m=1),
    )
    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model))

    # EXERCISE 5 -----------------------------------------------------------
    # Length scales are in units of each parameter's range.  One of a few
    # ranges or less means the objective changes noticeably across the box we
    # search: that parameter matters.  Tens of ranges means the model saw
    # essentially no dependence on it.  The labels are rules of thumb.
    lengthscales = model.covar_module.base_kernel.lengthscale.squeeze().detach()
    print("ARD length scales, in units of each parameter's range")
    print("(short = sensitive, long = inert):")
    for name, value in sorted(zip(NAMES, lengthscales.tolist()), key=lambda t: t[1]):
        verdict = "drives it" if value < 4 else ("some effect" if value < 30 else "inert")
        print(f"  {name:16s} {value:9.2f}   {verdict}")

    # EXERCISE 6 -----------------------------------------------------------
    # UpperConfidenceBound takes beta = kappa**2.  Optimize over the unit box,
    # then map the answer back to the physical values a CAD file needs.
    print("\nWhere should design 13 go?  (physical values; * = at a bound of the box)")
    for kappa in (0.0, 2.0, 10.0):
        candidate, _ = optimize_acqf(
            UpperConfidenceBound(model, beta=kappa ** 2),
            bounds=unit_box, q=1, num_restarts=16, raw_samples=256,
        )
        unit = candidate.squeeze()
        phys = lower + unit * (upper - lower)
        joined = "  ".join(
            f"{n}={v:.3f}{'*' if u < 1e-3 or u > 1 - 1e-3 else ' '}"
            for n, v, u in zip(NAMES, phys.tolist(), unit.tolist())
        )
        print(f"  kappa={kappa:5.1f}  ->  {joined}")


if __name__ == "__main__":
    main()
