#!/usr/bin/env python3
"""Exercises 5 and 6 -- YOUR version.

Fill in the four TODOs, then run

    apptainer exec --bind "$PWD/.." "$PADM_SIF" \
        python3 docs/exercises/surrogate_skeleton.py

Each unfinished TODO stops the script with a message saying what it needs.
No CFD is run: everything here uses the twelve designs the campaign already
paid for.  The reference answer is surrogate.py next to this file -- try
without it first.
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
CAMPAIGN = STUDY / "bayes_optimize_sequential.yaml"
NAMES = ["a_weak", "a_strong_ratio", "t_s", "t_m_ratio", "L_c", "L_s_ratio"]


def todo(what: str):
    raise NotImplementedError(f"TODO -- {what}")


def main() -> None:
    box = yaml.safe_load(CAMPAIGN.read_text())["bo_parameters"]
    lower = torch.tensor([float(box[n]["lower"]) for n in NAMES])
    upper = torch.tensor([float(box[n]["upper"]) for n in NAMES])

    rows = [r for r in csv.DictReader(open(SAMPLES))
            if r.get("metric_flux_weighted_mixing_index", "")]
    X_phys = torch.tensor([[float(r["bo_" + n]) for n in NAMES] for r in rows])
    Y = torch.tensor([[float(r["metric_flux_weighted_mixing_index"])] for r in rows])
    print(f"{len(rows)} completed designs, {X_phys.shape[1]} parameters\n")

    # TODO 1 -- the bo_* columns are PHYSICAL values with very different ranges
    # (look at `lower` and `upper`).  Map them to the unit box [0, 1]^6, so that
    # each length scale you read later is in units of that parameter's range.
    X = todo("1: X = map X_phys to the unit box using lower and upper")

    # TODO 2 -- build a SingleTaskGP with an ARD Matern-5/2 kernel (one length
    # scale per parameter, wrapped in a ScaleKernel), Standardize the outcome,
    # and fit it by maximizing the exact marginal log likelihood.
    model = todo("2: build and fit the GP")

    # EXERCISE 5 -----------------------------------------------------------
    # TODO 3 -- read the six length scales from the fitted kernel and print them
    # sorted, shortest first.  Which parameters does the objective depend on?
    todo("3: print the ARD length scales, shortest first")

    # EXERCISE 6 -----------------------------------------------------------
    # TODO 4 -- for kappa in 0, 2, 10, maximize UpperConfidenceBound (it takes
    # beta = kappa**2) over the unit box with optimize_acqf, map the candidate
    # back to physical values, and print it.  Which coordinates move as kappa
    # grows, and why those?
    todo("4: maximize UCB for kappa = 0, 2, 10 and print the designs")


if __name__ == "__main__":
    main()
