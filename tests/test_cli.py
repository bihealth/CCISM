"""Minimal test for CLI"""

import pytest
import tempfile
import os.path
import pandas as pd
import numpy as np

from scitcem import cli, scitcem


def test_run_help():
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cli.main(["--help"])
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


def test_load_data_from_indir():
    indir = os.path.join(os.path.dirname(__file__), "..", "test_data")
    res = scitcem.read_cellSNP(indir)
    assert "cells" in res
    assert "SNPs" in res
    assert "A" in res
    assert "D" in res
    assert res["A"].shape[1] == len(res["cells"])
    assert res["A"].shape[0] == len(res["SNPs"])
    assert res["D"].shape[1] == len(res["cells"])
    assert res["D"].shape[0] == len(res["SNPs"])


def test_run_scitcem():
    indir = os.path.join(os.path.dirname(__file__), "..", "test_data")
    outdir = tempfile.TemporaryDirectory().name
    scitcem.run_scitcem(
        indir, outdir, 3, 0.4, 1.e-4, False, True, 10, 0.5, False
    )

    assert os.path.isfile(os.path.join(outdir, "results.txt"))
    assert os.path.isfile(os.path.join(outdir, "parameters.txt"))
    assert os.path.isfile(os.path.join(outdir, "power_estimates.txt"))

    res = pd.read_csv(
        os.path.join(outdir, "results.txt"), sep="\t", header=0, index_col=0
    )

    res["label"] = res.index.str.rsplit("_", 1).str[-1] == "tumor"

    TP = np.sum(res["label"] & (res["p"] > 0.5))
    FP = np.sum(~res["label"] & (res["p"] > 0.5))
    TN = np.sum(~res["label"] & (res["p"] <= 0.5))
    FN = np.sum(res["label"] & (res["p"] <= 0.5))

    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)

    assert TPR > 0.9
    assert FPR < 0.05

    sim = pd.read_csv(
        os.path.join(outdir, "power_estimates.txt"),
        sep="\t",
        header=0,
        index_col=None,
    ).mean(0)

    assert sim["TPR"] > 0.9
    assert sim["FPR"] < 0.05
