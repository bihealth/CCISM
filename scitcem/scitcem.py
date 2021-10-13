import os
import sys
import scipy.sparse
import scipy.special
import scipy.io
import numpy as np
import pandas as pd


def logchoose(n, k):
    """helper function logchoose using scipy.special.loggamma"""

    return (
        scipy.special.loggamma(n + 1)
        - scipy.special.loggamma(n - k + 1)
        - scipy.special.loggamma(k + 1)
    )


def get_logbc(A, D):
    """helper function to get constant part of log likelihood"""
    nz = np.nonzero(D)
    tmp = D.copy()
    tmp[nz] = logchoose(D[nz], A[nz])
    return np.sum(tmp, axis=0).A1


def get_logP(A, D, th, logbc=None):
    """helper function to get logP given matrix A, D and parameter th"""

    if logbc is None:
        logbc = get_logbc(A, D)

    logP = (
        logbc
        + np.sum(A, axis=0).A1 * np.log(th)
        + np.sum(D - A, axis=0).A1 * np.log(1 - th)
    )

    return logP


def E_step(A, D, th, logbc=None):
    """E step in expectation maximization"""

    eps = np.finfo(float).eps

    log_pn = get_logP(A, D, th["normal"], logbc)
    log_pt = get_logP(A, D, th["tumor"], logbc)

    diff = np.minimum(
        -np.log(eps), np.maximum(np.abs(log_pn - log_pt), -np.log(1.0 - eps))
    )
    sign = np.sign(log_pn - log_pt)

    return 1.0 / (1.0 + np.exp(sign * diff))


def M_step(A, D, p, th_N):
    """M step in expectation maximization"""

    res = {}
    rsA = np.sum(A, axis=0).A1
    rsD = np.sum(D, axis=0).A1
    res["tumor"] = np.sum(p * rsA) / np.sum(p * rsD)
    if th_N is None:
        res["normal"] = np.sum((1 - p) * rsA) / np.sum((1 - p) * rsD)
    else:
        res["normal"] = th_N

    return res


def get_logL(A, D, th, logbc=None):
    """get log likelihood"""
    eps = np.finfo(float).eps

    log_pn = get_logP(A, D, th["normal"], logbc)
    log_pt = get_logP(A, D, th["tumor"], logbc)

    maxlogp = np.minimum(
        np.log(1.0 - eps), np.maximum(np.maximum(log_pn, log_pt), np.log(eps))
    )
    minlogp = np.minimum(
        np.log(1.0 - eps), np.maximum(np.minimum(log_pn, log_pt), np.log(eps))
    )
    return np.sum(maxlogp + np.log(1.0 + np.exp(minlogp - maxlogp)))


def EMoptimize(
    A, D, thT_0, thN_0, verbose=False, max_diff=1.0e-5, fit_normal=False
):
    """run expectation maximization"""
    logbc = get_logbc(A, D)

    th = {"normal": thN_0, "tumor": thT_0}

    logL = 0
    for i in range(10):

        p = E_step(A, D, th, logbc=logbc)
        if fit_normal:
            th_new = M_step(A, D, p)
        else:
            th_new = M_step(A, D, p, th["normal"])

        logL_new = get_logL(A, D, th_new, logbc=logbc)

        if np.isnan(logL_new) or np.isinf(logL_new):
            raise Exception("invalid value in logL")

        diff = np.sum(np.abs(logL_new - logL)) / (int(i == 0) + np.abs(logL))
        th = th_new
        logL = logL_new
        if verbose and i > 0:
            sys.stderr.write(
                "iter {0}, diff={1:.2g}, logL={2:.2g}\n".format(i, diff, logL)
            )

        if i > 1 and diff < max_diff:
            break

    if i < 100:
        if verbose:
            sys.stderr.write(
                "optimization done. th_N={0:.2g}, th_T={1:.2g}\n".format(
                    th["normal"], th["tumor"]
                )
            )
    else:
        sys.stderr.write("max number of iterations exceeded!\n")

    return {
        "p": p,
        "logL": logL,
        "thetaT": th["tumor"],
        "thetaN": th["normal"],
        "method": "scitcem",
    }


def vireo_optimize(A, D, verbose):
    """run vireo BinomMixtureVB with two donors"""
    import vireoSNP

    if verbose:
        sys.stderr.write("using vireoSNP.BinomMixtureVB\n")
    vo = vireoSNP.BinomMixtureVB(
        n_var=D.shape[0], n_cell=D.shape[1], n_donor=2
    )
    vo.fit(A, D)
    # determine which donor has higher variant allele frequency
    cells_0 = np.argmax(vo.ID_prob, axis=1) == 0
    cells_1 = np.argmax(vo.ID_prob, axis=1) == 1
    vaf_0 = A[:, cells_0].sum() / D[:, cells_0].sum()
    vaf_1 = A[:, cells_1].sum() / D[:, cells_1].sum()
    if vaf_0 > vaf_1:
        p = vo.ID_prob[:, 0]
    else:
        p = vo.ID_prob[:, 1]

    return {"p": p, "method": "vireo"}


def simulate_A(D, frac_tumor, thetaT, thetaN, seed=None):
    """get simulated matrix A given matrix D, tumor fraction,
       and parameters thetaT and thetaN"""
    if seed is not None:
        np.random.seed(seed)

    ncells = D.shape[1]
    cells = np.arange(ncells)

    tumor = np.random.rand(ncells) < frac_tumor

    st, ct = D[:, tumor].nonzero()
    At = np.random.binomial(D[:, tumor][(st, ct)].A1.astype(int), thetaT)
    sn, cn = D[:, ~tumor].nonzero()
    An = np.random.binomial(D[:, ~tumor][(sn, cn)].A1.astype(int), thetaN)

    A = scipy.sparse.csr_matrix(
        (
            np.concatenate([At, An]),
            (
                np.concatenate([st, sn]),
                np.concatenate([cells[tumor][ct], cells[~tumor][cn]]),
            ),
        ),
        shape=D.shape,
    )

    return A, tumor


def read_cellSNP(indir):
    """read cellSNP output"""

    assert os.path.isdir(indir), (
        "input directory " + indir + " does not exist!"
    )

    cell_file = os.path.join(indir, "cellSNP.samples.tsv")
    assert os.path.isfile(cell_file), (
        "cellSNP.samples.tsv missing from " + indir
    )
    cells = pd.read_csv(cell_file, sep="\t", squeeze=True, header=None).values

    vcf_file = os.path.join(indir, "cellSNP.base.vcf")
    assert os.path.isfile(vcf_file) or os.path.isfile(vcf_file + ".gz"), (
        "cellSNP.base.vcf(.gz.) file missing from " + indir
    )
    if os.path.isfile(vcf_file):
        vcf = pd.read_csv(
            vcf_file, skiprows=1, sep="\t", header=0, index_col=None
        )
    elif os.path.isfile(vcf_file + ".gz"):
        vcf = pd.read_csv(
            vcf_file + ".gz", skiprows=1, sep="\t", header=0, index_col=None
        )

    SNPs = vcf.apply(
        lambda x: str(x["#CHROM"])
        + ":"
        + str(x["POS"])
        + x["REF"]
        + ">"
        + x["ALT"],
        axis=1,
    ).values

    A_file = os.path.join(indir, "cellSNP.tag.AD.mtx")
    assert os.path.isfile(A_file), "cellSNP.tag.AD.mtx missing from " + indir
    A = scipy.io.mmread(A_file).tocsr()

    D_file = os.path.join(indir, "cellSNP.tag.DP.mtx")
    assert os.path.isfile(D_file), "cellSNP.tag.DP.mtx missing from " + indir
    D = scipy.io.mmread(D_file).tocsr()

    return {"cells": cells, "SNPs": SNPs, "A": A, "D": D}


def run_scitcem(
    indir,
    outdir,
    min_counts,
    thetaT,
    thetaN,
    use_vireo,
    estimate_power,
    nrep,
    frac_tumor,
    verbose,
):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    data = read_cellSNP(indir)

    take_SNPs = data["D"].sum(1).A1 > 0
    take_cells = data["D"].sum(0).A1 >= min_counts

    R = data["D"] - data["A"]
    ref_SNPs = [
        ",".join(data["SNPs"][data["D"][:, i].nonzero()[0]])
        for i in range(R.shape[1])
    ]
    alt_SNPs = [
        ",".join(data["SNPs"][data["A"][:, i].nonzero()[0]])
        for i in range(R.shape[1])
    ]

    df = pd.DataFrame(
        {
            "p": np.nan,
            "nSNPs_tot": (data["D"] > 0).sum(0).A1,
            "nSNPs_ref": (R > 0).sum(0).A1,
            "nSNPs_alt": (data["A"] > 0).sum(0).A1,
            "nUMI_tot": data["D"].sum(0).A1,
            "nUMI_ref": R.sum(0).A1,
            "nUMI_alt": data["A"].sum(0).A1,
            "ref_SNPs": ref_SNPs,
            "alt_SNPs": alt_SNPs,
        },
        index=data["cells"],
    )

    assert (
        take_SNPs.sum() > 0
        and take_cells.sum() > 0
        and data["A"][take_SNPs, :][:, take_cells].sum() > 0
    ), "not enough SNP information in this dataset!"

    A = data["A"][take_SNPs, :][:, take_cells]
    D = data["D"][take_SNPs, :][:, take_cells]

    if use_vireo:
        res = vireo_optimize(A, D, verbose)
    else:
        res = EMoptimize(A, D, thetaT, thetaN, verbose)

    df.loc[take_cells, "p"] = np.round(res["p"], 3)

    df.to_csv(
        os.path.join(outdir, "results.txt"), sep="\t", float_format="%.4f"
    )
    res.pop("p", None)
    pd.Series(res).to_csv(
        os.path.join(outdir, "parameters.txt"),
        sep="\t",
        header=False,
        float_format="%.4f",
    )

    if estimate_power:

        TPR = []
        FPR = []
        if verbose:
            sys.stderr.write("estimating power\n")
        for k in range(nrep):
            At, label = simulate_A(D, frac_tumor, thetaT, thetaN, seed=k)
            res = (
                vireo_optimize(At, D, verbose)
                if use_vireo
                else EMoptimize(At, D, thetaT, thetaN, verbose)
            )
            TP = np.sum(label & (res["p"] > 0.5))
            FP = np.sum(~label & (res["p"] > 0.5))
            TN = np.sum(~label & (res["p"] <= 0.5))
            FN = np.sum(label & (res["p"] <= 0.5))
            TPR.append(TP / (TP + FN))
            FPR.append(FP / (FP + TN))

        pd.DataFrame({"TPR": TPR, "FPR": FPR}).to_csv(
            os.path.join(outdir, "power_estimates.txt"), sep="\t", index=False
        )

    return
