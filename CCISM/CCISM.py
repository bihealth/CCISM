import os
import sys
import scipy.sparse
import scipy.special
import scipy.io
import numpy as np
import pandas as pd

DTYPE = float
EPS = np.finfo(DTYPE).eps


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
    tmp = D.copy().astype(DTYPE)
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

    return logP.astype(DTYPE)


def E_step(A, D, th, logbc=None, check_overflow=True):
    """E step in expectation maximization"""

    log_pn = get_logP(A, D, th["normal"], logbc)
    log_pt = get_logP(A, D, th["tumor"], logbc)

    if check_overflow:
        diff = np.minimum(
            -np.log(EPS),
            np.maximum(np.abs(log_pn - log_pt), -np.log(1.0 - EPS)),
        )
    else:
        diff = np.abs(log_pn - log_pt)
    sign = np.sign(log_pn - log_pt)

    return 1.0 / (1.0 + np.exp(sign * diff))


def M_step(A, D, p, th_N=None):
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


def get_logL(A, D, th, logbc=None, check_overflow=True):
    """get log likelihood"""

    log_pn = get_logP(A, D, th["normal"], logbc)
    log_pt = get_logP(A, D, th["tumor"], logbc)

    if check_overflow:
        maxlogp = np.minimum(
            np.log(1.0 - EPS),
            np.maximum(np.maximum(log_pn, log_pt), np.log(EPS)),
        )
        minlogp = np.minimum(
            np.log(1.0 - EPS),
            np.maximum(np.minimum(log_pn, log_pt), np.log(EPS)),
        )
    else:
        maxlogp = np.maximum(log_pn, log_pt)
        minlogp = np.minimum(log_pn, log_pt)

    return np.sum(maxlogp + np.log(1.0 + np.exp(minlogp - maxlogp)))


def EMoptimize(
        A, D, thT_0, thN_0, verbose=False, max_diff=1.0e-5, fit_normal=False, max_iter=1000, min_iter=10
):
    """run expectation maximization"""
    logbc = get_logbc(A, D)

    th = {"normal": thN_0, "tumor": thT_0}

    logL = 0
    for i in range(max_iter):

        p = E_step(A, D, th, logbc=logbc)
        if fit_normal:
            th_new = M_step(A, D, p)
        else:
            th_new = M_step(A, D, p, th_N=th["normal"])

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

        if i > min_iter and diff < max_diff:
            break

    if i < max_iter:
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
        "method": "CCISM",
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


def read_cellSNP(indir, read_other=False):
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

    if read_other:
        O_file = os.path.join(indir, "cellSNP.tag.OTH.mtx")
        assert os.path.isfile(D_file), "cellSNP.tag.OTH.mtx missing from " + indir
        O = scipy.io.mmread(O_file).tocsr()
    else:
        O = None

    return {"cells": cells, "SNPs": SNPs, "A": A, "D": D, "O": O}


def read_SNVs(SNV_file):

    assert os.path.isfile(SNV_file), "SNV file {0} missing".format(SNV_file)

    if SNV_file.endswith('.vcf') or SNV_file.endswith('.vcf.gz'):
        vcf = pd.read_csv(
            SNV_file, skiprows=1, sep="\t", header=0, index_col=None
        )
        SNVs = vcf.apply(
        lambda x: str(x["#CHROM"])
        + ":"
        + str(x["POS"])
        + x["REF"]
        + ">"
        + x["ALT"],
        axis=1,
    ).values
    else:
        SNVs = np.array([line.strip() for line in open(SNV_file)])

    return SNVs


def run_CCISM(
    indir,
    outdir,
    min_counts,
    thetaT,
    thetaN,
    SNV_file,
    use_vireo,
    estimate_power,
    nrep,
    frac_tumor,
    verbose,
):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if verbose:
        sys.stderr.write('reading cellSNP output from '+indir+'\n')
    data = read_cellSNP(indir)

    take_SNVs = data["D"].sum(1).A1 > 0
    take_cells = data["D"].sum(0).A1 >= min_counts

    if SNV_file is not None:
        if verbose:
            sys.stderr.write('restricting to SNVs in '+SNV_file+'\n')
        use_SNVs = read_SNVs(SNV_file)
        take_SNVs = take_SNVs & (np.isin(data['SNPs'],use_SNVs))

    assert (
        take_SNVs.sum() > 0
        and take_cells.sum() > 0
        and data["A"][take_SNVs, :][:, take_cells].sum() > 0
    ), "not enough SNV information in this dataset!"

    A = data["A"][take_SNVs, :]
    D = data["D"][take_SNVs, :]
    SNVs = data["SNPs"][take_SNVs]

    all_SNVs = [
        ",".join(SNVs[D[:, i].nonzero()[0]])
        for i in range(D.shape[1])
    ]
    alt_SNVs = [
        ",".join(SNVs[A[:, i].nonzero()[0]])
        for i in range(D.shape[1])
    ]

    df = pd.DataFrame(
        {
            "p": np.nan,
            "nSNVs_tot": (D > 0).sum(0).A1,
            "nSNVs_alt": (A > 0).sum(0).A1,
            "nUMI_tot": D.sum(0).A1,
            "nUMI_alt": A.sum(0).A1,
            "all_SNVs": all_SNVs,
            "alt_SNVs": alt_SNVs,
        },
        index=data["cells"],
    )

    if use_vireo:
        res = vireo_optimize(A[:,take_cells],
                             D[:,take_cells],
                             verbose)
    else:
        res = EMoptimize(A[:,take_cells],
                         D[:,take_cells],
                         thetaT, thetaN, verbose)

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
            At, label = simulate_A(D[:, take_cells], frac_tumor, thetaT, thetaN, seed=k)
            res = (
                vireo_optimize(At, D[:, take_cells], verbose)
                if use_vireo
                else EMoptimize(At, D[:, take_cells], thetaT, thetaN, verbose)
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
