"""Command line interface for CCISM.

This is the main program entry point for the ``CCISM`` executable
"""

import argparse
from .CCISM import run_CCISM


def main(argv=None):
    """Main entry point before parsing command line arguments."""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--verbose", dest="verbose", action="store_true", default=False
    )

    required = parser.add_argument_group("required arguments")
    additional = parser.add_argument_group("additional arguments")

    required.add_argument(
        "-i",
        "--indir",
        dest="indir",
        help="""input folder (cellSNP output folder)""",
    )
    required.add_argument(
        "-o", "--outdir", dest="outdir", help="""output directory"""
    )
    additional.add_argument(
        "-m",
        "--min_counts",
        dest="min_counts",
        default=1,
        type=int,
        help="""min number of SNP-covering counts to include a cell""",
    )
    additional.add_argument(
        "--thetaT",
        dest="thetaT",
        default=0.4,
        type=float,
        help="""initial estimate for theta_T""",
    )
    additional.add_argument(
        "--thetaN",
        dest="thetaN",
        default=0.01,
        type=float,
        help="""fixed estimate for theta_N""",
    )
    additional.add_argument(
        "--use_vireo",
        dest="use_vireo",
        action="store_true",
        default=False,
        help="""use vireo binomial mixture VB with two clones""",
    )
    additional.add_argument(
        "--use_SNVs",
        dest="use_SNVs",
        default=None,
        help="""restrict to SNVs from a file (either vcf or text file with each line like so [chrom]:[pos][ref]>[alt]"""
    )
    additional.add_argument(
        "--estimate_power",
        dest="estimate_power",
        action="store_true",
        default=False,
        help="""use simulations to estimate statistical power""",
    )
    additional.add_argument(
        "--n",
        dest="nrep",
        default=10,
        type=int,
        help="""number of replicates for power estimation""",
    )
    additional.add_argument(
        "--frac_tumor",
        dest="frac_tumor",
        default=0.5,
        type=float,
        help="""tumor fraction for power estimation""",
    )

    args = parser.parse_args()

    if args.indir and args.outdir:
        run_CCISM(
            args.indir,
            args.outdir,
            args.min_counts,
            args.thetaT,
            args.thetaN,
            args.use_SNVs,
            args.use_vireo,
            args.estimate_power,
            args.nrep,
            args.frac_tumor,
            args.verbose,
        )
    else:
        """no input/output directories given, print help and ``exit(1)``."""
        parser.print_help()
        parser.exit(1)
