"""Command line interface for scitcem.

This is the main program entry point for the ``scitcem`` executable
"""

import argparse
from .scitcem import run_scitcem

def run_nocmd(_, parser):
    """No command given, print help and ``exit(1)``."""
    parser.print_help()
    parser.exit(1)

def main(argv=None):
    """Main entry point before parsing command line arguments."""

    parser = argparse.ArgumentParser()

    parser.add_argument('-i','--indir',dest='indir',
                        required=True,
                        help="""input folder (cellSNP output folder)""")
    parser.add_argument('-o','--outdir',dest='outdir',
                        required=True,
                        help="""output directory""")
    parser.add_argument('-m','--min_counts',dest='min_counts',
                        default=3,type=int,
                        help="""min number of SNP-covering counts to include a cell [3]""")
    parser.add_argument('--thetaT',dest='thetaT',
                        default=0.4,
                        help="""initial estimate for theta_T""")
    parser.add_argument('--thetaN',dest='thetaN',
                        default=0.01,
                        help="""fixed estimate for theta_N""")
    parser.add_argument('--use_vireo',dest='use_vireo',
                        action='store_true',default=False)
    parser.add_argument('--verbose',dest='verbose',
                        action='store_true',default=False)
    parser.add_argument("--estimate_power",dest='estimate_power',
                        action='store_true',default=False)
    parser.add_argument("--n",dest='nrep',
                        default=10,type=int,
                        help="""number of replicates for power estimation""")
    parser.add_argument('--frac_tumor',dest='frac_tumor',
                        default=.5,
                        help="""tumor fraction for power estimation""")
 
    args = parser.parse_args()
    run_scitcem(args.indir,
                args.outdir,
                args.min_counts,
                args.thetaT,
                args.thetaN,
                args.use_vireo,
                args.estimate_power,
                args.nrep,
                args.frac_tumor,
                args.verbose)

