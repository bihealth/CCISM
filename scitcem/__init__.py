import os
import sys
import gzip
import scipy.sparse
import scipy.special
import scipy.io
import numpy as np
import pandas as pd
from argparse import ArgumentParser

parser=ArgumentParser()
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

def logchoose (n, k):

    return (scipy.special.loggamma(n + 1) - 
            scipy.special.loggamma(n - k + 1) - 
            scipy.special.loggamma(k + 1))

def get_logbc (A, D):
    nz = np.nonzero(D)
    tmp = D.copy()
    tmp[nz] = logchoose(D[nz], A[nz])
    return np.sum(tmp, axis=0).A1

def get_logP (A, D, th, logbc=None):

    # restrict P to [eps, 1-eps] to prevent under/overflow
    # eps = np.finfo(float).eps

    if logbc is None:
        logbc = get_logbc(A, D)

    logP = (logbc + np.sum(A, axis=0).A1*np.log(th) + 
            np.sum(D-A, axis=0).A1*np.log(1-th))

    return logP
    #return np.minimum(np.log(1.-eps),
    #                  np.maximum(np.log(eps), logP))
    

def E_step (A, D, th, logbc=None):

    # restrict abs differences in P to [eps, 1-eps] to prevent under/overflow

    eps = np.finfo(float).eps

    log_pn = get_logP(A, D, th['normal'], logbc)
    log_pt = get_logP(A, D, th['tumor'], logbc)

    diff = np.minimum(-np.log(eps),
                      np.maximum(np.abs(log_pn-log_pt),-np.log(1.-eps)))
    sign = np.sign(log_pn-log_pt)
    
    return 1./(1.+np.exp(sign*diff))
    #return 1./(1.+np.exp(log_pn-log_pt))


def M_step (A, D, p, th_N): 

    res = {}
    rsA = np.sum(A, axis=0).A1
    rsD = np.sum(D, axis=0).A1
    res['tumor'] = np.sum(p*rsA)/np.sum(p*rsD)
    if th_N is None:
        res['normal'] = np.sum((1-p)*rsA)/np.sum((1-p)*rsD)
    else:
        res['normal'] = th_N

    return res


def get_logL (A, D, th, logbc=None):

    # restrict max and min of P to [eps, 1-eps] to prevent under/overflow
    eps = np.finfo(float).eps

    log_pn = get_logP(A, D, th['normal'], logbc)
    log_pt = get_logP(A, D, th['tumor'], logbc)

    maxlogp = np.minimum(np.log(1.-eps),np.maximum(np.maximum(log_pn, log_pt),np.log(eps)))
    minlogp = np.minimum(np.log(1.-eps),np.maximum(np.minimum(log_pn, log_pt),np.log(eps)))
    return np.sum(maxlogp+np.log(1.+np.exp(minlogp-maxlogp)))
    #return np.sum(np.log(np.exp(log_pn) + np.exp(log_pt)))


def EMoptimize (A, D, thT_0, thN_0, 
                verbose=0, max_diff=1.e-5, fit_normal=False):

    logbc = get_logbc(A, D)

    th = {'normal': thN_0,
          'tumor': thT_0}

    logL = 0
    for i in range(10):
        
        p = E_step (A, D, th, logbc=logbc)
        if fit_normal:
            th_new = M_step(A, D, p)
        else:
            th_new = M_step(A, D, p, th['normal'])

        logL_new = get_logL (A, D, th_new, logbc=logbc)

        if np.isnan(logL_new) or np.isinf(logL_new):
            raise Exception("invalid value in logL")

        diff = np.sum(np.abs(logL_new-logL))/(int(i==0) + np.abs(logL))
        th = th_new
        logL = logL_new
        if verbose > 0 and i > 0:
            sys.stderr.write("iter {0}, diff={1:.2g}, logL={2:.2g}\n".format(i,diff, logL))

        if i > 1 and diff < max_diff:
            break

    if i < 100:

        if verbose > 0:
            sys.stderr.write("optimization done. th_N={0:.2g}, th_T={1:.2g}\n".format(th['normal'],th['tumor']))
        return p

    else:

        sys.stderr.write("max number of iterations exceeded. exiting ...\n")
        return np.empty(D.shape[1])

def vireo_optimize (A, D, verbose):
    import vireoSNP
    if (verbose):
        sys.stderr.write('using vireoSNP.BinomMixtureVB\n')
    vo = vireoSNP.BinomMixtureVB(n_var=D.shape[0],
                                 n_cell=D.shape[1],
                                 n_donor=2)
    vo.fit(A, D)
    # determine which donor has higher variant allele frequency
    cells_0 = np.argmax(vo.ID_prob,axis=1)==0
    cells_1 = np.argmax(vo.ID_prob,axis=1)==1
    vaf_0 = A[:,cells_0].sum()/D[:,cells_0].sum()
    vaf_1 = A[:,cells_1].sum()/D[:,cells_1].sum()
    if vaf_0 > vaf_1:
        p = vo.ID_prob[:,0]
    else:
        p = vo.ID_prob[:,1]

    return p
                

def simulate_A (D, frac_tumor, thetaT, thetaN, seed=None):

    if seed is not None:
        np.random.seed(seed)

    ncells = D.shape[1]
    cells = np.arange(ncells)

    tumor = np.random.rand(ncells) < frac_tumor

    st,ct = D[:,tumor].nonzero()
    At = np.random.binomial(D[:,tumor][(st,ct)].A1.astype(int), thetaT)
    sn,cn = D[:,~tumor].nonzero()
    An = np.random.binomial(D[:,~tumor][(sn,cn)].A1.astype(int), thetaN)

    A = scipy.sparse.csr_matrix((np.concatenate([At,An]),
                                 (np.concatenate([st,sn]),
                                  np.concatenate([cells[tumor][ct],
                                                  cells[~tumor][cn]]))),
                                shape=D.shape)

    return A, tumor
    

def read_cellSNP (indir):

    cell_file = os.path.join(indir,'cellSNP.samples.tsv')
    cells = pd.read_csv(cell_file,sep='\t',squeeze=True,header=None).values

    vcf_file = os.path.join(indir,'cellSNP.base.vcf')
    if os.path.isfile(vcf_file):
        vcf = pd.read_csv(vcf_file,skiprows=1,sep='\t',header=0,index_col=None)
    elif os.path.isfile(vcf_file+'.gz'):
        vcf = pd.read_csv(vcf_file+'.gz',skiprows=1,sep='\t',header=0,index_col=None)
    else:
        raise Exception("cannot open vcf file!")
    
    SNPs = vcf.apply(lambda x: str(x['#CHROM']) + 
                     ':g' + str(x['POS']) + x['REF'] + 
                     '>' + x['ALT'], axis=1).values

    A_file = os.path.join(indir,'cellSNP.tag.AD.mtx')
    A = scipy.io.mmread(A_file).tocsr()

    D_file = os.path.join(indir,'cellSNP.tag.DP.mtx')
    D = scipy.io.mmread(D_file).tocsr()

    return {'cells': cells,
            'SNPs': SNPs,
            'A': A,
            'D': D}

def main(indir, output, min_counts, thetaT, thetaN, use_vireo, estimate_power, nrep,
         frac_tumor, dataset_label, verbose):

    data = read_cellSNP(indir)

    take_SNPs = data['D'].sum(1).A1 > 0
    take_cells = data['D'].sum(0).A1 >= min_counts

    R = data['D']-data['A']
    ref_SNPs = [','.join(data['SNPs'][data['D'][:,i].nonzero()[0]]) for i in range(R.shape[1])]
    alt_SNPs = [','.join(data['SNPs'][data['A'][:,i].nonzero()[0]]) for i in range(R.shape[1])]

    df = pd.DataFrame({'p': np.nan,
                       'nSNPs_tot': (data['D'] > 0).sum(0).A1,
                       'nSNPs_ref': (R > 0).sum(0).A1,
                       'nSNPs_alt': (data['A'] > 0).sum(0).A1,
                       'nUMI_tot': data['D'].sum(0).A1,
                       'nUMI_ref': R.sum(0).A1,
                       'nUMI_alt': data['A'].sum(0).A1,
                       'ref_SNPs': ref_SNPs,
                       'alt_SNPs': alt_SNPs},
                      index=data['cells'])

    if take_SNPs.sum() > 0 and take_cells.sum() > 0 and data['A'][take_SNPs,:][:,take_cells].sum() > 0:

        A = data['A'][take_SNPs,:][:,take_cells]
        D = data['D'][take_SNPs,:][:,take_cells]

        if use_vireo:
            p = vireo_optimize(A, D, verbose)
        else:
            p = EMoptimize(A, D, thetaT, thetaN, verbose)

        df.loc[take_cells,'p'] = np.round(p,3)

    else:

        sys.stderr.write("not enough SNP information in this sample!\n")

    df.to_csv(output if output else sys.stdout, sep='\t')

    if args.estimate_power:

        TPR = []
        FPR = []
        if (verbose > 0):
            sys.stderr.write('estimating power\n')
        sys.stdout.write('label\tTPR\tFPR\n')
        for k in range(nrep):
            At, label = simulate_A (D, frac_tumor, thetaT, thetaN, seed=k)
            p = vireo_optimize(At, D, verbose) if use_vireo else EMoptimize(At, D, thetaT, thetaN, verbose)
            TP = np.sum(label & (p > .5))
            FP = np.sum(~label & (p > .5))
            TN = np.sum(~label & (p <= .5))
            FN = np.sum(label & (p <= .5))
            sys.stdout.write('{0}\t{1:.3g}\t{2:.3g}\n'.format(dataset_label,TP/(TP+FN),FP/(FP+TN)))

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.indir,
         args.out,
         args.min_counts,
         args.thetaT,
         args.thetaN,
         args.use_vireo,
         args.estimate_power,
         args.nrep,
         args.frac_tumor,
         args.label,
         args.verbose)
