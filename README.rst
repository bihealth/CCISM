.. image:: https://img.shields.io/pypi/pyversions/CCISM.svg
       :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/CCISM.svg
       :target: https://pypi.python.org/pypi/CCISM

.. image:: https://api.codacy.com/project/badge/Grade/9ee0ec1424c143dfad9977a649f917f7
       :target: https://www.codacy.com/app/bihealth/CCISM?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/CCISM&amp;utm_campaign=Badge_Grade

.. image:: https://api.codacy.com/project/badge/Coverage/9ee0ec1424c143dfad9977a649f917f7
       :target: https://www.codacy.com/app/bihealth/CCISM?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/CCISM&amp;utm_campaign=Badge_Coverage

.. image:: https://travis-ci.org/bihealth/CCISM.svg?branch=master
       :target: https://travis-ci.org/bihealth/CCISM


=======
CCISM
=======

CCISM is a python tool to determine tumor cells using expressed somatic variants

---------------
Installation
---------------

the recommended way is to create a conda enviroment containing CCISM 

.. code-block:: shell

	$ conda create -n CCISM python=3.7
	$ conda activate CCISM
        $ pip install git+https://github.com/bihealth/CCISM.git
	
------------------------
Running CCISM
------------------------

use `cellsnp-lite <https://github.com/single-cell-genetics/cellsnp-lite>`_ to count expressed somatic variants in a bam file with reads from your single-cell experiment:

.. code-block:: shell

	$ cellsnp-lite -s $BAM_FILE -b $BARCODE_FILE -O $CELLSNP_OUT -R $SOMATIC_VARIANT_VCF -p 8 --minMAF 0.001 --minCOUNT 1 --gzip


this will create ``$CELLSNP_OUT`` with the following contents

.. code-block:: shell

	$ ls $CELLSNP_OUT
	cellSNP.base.vcf.gz  cellSNP.samples.tsv  cellSNP.tag.AD.mtx  cellSNP.tag.DP.mtx  cellSNP.tag.OTH.mtx


next, run CCISM like this

.. code-block:: shell

	$ CCISM -i $CELLSNP_OUT -o $SCITCEM_OUT

which will create ``$SCITCEM_OUT`` with a ``results.txt`` file and a ``parameters.txt`` file.

--------------------
additional arguments
--------------------

``--min_counts``
    filter cells by min number of SNV-covering counts (default: 1)

``--thetaT``
    initial estimate for allelic fraction in tumor cells (default: 0.4)

``--thetaN``
    fixed allelic fraction in normal cells (default 1.e-4)

``--use_vireo``
    use `vireo <https://github.com/single-cell-genetics/vireo>`_ instead (requires ``vireoSNP`` to be installed)

``--estimate_power``
    use simulations to estimate statistical power (use ``--n`` to set the number of replicates and ``--frac_tumor`` to set the fraction of tumor cells)
    
--------------------
Understanding output
--------------------

``results.txt``
    tab-separated file with one row for each cell barcode and the following columns
    
    1. p: probability that this cell is a tumor cell
       
    2. nSNPs_tot: tot # of SNVs that were detected in this cell
       
    3. nSNPs_ref: # of SNVs with reference allele
       
    4. nSNPs_alt: # of SNVs with alternative allele
       
    5. nUMI_tot: tot UMI count for SNVs in this cell
       
    6. nUMI_ref: UMI count for SNVs with reference allele
    
    7. nUMI_alt: UMI count for SNVs with alternative allele
       
    8. ref_SNPs: comma-separated list of SNVs found with reference allele
       
    9. alt_SNPs: comma-separated list of SNVs found with alternative allele


``parameters.txt``
   summarises final estimates for log likelihood, ``thetaT`` and ``thetaN``.

``power_estimates.txt`` (optional)
   contains estimated true positive and false positive rates (TPR / FPR) in simulations
