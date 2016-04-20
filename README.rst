# snpmatch

snpmatch is a simple library to compare the given SNPs to that of database SNP matrix and identify the right accession. It is used ot genotype a sample from the low-coverage sequencing data. 
It calculates a likelihood score with a each accession and performs a likelihood ratio test of top accession with rest. 

Requirements:
The snpmatch uses various python packages (numpy, pandas, pygwas, vcfnp)
The main SNP database should be a HDF5 file with specific keys. Mainly read using pygwas package. 



