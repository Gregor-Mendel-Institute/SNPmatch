SNPmatch

Genotyping a sample with high-coverage sequencing is expensive, whereas low-coverage is highly error prone. Now with more than 1000 Arabidopsis thaliana high-quality genomes available, genotyping can be simple by identifying the natural lines (accessions) that matches to the sample of interest. Hence we have developed SNPmatch, a tool that identifies the most likely accession with a minimum number of SNPs. This makes it an efficient and reliable tool for low-coverage sequencing data. SNPmatch can identify inbreds but can also be applied to identify the parents in a F1 or F2 population. If you wish to identify your plant you just need to upload your vcf file.

SNPmatch can genotype samples efficiently and economically using a simple likelihood approach. This approach allows SNPmatch to genotype samples with as low as 5000 informative SNPs. SNPmatch is a python script taking about 18 seconds for the entire analysis, to compare 500000 SNPs to the 10.7 million markers in the reference database (1001 SNP matrix) 

Requirements:
The snpmatch uses various python packages (numpy, pandas, pygwas, vcfnp)
The main SNP database should be a HDF5 file with specific keys. Mainly read using pygwas package. 

Input File:
Presently, SNPmatch takes input files with two different formats, bed and vcf. The bed file is a three column file with chr, position and genotype tab-separated. An example is given below

1 125 0/0
1 284 0/0
1 336 0/0
1 346 1/1
1 353 0/0
1 363 0/0
1 465 0/0
1 471 0/1
1 540 0/0
1 564 0/0
1 597 0/0
1 612 1/1
1 617 0/1

The VCF file has a default format detailed in the link below. The main arguments required for SNPmatch are CHROM and POS in header and GT in the INFO column. PL (Normalized Phred-scaled likelihoods of the possible genotypes), if present improves the efficiency of SNPmatch.
http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it

Template files are given in sample_files


SNPmatch as a webtool for A. thaliana, AraGeno is here: http://arageno.gmi.oeaw.ac.at
