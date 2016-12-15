# SNPmatch

SNPmatch is a python toolkit which can be used to genotype a sample from as-low-as as 4000 markers from the database lines. SNPmatch can genotype samples efficiently and economically using a simple likelihood approach.

## Installation & Usage

The below steps deal with running SNPmatch on a local machine.

### Using pip: 

```bash 
pip install -e git+https://github.com/Gregor-Mendel-Institute/SNPmatch.git
## or just
pip install SNPmatch
```

### Requirements

The SNPmatch uses various python packages (numpy, pandas, pygwas, vcfnp). The main SNP database should be a HDF5 file with specific keys. Mainly read using pygwas package.

### Input files

Database SNPs need to be formatted as HDF5 file using pygwas. Currently, SNPmatch takes two file formats, BED and VCF for sample markers. Example input files are given in sample_files folder in git repo. 
Briefly, BED files should be three tab-separated column with chr, position and genotype given below.

```
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
```
VCF file in a default format in the [link](http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it). The main arguments required for SNPmatch are CHROM and POS in header and GT in the INFO column. PL (Normalized Phred-scaled likelihoods of the possible genotypes), if present improves the efficiency of SNPmatch.

### Commands

SNPmatch can be run as bash commands given below. A detailed manual for each command with -h.

```bash
snpmatch inbred -i input_file -d db.hdf5 -e db.acc.hdf5 -o output_file
snpmatch cross -i input_file -d db.hdf5 -e db.acc.hdf5 -o output_file -s score_file
# the SNPmatch parser
snpmatch parser -i intput_file -o output_npz 
```


### AraGENO

SNPmatch can be run directly for A. thaliana researchers as a web tool, [AraGeno](http://arageno.gmi.oeaw.ac.at)

## Contributing
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D


## History

- 1.7.1: Stable version, 14-12-2016


## Credits

- Rahul Pisupati (rahul.pisupati[at]gmi.oeaw.ac.at) 
- Ümit Seren (uemit.seren[at]gmi.oeaw.ac.at)
