# SNPmatch

SNPmatch is a python toolkit which can be used to genotype a sample from as-low-as as 4000 markers from the database lines. SNPmatch can genotype samples efficiently and economically using a simple likelihood approach.

## Installation & Usage

The below steps deal with running SNPmatch on a local machine. This was also tested in Python 3.

### Using pip:

```bash
pip install -e git+https://github.com/Gregor-Mendel-Institute/SNPmatch.git
## or just
pip install SNPmatch
```

### Requirements

The SNPmatch uses various python packages (Cython, numpy, pandas, [PyGWAS](https://github.com/timeu/PyGWAS), vcfnp). The SNP database should be a HDF5 file can be generated with the scripts given in the scripts [folder](https://github.com/Gregor-Mendel-Institute/SNPmatch/tree/master/scripts).

This is read using pygwas package. Database SNPs for the Regmap and 1001genomes dataset for *Arabidopsis thaliana* can be downloaded [here](https://gmioncloud-my.sharepoint.com/personal/uemit_seren_gmi_oeaw_ac_at/_layouts/15/guestaccess.aspx?folderid=0ca806e676c154094992a9e89e5341d43&authkey=AXJPl6GkD8vNPDZJwheb6uk).

### Input files

Database files should be HDF5 file formatted. The database file is generated twice chuked rowwise and column wise in another to increase the efficiency. For *A. thaliana*, users can download the files from the link above.

As the input files, SNPmatch takes two file formats for the markers (BED and VCF). Example input files are given in the folder [sample_files](https://github.com/Gregor-Mendel-Institute/SNPmatch/tree/master/sample_files). Briefly, BED files should be three tab-separated column with chromosome, position and genotype shown below.

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

SNPmatch can be run directly for *A. thaliana* researchers as a web tool, [AraGeno](http://arageno.gmi.oeaw.ac.at)

## Contributing
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D


## History

- 1.7.2: Stable version, 15-12-2016
- 1.8.2: Stable version, 16-02-2017
- 1.9.1: Stable version, 24-08-2017


## Credits

- Rahul Pisupati (rahul.pisupati[at]gmi.oeaw.ac.at)
- Ümit Seren (uemit.seren[at]gmi.oeaw.ac.at)

## Citation

Pisupati et. al. Verification of Arabidopsis stock collection using SNPmatch - an algorithm to genotype high-plexed samples, manuscript under preparation.
