# SNPmatch

SNPmatch is a Python toolkit which can be used to genotype a sample from as-low-as as 4000 markers from the database lines. SNPmatch can genotype samples efficiently and economically using a simple likelihood approach.

## Installation & Usage

The below steps deal with running SNPmatch on a local machine. This was also tested in Python 3.

### Dependencies
The SNPmatch uses various python packages (Cython, numpy, pandas, [PyGWAS](https://github.com/timeu/PyGWAS), vcfnp). Most of which are automatically downloaded and installed with pip. Cython has to be installed before hand as given below.
```bash
pip install Cython
```

### Installation using pip:

SNPmatch can be easily installed with the help of pip. It can be installed from the git repo or through PyPi. The requried pip commands are given below.

```bash
## installing SNPmatch from git hub repository
pip install -e git+https://github.com/Gregor-Mendel-Institute/SNPmatch.git
## or PyPi
pip install SNPmatch
```

### Database files

Database files containing the known genotype information for many strains have to be provided as HDF5 formatted file. These can be generated with given markers or variants present in a VCF file. The database files can be generated with the help scripts given in the git [folder](https://github.com/Gregor-Mendel-Institute/SNPmatch/tree/master/scripts). A detailed README is also provided in the folder.
These files are read using PyGWAS package. So we recommend you to generate files as mentioned.

For *Arabidopsis thaliana* users, we have made SNP database files for the `RegMap` and `1001Genomes` panel available and can be downloaded [here](https://gmioncloud-my.sharepoint.com/personal/uemit_seren_gmi_oeaw_ac_at/_layouts/15/guestaccess.aspx?folderid=0ca806e676c154094992a9e89e5341d43&authkey=AXJPl6GkD8vNPDZJwheb6uk).

### Input file

As the input file, SNPmatch takes genotype information in two file formats (BED and VCF). Example input files are given in the folder [sample_files](https://github.com/Gregor-Mendel-Institute/SNPmatch/tree/master/sample_files). Briefly, BED files should be three tab-separated column with chromosome, position and genotype shown below.

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

### Usage

SNPmatch can be run as bash commands given below. A detailed manual for each command with -h.

```bash
snpmatch inbred -i input_file -d db.hdf5 -e db.acc.hdf5 -o output_file
# or
snpmatch parser -i input_file -o input_npz
snpmatch inbred -i input_npz -d db.hdf5 -e db.acc.hdf5 -o output_file
```

### AraGeno

SNPmatch can be run directly for *A. thaliana* researchers as a web tool, [AraGeno](http://arageno.gmi.oeaw.ac.at)

## Genotyping a hybrid

SNPmatch can be used to identify hybrid individuals when parental strains are present in database. For such individuals, SNPmatch can be run in windows across the genome. The commands used to run are given below

```bash
snpmatch cross -d db.hdf5 -e db.acc.hdf5 -i input_file -b window_size_in_bp -o output_file
#to get a genetic map for the hybrid
snpmatch genotype_cross -e db.acc.hdf5 -p parent1xparent2 -i input_file -o output_file
# or if parents have VCF files individually
snpmatch genotype_cross -p parent1.vcf -q parent2.vcf -i input_file -o output_file
```

These scripts are implemented based on the *A. thaliana* genome sizes. But the global variable in csmatch [script](https://github.com/Gregor-Mendel-Institute/SNPmatch/blob/master/snpmatch/core/csmatch.py#L19) can be modified to the corresponding genome sizes.


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
