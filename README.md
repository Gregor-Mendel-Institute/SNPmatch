# SNPmatch

SNPmatch is a Python toolkit which can be used to genotype a sample from as-low-as as 4000 markers from the database lines. SNPmatch can genotype samples efficiently and economically using a simple likelihood approach.

## Installation & Usage

The below steps deal with running SNPmatch on a local machine. This package is only tested in Python 2. Also consider using Nextflow when implementing it on your cluster, we have provided best practices scripts [here](https://github.com/rbpisupati/nf-snpmatch).

### Installation using pip

SNPmatch can be easily installed with the help of pip. SNPmatch uses various python packages (NumPy, pandas, [PyGWAS](https://github.com/timeu/PyGWAS), [scikit-allel](https://github.com/cggh/scikit-allel)), which are automatically downloaded and installed while using pip. Follow the commands below for successful installation.

```bash
## installing SNPmatch from git hub repository
pip install git+https://github.com/Gregor-Mendel-Institute/SNPmatch.git
## or PyPi
pip install SNPmatch
```
SNPmatch can be installed either from the git repo or through PyPi. In case of installation errors, please install these dependencies using the commands below (for a Debian based system).
```bash
sudo apt-get install python-dev libfreetype6-dev libxft-dev libblas-dev liblapack-dev libatlas-base-dev libhdf5-dev gfortran
sudo pip install NumPy
```
Mac users can install these packages using [Homebrew](https://brew.sh/). These packages should be enough to install SNPmatch correctly. Please raise an issue in the Github repo if you still have trouble installing.

### Database files

Database files containing the known genotype information for many strains have to be provided as HDF5 formatted file. These can be generated with given markers or variants present in a VCF file. The database files can be generated with the functions given in SNPmatch. They are generated using the commands given below.

The below commands require BCFtools executable in the path environment. The database files are read using PyGWAS package. So the VCF files need to have biallelic SNPs only for now.

```bash
snpmatch makedb -i input_database.vcf -o db
```

The above command generates three files,
  * db.csv
  * db.hdf5
  * db.acc.hdf5

The two hdf5 files are the main database files used for further analysis. The files have the same information but are chunked for better efficiency. The files db.hdf5 and db.acc.hdf5 are given to the SNPmatch command under -d and -e options respectively.

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

### Output files for inbred

SNPmatch outputs two file,
1. `output_file.scores.txt` --- tab-separated file


|1|2|3|4|5|6|7|8|
|---|---|---|---|---|---|---|---|
|8426|4946|4987|0.99|517.57|1.0|5525|4.55|
|8427|4861|5194|0.93|4897.21|9.46|5525|4.55|
|6191|4368|4933|0.88|8652.07|16.72|5525|4.55|

The column information for the table is strain ID, Number of matched SNPs, Total informative SNPs, Probability of match, Likelihood, Likelihood ratio against best hit, Number of SNPs, Average depth of SNPs respectively.
You can filter the strains that do not match the sample by the likelihood ratio value which is chi-square distributed.

2. `output_file.matches.json` --- JSON file
It might be easier to parse this file using [json editor](https://docs.python.org/2/library/json.html). All the strains that might match to the sample is present in this file. Also a simple interpretation based on the SNPmatch result is provided.

## Genotyping a hybrid

SNPmatch can be used to identify hybrid individuals when parental strains are present in database. For such individuals, SNPmatch can be run in windows across the genome. The commands used to run are given below

```bash
snpmatch cross -d db.hdf5 -e db.acc.hdf5 -i input_file -b window_size_in_bp -o output_file
#to identify the windows matching to each parent in a hybrid
snpmatch genotype_cross -e db.acc.hdf5 -p parent1xparent2 -i input_file -o output_file
# or if parents have VCF files individually
snpmatch genotype_cross -p parent1.vcf -q parent2.vcf -i input_file -o output_file
```

These scripts are implemented based on the *A. thaliana* genome sizes. But the global variable in csmatch [script](https://github.com/Gregor-Mendel-Institute/SNPmatch/blob/master/snpmatch/core/csmatch.py#L19) can be modified to the corresponding genome sizes.

### Output files for cross
We have three output files for the cross match in SNPmatch
1. `output_file.scores.txt` --- tab-separated file
The file is exactly same as explained [before](#Output-files-for-inbred). Additionally, F1 simulated results are appended to the file.
2. `output_file.windowscore.txt` --- tab-separated file
The file provides information on the strains that match to input sample for each window across the genome.

|1|2|3|4|5|6|7|8|
|---|---|---|---|---|---|---|---|
|1006|11|11|1.0|1.0|14.413599933851167|222|1|
|1158|11|11|1.0|1.0|14.413599933851167|222|1|
|1166|11|11|1.0|1.0|14.413599933851167|222|1|

Here columns are strain ID, number of SNPs matched, Informative SNPs, Probability of match, Likelihood, Likelihood ratio for the strain that matches alternatively, Number of strains that match at this window, window ID (number starting for 1 covering genome linearly).
Filtering this table by column 7 having 1 would result in homozygous windows.

3. `output_file.scores.txt.matches.json` --- JSON file

The file containing the list of matched strains, list of homozygous windows and strains matched to them and along with a simple interpretation.


## Contributing
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## History

- 1.9.2: Stable version, 24-08-2017
- 2.0.0: Stable version, 26-01-2018
- 2.1.0: Stable version, 09-08-2018


## Credits

- Rahul Pisupati (rahul.pisupati[at]gmi.oeaw.ac.at)
- Ümit Seren (uemit.seren[at]gmi.oeaw.ac.at)

## Citation

Pisupati, R. *et al.*. Verification of *Arabidopsis* stock collections using SNPmatch, a tool for genotyping high-plexed samples.  *Nature Scientific Data*  **4**, 170184 (2017).
[doi:10.1038/sdata.2017.184](https://www.nature.com/articles/sdata2017184)
