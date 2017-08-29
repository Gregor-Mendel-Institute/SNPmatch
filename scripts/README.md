# Converting VCF file to HDF5 database files

The database files are required for SNPmatch analysis. These can be generated from the joint VCF files which contain genotype information for multiple strains.
The scripts given in this folder can be used to convert a joint VCF file into HDF5 database file. These files are later used in SNPmatch analysis.

## Dependencies for the scripts

These are tested and work only on Linux based systems. The below dependencies should be present in either bash environment or python path.

1) BCFtools
3) h5py
4) numpy
5) scipy
6) PyGWAS


## Input files and Usage

The input file is a VCF file either zipped or plain text file. We use zgrep and bcftools and other bash commands to parse the VCF file. The bash script in the folder can be run as given below.

```bash
bash ./01.convertVCF_SNPmat_csv.sh input_file db
```

The output to the above command are three database files.

1) db.csv
2) db.hdf5
3) db.acc.hdf5

The two hdf5 files present are the main database files which are used in the SNPmatch analysis. The files have the same information but are chunked for better efficiency.
These files db.hdf5 and db.acc.hdf5 files are given to the SNPmatch command under -d and -e options respectively.
