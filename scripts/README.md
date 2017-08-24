

# README for converting VCF file to HDF5 database file

The scripts given in this folder can be used to convert a joint VCF file into a HDF5 database file. These files are later used for genotyping oterh samples using SNPmatch.

# Dependencies for these scripts

These are tested and work on Linux based systems only.

1) BCFtools
2) python
3) h5py
4) numpy
5) scipy
6) pygwas


# Input files and usage

The input file is a VCF file which can be read through bcftools. The script generates couple of files as output.

```bash
bash ./01.convertVCF_SNPmat_csv.sh input_file output_file_id
```

Given the above command it generates files output_file_id.csv, output_file_id.hdf5 and output_file_id.acc.hdf5

The two hdf5 files have the same information but chunked in a different way for efficient readability.
