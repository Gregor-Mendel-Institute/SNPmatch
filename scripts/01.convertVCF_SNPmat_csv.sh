
#######
# Input File: VCF file
inFile=$1
outFile=$2

script_path=`dirname $0`

# This script needs BCFtools to be install in the system

echo "Starting to convert a gVCF to binary CSV file"
zgrep -m 1 "^#CHROM" ${inFile} | cut -f 1,2,10- | sed 's/\t/,/g' | sed 's/^#CHROM,POS/Chromosome,Position/' > ${outFile}.csv
bcftools query -f "%CHROM,%POS[,%GT]\n" ${inFile} | sed 's/0\/0/0/g' | sed 's/1\/1/1/g' | sed 's/0\/1/2/g' | sed 's/\.\/\./-1/g' >> ${outFile}.csv

## For the below command we need python installed

echo "Converting CSV file to an HDF5 files"
python $script_path/02.convertSNPmat_csv_HDF5.py -c ${outFile}.csv -d ${outFile}.hdf5 -e ${outFile}e.acc.hdf5
