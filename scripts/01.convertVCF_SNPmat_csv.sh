#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N convertVCF
#PBS -V
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -o log.o.convertVCF
#PBS -e log.e.convertVCF


cd $PBS_O_WORKDIR
inFile=`ls 1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.vcf | sed 's/\.vcf$//' `

# Convert VCF file into hetfiltered CSV file

echo "Starting to convert a gVCF to binary CSV file"
grep -v "^##" ${inFile}.vcf | cut -f 1,2,10- | sed 's/^Chr//' | awk ' {for (i = 3; i <= NF; i++)
if(substr($i,0,3) == "0/0")
$i = 0
else if (substr($i, 0,3) == "1/1")
$i = 1
else if(substr($i, 0,3) == "1/0" || substr($i, 0,3) == "0/1")
$i = -1
else if(substr($i,0,3) == "./.")
$i = -1
print $0}' OFS=',' |sed 's/\t/,/g' | sed 's/^#CHROM,POS/Chromosome,Position/' > $inFile.hetfiltered.csv


## For the compressed vcf file

(zgrep -m 1 "^#CHROM" 1135g_SNP_BIALLELIC.vcf.gz | cut -f 1,2,10- && zgrep -v "^#" 1135g_SNP_BIALLELIC.vcf.gz | cut -f 1,2,10- | sed 's/^Chr//' | awk ' {for (i = 3; i <= NF; i++)
if(substr($i,1,3) ~ /0.0/)
$i = 0
else if (substr($i,1,3) ~ /1.1/)
$i = 1
else if(substr($i,1,3) ~ /1.0/ || substr($i, 1,3) ~ /0.1/)
$i = 2
else if(substr($i,1,3) ~ /\..\./)
$i = -1
print $0}' OFS=',' ) | sed 's/\t/,/g' | sed 's/^#CHROM,POS/Chromosome,Position/' > 1135g_SNP_BIALLELIC.SNPmatrix_11-Dec-2015.csv


#grep -v "^##" ${inFile}.vcf | cut -f 1,2,10- | sed 's/^Chr//' | awk '{for (i = 3; i <= NF; i++)
#if(substr($i,1,3) == "0/0")
#$i = 0
#else if (substr($i, 1,3) == "1/1")
#$i = 1
#else if(substr($i, 1,3) == "1/0" || substr($i, 1,3) == "0/1")
#$i = -1
#else if(substr($i,1,3) == "./.")
#$i = -2
#print $0}' OFS=',' |sed 's/\t/,/g' | sed 's/^#CHROM,POS/Chromosome,Position/' > ${inFile}.csv
## Done with awk

echo "Saving as HDF5 files"

module load h5py
module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas

python ./scripts/02.convertSNPmat_csv_HDF5.py -c $inFile.hetfiltered.csv -d $inFile.hetfiltered.hdf5 -e $inFile.hetfiltered.acc.hdf5
