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

# Convert VCF file into hetfiltered CSV file

module load BCFtools/1.3-foss-2015b

echo "Starting to convert a gVCF to binary CSV file"

inFile=`ls *vcf | head -n 1 | sed 's/\.vcf$//'`
grep -m 1 "^#CHROM" ${inFile}.vcf | cut -f 1,2,10- | sed 's/\t/,/g' | sed 's/^#CHROM,POS/Chromosome,Position/' > ${inFile}.csv
#bcftools query -f "%CHROM\t%POS[\t%GT]\n" ${inFile}.vcf | sed 's/\b0\/0/0/g' | sed 's/\b1\/1\b/1/g' | sed 's/\b0\/1\b/2/g' | sed 's/\.\/\./-1/g' | sed 's/\t/,/g'  >> ${inFile}.csv
bcftools query -f "%CHROM,%POS[,%GT]\n" ${inFile}.vcf | sed 's/0\/0/0/g' | sed 's/1\/1/1/g' | sed 's/0\/1/2/g' | sed 's/\.\/\./-1/g' >> ${inFile}.csv


## For the compressed vcf file
#inFile=`ls *vcf.gz | head -n 1 | sed 's/\.vcf\.gz$//'`
#zgrep -m 1 "^#CHROM" ${inFile}.vcf.gz | cut -f 1,2,10- | sed 's/\t/,/g' | sed 's/^#CHROM,POS/Chromosome,Position/' > ${inFile}.csv
#bcftools query -f "%CHROM\t%POS[\t%GT]\n" ${inFile}.vcf.gz | sed 's/\b0\/0/0/g' | sed 's/\b1\/1\b/1/g' | sed 's/\b0\/1\b/2/g' | sed 's/\.\/\./-1/g' | sed 's/\t/,/g'  >> ${inFile}.csv

#(zgrep -m 1 "^#CHROM" 1135g_SNP_BIALLELIC.vcf.gz | cut -f 1,2,10- && zgrep -v "^#" 1135g_SNP_BIALLELIC.vcf.gz | cut -f 1,2,10- | sed 's/^Chr//' | awk ' {for (i = 3; i <= NF; i++)
#if(substr($i,1,3) ~ /0.0/)
#$i = 0
#else if (substr($i,1,3) ~ /1.1/)
#$i = 1
#else if(substr($i,1,3) ~ /1.0/ || substr($i, 1,3) ~ /0.1/)
#$i = 2
#else if(substr($i,1,3) ~ /\..\./)
#$i = -1
#print $0}' OFS=',' ) | sed 's/\t/,/g' | sed 's/^#CHROM,POS/Chromosome,Position/' > 1135g_SNP_BIALLELIC.SNPmatrix_11-Dec-2015.csv

echo "Saving as HDF5 files"

module load h5py
module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas

python ./scripts/02.convertSNPmat_csv_HDF5.py -c $inFile.csv -d $inFile.hdf5 -e $inFile.acc.hdf5
