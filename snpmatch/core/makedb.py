import logging
import h5py
import numpy as np
from pygwas.core import genotype
import sys
import os
import os.path
import json
import gzip
from subprocess import Popen, PIPE
import shlex

log = logging.getLogger(__name__)
def die(msg):
    sys.stderr.write('Error: ' + msg + '\n')
    sys.exit(1)

def getSamples(inVCF):
    vname, vext = os.path.splitext(inVCF)
    if vext == '.vcf':
        v = open(inVCF, 'r')
    elif vext == '.gz':
        v = gzip.open(inVCF)
    else:
        die('please provide VCF file either gunzipped or text!')
    eline = v.readline()
    while eline[0] == '#':
        if eline[0:6] == '#CHROM':
            info = eline.rstrip().split('\t')
        eline = v.readline()
    samples = np.array(info[9:])
    return samples

def getCSV(inVCF, outFile, bcfpath = ''):
    outcsv = open(outFile, "w")
    samples = getSamples(inVCF)
    outcsv.write('Chromosome,Position')
    for es in samples:
        outcsv.write(',%s' % es)
    outcsv.write('\n')
    outcsv.truncate()
    log.info('running bcftools!')
    if bcfpath == '':
        bcftool_command = "bcftools query -f \"%CHROM,%POS[,%GT]\n\" " + inVCF
    else:
        bcftool_command = bcfpath + '/' + "bcftools query -f \"%CHROM,%POS[,%GT]\n\" " + inVCF
    sed_command = "| sed 's/,0[\/|]0/,0/g' | sed 's/,1[\/|]1/,1/g' | sed 's/,0[\/|]1/,2/g' | sed 's/,\.[\/|]\./,-1/g'"
    full_command = bcftool_command + sed_command
    convertcsv = Popen(full_command, shell=True, stdout = outcsv)
    convertcsv.wait()
    outcsv.close()
    log.info('done!')

def save_as_hdf5_acc(g, outHDF5):
    NumAcc = len(g.accessions)
    log.info("Writing into HDF5 file acc wise")
    h5file = h5py.File(outHDF5, 'w')
    NumSNPs = len(g.snps)
    h5file.create_dataset('accessions', data=g.accessions, shape=(NumAcc,))
    h5file.create_dataset('positions', data=g.positions, shape=(NumSNPs,),dtype='i4')
    h5file['positions'].attrs['chrs'] = g.chrs
    h5file['positions'].attrs['chr_regions'] = g.chr_regions
    h5file.create_dataset('snps', shape=(NumSNPs, NumAcc), dtype='int8', compression="gzip", chunks=((NumSNPs, 1)))
    for i in range(NumAcc):
        h5file['snps'][:,i] = np.array(g.snps)[:,i]
        if i+1 % 10 == 0:
            log.info("written SNP info for %s accessions", i+1)
    h5file['snps'].attrs['data_format'] = g.data_format
    h5file['snps'].attrs['num_snps'] = NumSNPs
    h5file['snps'].attrs['num_accessions'] = NumAcc
    h5file.close()

def makeHDF5s(csvFile, outFile):
    GenotypeData = genotype.load_csv_genotype_data(csvFile)
    log.info("saving CSV file into HDF5 file chunked rowwise")
    GenotypeData.save_as_hdf5(outFile + '.hdf5')
    log.info("done!")
    log.info("saving CSV file into HDF5 file chunked accession wise")
    save_as_hdf5_acc(GenotypeData, outFile + '.acc.hdf5')
    logging.info("done!")

def makedb_from_vcf(args):
    _,inType = os.path.splitext(args['inFile'])
    if inType == '.vcf':
        log.info("converting VCF to CSV")
        getCSV(args['inFile'], args['db_id'] + '.csv', args['bcfpath'])
        log.info("converting CSV to hdf5!")
        makeHDF5s(args['db_id'] + '.csv', args['db_id'])
        log.info('done!')
    elif inType == '.csv':
        log.info("converting CSV to hdf5!")
        makeHDF5s(args['inFile'], args['db_id'])
        log.info('done!')
    else:
        die("please provide either a VCF file or a CSV!")
