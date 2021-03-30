import logging
import h5py
import numpy as np
from snpmatch.pygwas import genotype
import sys
import os
import os.path
import json
import re
from subprocess import Popen, PIPE, check_output

log = logging.getLogger(__name__)
def die(msg):
    sys.stderr.write('Error: ' + msg + '\n')
    sys.exit(1)

def run_command(command, stdout=""):
    if stdout == "":
        out_cmd = check_output( command, shell=True )
        return(out_cmd.rstrip().split("\n"))
    out_cmd = Popen(command, shell=True, stdout = stdout)
    out_cmd.wait()

def get_contigs(vcf_header):
    chr_names = []
    chr_len = []
    for eh in vcf_header:
        if eh[0:8] == '##contig':
            chr_names.append( eh.replace(">", "").replace("<","").split("ID=")[1].split(",")[0] )
            chr_len.append( int(eh.replace(">", "").replace("<","").split("length=")[1].split(",")[0]) )
    return({"ref_chrs": chr_names, "ref_chrlen": chr_len})

def getCSV(inVCF, outFile, bcfpath = ''):
    try:
        full_bcfpath = ""
        if bcfpath  == '':
            check_output('bcftools -h', shell=True)
            full_bcfpath = "bcftools "
        else:
            check_output(bcfpath + '/bcftools -h', shell=True)
            full_bcfpath = bcfpath + '/bcftools '
    except:
        die("please provide bcftools installation path. '%s' is not found!" % (bcfpath + '/bcftools'))
    outcsv = open(outFile + '.csv', "w")
    vcf_header = run_command( full_bcfpath + " view -h " + inVCF )
    samples = np.array(run_command( full_bcfpath + " query -l " + inVCF ))
    contigs = get_contigs(vcf_header)
    log.info("Number of contigs found: %s" % len(contigs['ref_chrs']))
    with open(outFile + ".json", "w") as out_stats:
        out_stats.write(json.dumps(contigs, sort_keys=True, indent=4))
    outcsv.write('Chromosome,Position')
    for es in samples:
        outcsv.write(',%s' % es)
    outcsv.write('\n')
    outcsv.truncate()
    log.info('running bcftools!')
    bcftool_command = full_bcfpath + " query -f \"%CHROM,%POS[,%GT]\n\" " + inVCF
    sed_command = "| sed 's/,0[\/|]0/,0/g' | sed 's/,1[\/|]1/,1/g' | sed 's/,0[\/|]1/,2/g'| sed 's/,1[\/|]0/,2/g' | sed 's/,\.[\/|]\./,-1/g'"
    run_command( bcftool_command + sed_command, outcsv)
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
    if inType == '.vcf' or len(re.compile(".vcf.gz$").findall(os.path.basename(args['inFile']))) > 0 :
        log.info("converting VCF to CSV")
        getCSV(args['inFile'], args['db_id'], args['bcfpath'])
        log.info("converting CSV to hdf5!")
        makeHDF5s(args['db_id'] + '.csv', args['db_id'])
        log.info('done!')
    elif inType == '.csv':
        log.info("converting CSV to hdf5!")
        makeHDF5s(args['inFile'], args['db_id'])
        log.info('done!')
    else:
        die("please provide either a VCF file or a CSV!")
