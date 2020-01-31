"""
    SNPmatch
    ~~~~~~~~~~~~~
    The main module for running SNPmatch
    :copyright: year by my name, see AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""
import os
import os.path
import argparse
import sys
from snpmatch.core import snpmatch
from snpmatch.core import csmatch
from snpmatch.core import genotype_cross as gtm
from snpmatch.core import makedb
from snpmatch.core import simulate
import logging, logging.config

__version__ = '4.0.0'
__updated__ = "31.01.2020"
__date__ = "25.10.2016"

def setLog(logDebug):
  log = logging.getLogger()
  if logDebug:
    numeric_level = getattr(logging, "DEBUG", None)
  else:
    numeric_level = getattr(logging, "ERROR", None)
  log_format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
  lch = logging.StreamHandler()
  lch.setLevel(numeric_level)
  lch.setFormatter(log_format)
  log.setLevel(numeric_level)
  log.addHandler(lch)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def get_options(program_license,program_version_message):
  inOptions = argparse.ArgumentParser(description=program_license)
  inOptions.add_argument('-V', '--version', action='version', version=program_version_message)
  subparsers = inOptions.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')
  inbred_parser = subparsers.add_parser('inbred', help="SNPmatch on the inbred samples")
  inbred_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  inbred_parser.add_argument("-d", "--hdf5_file", default = None, dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  inbred_parser.add_argument("-e", "--hdf5_acc_file", default = None, dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  inbred_parser.add_argument("--refine", action="store_true", dest="refine", default=False, help="Refine scores for indistinguishable lines")
  inbred_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  inbred_parser.add_argument("-o", "--output", dest="outFile", default="identify_inbred", help="Output file with the probability scores")
  inbred_parser.set_defaults(func=snpmatch_inbred)

  cross_parser = subparsers.add_parser('cross', help="SNPmatch on the crosses (F2s and F3s) of A. thaliana")
  cross_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  cross_parser.add_argument("-d", "--hdf5_file", default = None, dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  cross_parser.add_argument("-e", "--hdf5_acc_file",  default = None, dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  cross_parser.add_argument("-b", "--binLength", dest="binLen", help="Length of bins to calculate the likelihoods", default=300000, type=int)
  cross_parser.add_argument("--genome", dest="genome", default="athaliana_tair10", help="Path to Reference JSON file, if you are working with non-thaliana tair10 assembly")
  cross_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  cross_parser.add_argument("-o", "--output", dest="outFile", default="identify_cross", help="Output files with the probability scores and scores along windows")
  cross_parser.set_defaults(func=snpmatch_cross)

  genocross_parser = subparsers.add_parser('genotype_cross', help="Genotype the crosses by windows given parents")
  genocross_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF file for the variants in the sample")
  genocross_parser.add_argument("-a", "--all_samples", action="store_true", dest="all_samples", default=False, help="Provide this argument if input file is a matrix generated from GVCF containing all samples, needed to be genotyped.")
  genocross_parser.add_argument("-d", "--hdf5_file", default = None, dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  genocross_parser.add_argument("-e", "--hdf5_acc_file", default = None, dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  genocross_parser.add_argument("-p", "--parents", dest="parents", help="Parents for the cross, parent1 x parent2")
  genocross_parser.add_argument("-q", "--father", dest="father", help="If given this should be VCF file for the mother (ex., if the cross is parent1 x parent2, parent2.vcf should be the input file. Also -p should be parent1.vcf")
  genocross_parser.add_argument("-b", "--binLength", dest="binLen", help="bin length",type=int, default=200000)
  genocross_parser.add_argument("--good_samples", dest="good_samples", help="list of sample IDs which need to be present, only works with '-a' option", default=None)
  genocross_parser.add_argument("--lr_thres", dest="lr_thres", default=1.5, help="Likelihood ratio threshold for genotype calling.")
  genocross_parser.add_argument("--hmm", dest="hmm", action="store_true", help="Use HMM Viterbi method to determine underlying genotypes.")
  genocross_parser.add_argument("--genome", dest="genome", default="athaliana_tair10", help="Path to Reference JSON file, if you are working with non-thaliana tair10 assembly")
  genocross_parser.add_argument("-o", "--output", dest="outFile", default="genotype_cross", help="output file")
  genocross_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  genocross_parser.set_defaults(func=genotype_cross)

  parser = subparsers.add_parser('parser', help="parse the input file")
  parser.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  parser.add_argument("-o", "--output", dest="outFile", help="output + .npz file is generater required for SNPmatch")
  parser.set_defaults(func=snpmatch_parser)

  pairparser = subparsers.add_parser('pairsnp', help="pairwise comparison of two snp files")
  pairparser.add_argument("-i", "--input_file_1", dest="inFile_1", help="VCF/BED file for the variants in the sample one")
  pairparser.add_argument("-j", "--input_file_2", dest="inFile_2", help="VCF/BED file for the variants in the sample two")
  pairparser.add_argument("-d", "--hdf5_file", dest="hdf5File", default=None, help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  pairparser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  pairparser.add_argument("-o", "--output", dest="outFile", default="pairsnp", help="output json file")
  pairparser.set_defaults(func=snpmatch_paircomparions)

  makedbparser = subparsers.add_parser('makedb', help="Create database files from given VCF, only give biallelic SNPs")
  makedbparser.add_argument("-i", "--input_vcf", dest="inFile", help="input VCF file for the known strains. You can also provide a CSV file which is an intermediate file in the process.")
  makedbparser.add_argument("-p", "--bcftools_path", dest="bcfpath", help="path to the bcftools executable. Not necessary if present in BASH PATH", default='')
  makedbparser.add_argument("-o", "--out_db_id", dest="db_id", help="output id for database files")
  makedbparser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  makedbparser.set_defaults(func=makedb_vcf_to_hdf5)

  simparser = subparsers.add_parser('simulate', help="Given SNP database, check the genotyping efficiency randomly selecting 'n' number of SNPs")
  simparser.add_argument("-d", "--hdf5_file",  default = None, dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  simparser.add_argument("-e", "--hdf5_acc_file",  default = None, dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  simparser.add_argument("-a", "--ecotype_id", dest="AccID", help= "Ecotype ID you want draw the SNPs")
  simparser.add_argument("-n", "--number_of_snps", dest="numSNPs", help= "number of SNPs to draw in random to genotype the sample", type=int)
  simparser.add_argument("-p", "--error_rate", dest="err_rate", help= "error rate while matching the SNPs, error rate of 0 gives perfect match to the accession", default = 0.001, type=float)
  simparser.add_argument("--f1", action="store_true", dest="simF1", default=False, help="Simulate SNPs for an F1, give parents as 1061x1062 in argument '-a'")
  simparser.add_argument("-o", "--output", dest="outFile", help="Output file with scores")
  simparser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  simparser.set_defaults(func=simulate_snps)
  return inOptions

def check_file(inFile):
  if not inFile:
    die("file: %s not specified" % inFile)
  if not os.path.isfile(inFile):
    die("input file does not exist: " + inFile)

def snpmatch_inbred(args):
  check_file(args['inFile'])
  snpmatch.potatoGenotyper(args)

def snpmatch_cross(args):
  check_file(args['inFile'])
  csmatch.potatoCrossIdentifier(args)

def snpmatch_parser(args):
  check_file(args['inFile'])
  if not args['outFile']:
    if os.path.isfile(args['inFile'] + ".snpmatch.npz"):
      os.remove(args['inFile'] + ".snpmatch.npz")
  from snpmatch.core import parsers
  parsers.potatoParser(inFile = args['inFile'], logDebug =  args['logDebug'], outFile = args['outFile'])

def genotype_cross(args):
    #checkARGs(args)
    if not args['parents']:
        die("parents not specified")
    gtm.potatoCrossGenotyper(args)

def snpmatch_paircomparions(args):
    check_file(args['inFile_1'])
    check_file(args['inFile_2'])
    snpmatch.pairwiseScore(args['inFile_1'], args['inFile_2'], args['logDebug'], args['outFile'], args['hdf5File'])

def makedb_vcf_to_hdf5(args):
    check_file(args['inFile'])
    makedb.makedb_from_vcf(args)

def simulate_snps(args):
    simulate.potatoSimulate(args)


def main():
  ''' Command line options '''
  program_version = "v%s" % __version__
  program_build_date = str(__updated__)
  program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
  program_shortdesc = "The main module for SNPmatch"
  program_license = '''%s
  Created by Rahul Pisupati on %s.
  Copyright 2016 Gregor Mendel Institute. All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
USAGE
''' % (program_shortdesc, str(__date__))

  parser = get_options(program_license,program_version_message)
  args = vars(parser.parse_args())
  setLog(args['logDebug'])
  if 'func' not in args:
    parser.print_help()
    return(0)
  try:
    args['func'](args)
    return(0)
  except KeyboardInterrupt:
    return(0)
  except Exception as e:
    logging.exception(e)
    return(2)

if __name__=='__main__':
  sys.exit(main())
