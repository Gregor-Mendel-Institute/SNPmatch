"""
    snpmatch
    ~~~~~~~~~~~~~
    The main module for running snpmatch
    :copyright: year by my name, see AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""
import os
import os.path
import argparse
import sys
from snpmatch.core import snpmatch
from snpmatch.core import csmatch
import logging, logging.config

__version__ = '1.5.0'
__updated__ = "16.11.2016"
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
  inbred_parser = subparsers.add_parser('inbred', help="snpmatch on the inbred samples")
  inbred_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  inbred_parser.add_argument("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  inbred_parser.add_argument("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  inbred_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  inbred_parser.add_argument("-o", "--output", dest="outFile", help="Output file with the probability scores")
  inbred_parser.set_defaults(func=snpmatch_inbred)
  cross_parser = subparsers.add_parser('cross', help="snpmatch on the crosses (F2s and F3s) of A. thaliana")
  cross_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  cross_parser.add_argument("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  cross_parser.add_argument("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  cross_parser.add_argument("-b", "--binLength", dest="binLen", help="Length of bins to calculate the likelihoods", default=300000)
  cross_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  cross_parser.add_argument("-c", "--crossf1", action="store_true", dest="crossf1", default=False, help="Option to identify parents in an F1 cross")
  cross_parser.add_argument("-o", "--output", dest="outFile", help="Output file with the probability scores")
  cross_parser.add_argument("-s", "--scoreFile", dest="scoreFile", help="Output of score files in each windows")
  cross_parser.set_defaults(func=snpmatch_cross)
  genocross_parser = subparsers.add_parser('genotype_cross', help="Genotype the crosses by windows given parents")
  genocross_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF file for the variants in the sample")
  genocross_parser.add_argument("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  genocross_parser.add_argument("-p", "--parents", dest="parents", help="Parents for the cross, parent1 x parent2")
  genocross_parser.add_argument("-b", "--binLength", dest="binLen", help="bin length", default=200000)
  genocross_parser.add_argument("-o", "--output", dest="outFile", help="output file")
  genocross_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  genocross_parser.set_defaults(func=genotype_cross)
  parser = subparsers.add_parser('parser', help="parse the input file")
  parser.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  parser.add_argument("-o", "--output", dest="outFile", help="output + .npz file is generater required for SNPmatch")
  parser.set_defaults(func=snpmatch_parser)
  genotyper = subparsers.add_parser('genotyper', help="SNPmatch genotyper for inbred samples")
  genotyper.add_argument("-i", "--input_file", dest="inFile", help="input file from parser -- .npz")
  genotyper.add_argument("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  genotyper.add_argument("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  genotyper.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  genotyper.add_argument("-o", "--output", dest="outFile", help="Output file with the probability scores")
  genotyper.set_defaults(func=snpmatch_genotyper)
  return inOptions

def checkARGs(args):
  if not args['hdf5File']:
    die("hdf5_file not specified")
  if not args['hdf5accFile']:
    die("hdf5accFile not specified")
  if not args['inFile']:
    die("input file not specified")
  if not os.path.isfile(args['hdf5File']):
    die("hdf5_file does not exist: " + args['hdf5File'])
  if not os.path.isfile(args['hdf5accFile']):
    die("hdf5accFile does not exist: " + args['hdf5accFile'])
  if not os.path.isfile(args['inFile']):
    die("input file does not exist: " + args['inFile'])

def snpmatch_inbred(args):
  checkARGs(args)
  _,inType = os.path.splitext(args['inFile'])
  if inType == '.vcf':
    snpmatch.match_vcf_to_acc(args)
  elif inType == '.bed':
    snpmatch.match_bed_to_acc(args)
  else:
    die("input file type %s not supported", inType)

def snpmatch_cross(args):
  checkARGs(args)
  _,inType = os.path.splitext(args['inFile'])
  if args['crossf1']:
    csmatch.crossF1genotyper(args)
  else:
    if not args['outFile']:
      die("specify an output file")
    if not args['scoreFile']:
      die("file to give out scores is not specified")
    if inType == '.vcf':
      csmatch.match_vcf_to_acc(args)
    elif inType == '.bed':
      csmatch.match_bed_to_acc(args)
    else:
      die("input file type %s not supported", inType)

def snpmatch_parser(args):
  if not args['inFile']:
    die("input file not specified")
  if not args['outFile']:
    args['outFile'] = args['inFile']
  snpmatch.parseInput(args)

def snpmatch_genotyper(args):
  checkARGs(args)
  _,inType = os.path.splitext(args['inFile'])
  if inType == ".npz":
    snpmatch.potatoGenotyper(args)
  else:
    die("input for the genotyper should be .npz file generated from SNPmatch parser")

def genotype_cross(args):
  #checkARGs(args)
  if not args['parents']:
    die("parents not specified")
  csmatch.crossGenotyper(args)

def main():
  ''' Command line options '''
  program_version = "v%s" % __version__
  program_build_date = str(__updated__)
  program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
  program_shortdesc = "The main module for snpmatch"
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
    return 0
  try:
    args['func'](args)
    return 0
  except KeyboardInterrupt:
    return 0
  except Exception as e:
    logging.exception(e)
    return 2

if __name__=='__main__':
  sys.exit(main())

