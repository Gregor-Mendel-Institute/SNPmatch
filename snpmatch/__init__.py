"""
    snpmatch
    ~~~~~~~~~~~~~
    The main module for running snpmatch
    :copyright: year by my name, see AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""
import os
import argparse
import sys
from snpmatch.core import snpmatch
from snpmatch.core import csmatch

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)


def get_options():
  inOptions = argparse.ArgumentParser()
  subparsers = inOptions.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')

  inbred_parser = subparsers.add_parser('inbred', help="snpmatch on the inbred samples")
  inbred_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  inbred_parser.add_argument("-t", "--input_type", dest="inType", help="Type of the input file given. Possible inputs: 'vcf', 'bed'")
  inbred_parser.add_argument("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  inbred_parser.add_argument("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  inbred_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  inbred_parser.add_argument("-o", "--output", dest="outFile", help="Output file with the probability scores")
  inbred_parser.set_defaults(func=snpmatch_inbred)
  cross_parser = subparsers.add_parser('cross', help="snpmatch on the crosses (F2s and F3s) of A. thaliana")
  cross_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  cross_parser.add_argument("-t", "--input_type", dest="inType", help="Type of the input file given. Possible inputs: 'vcf', 'bed'")
  cross_parser.add_argument("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  cross_parser.add_argument("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  cross_parser.add_argument("-b", "--binLength", dest="binLen", help="Length of bins to calculate the likelihoods", default=300000)
  cross_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
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
  if not args['inType']:
    die("not mentioned the type of input")
  if args['inType'] == "vcf":
    snpmatch.match_vcf_to_acc(args)
  elif args['inType'] == "bed":
    snpmatch.match_bed_to_acc(args)

def snpmatch_cross(args):
  checkARGs(args)
  if not args['inType']:
    die("not mentioned the type of input")
  if not args['output']:
    die("specify an output file")
  if not args['scoreFile']:
    die("file to give out scores is not specified")
  if args['inType'] == "vcf":
    csmatch.match_vcf_to_acc(args)
  elif args['inType'] == "bed":
    csmatch.match_bed_to_acc(args)

def genotype_cross(args):
  #checkARGs(args)
  if not args['parents']:
    die("parents not specified")
  csmatch.genotypeCross(args)

def main():
  parser = get_options()
  args = vars(parser.parse_args())
  if 'func' not in args:
    parser.print_help()
    return 0
  try:
    args['func'](args)
    return 0
  except KeyboardInterrupt:
    return 0
  except Exception as e:
    return 2


if __name__=='__main__':
  sys.exit(main())

