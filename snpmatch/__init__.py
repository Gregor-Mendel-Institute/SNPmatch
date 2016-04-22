"""
    snpmatch
    ~~~~~~~~~~~~~
    The main module for running SNPmatch
    :copyright: year by my name, see AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""

import argparse
import sys
from snpmatch.core import snpmatch
from snpmatch.core import csmatch


def get_options():
  inOptions = argparse.ArgumentParser()
  subparsers = inOptions.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')

  inbred_parser = subparsers.add_parser('inbred', help="Run the SNPmatch on the inbred samples")
  inbred_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  inbred_parser.add_argument("-t", "--input_type", dest="inType", help="Type of the input file given. Possible inputs: 'vcf', 'bed'")
  inbred_parser.add_argument("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  inbred_parser.add_argument("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  inbred_parser.add_argument("-o", "--output", dest="outFile", help="Output file with the probability scores")
  inbred_parser.set_defaults(func=snpmatch_inbred)
  
  cross_parser = subparsers.add_parser('cross', help="Run the SNPmatch on the crosses (F2s and F3s), works on A. thaliana due to Chr length")
  cross_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  cross_parser.add_argument("-t", "--input_type", dest="inType", help="Type of the input file given. Possible inputs: 'vcf', 'bed'")
  cross_parser.add_argument("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  cross_parser.add_argument("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  cross_parser.add_argument("-o", "--output", dest="outFile", help="Output file with the probability scores")
  cross_parser.add_argument("-s", "--scoreFile", dest="scoreFile", help="Output of score files in each windows")
  cross_parser.add_argument("-b", "--binLength", dest="binLen", help="Length of bins to calculate the likelihoods", default=300000)
  cross_parser.set_defaults(func=snpmatch_cross)
  
  return inOptions

def snpmatch_inbred(args):
  if args['inType'] == "vcf":
    snpmatch.match_vcf_to_acc(args)
  elif args['inType'] == "bed":
    snpmatch.match_bed_to_acc(args)

def snpmatch_cross(args):
  if args['inType'] == "vcf":
    csmatch.match_vcf_to_acc(args)
  elif args['inType'] == "bed":
    csmatch.match_bed_to_acc(args)

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

