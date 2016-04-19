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


def get_options():
  inOptions = argparse.ArgumentParser()
  inOptions.add_argument("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample")
  inOptions.add_argument("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise")
  inOptions.add_argument("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise")
  inOptions.add_argument("-o", "--output", dest="outFile", help="Output file with the probability scores")
  inOptions.add_argument("-t", "--input_type", dest="inType", help="Type of the Input given. Possible inputs: 'vcf', 'bed'")
  return inOptions

def main():
  parser = get_options()
  args = vars(parser.parse_args())
  if args['inType'] is None:
    parser.print_help()
    return 0
  try:
    snpmatch(args)
    return 0
  except KeyboardInterrupt:
    ### handle keyboard interrupt ###
    return 0
 

def snpmatch(args):
  if args['inType'] is "vcf":
    snpmatch.match_vcf_to_acc(args)
  elif args['inType'] is "bed":
    snpmatch.match_bed_to_acc(args)
  else:
    parser.print_help()


if __name__=='__main__':
  sys.exit(main())

