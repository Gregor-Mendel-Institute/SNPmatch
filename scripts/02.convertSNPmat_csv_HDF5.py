#!/usr/bin/python
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas
import logging
import h5py
import numpy
from pygwas.core import genotype
import scipy

def save_as_hdf5_acc(g, outHDF5):
  NumAcc = len(g.accessions)
  logging.info("Writing into HDF5 file acc wise")
  h5file = h5py.File(outHDF5, 'w')
  NumSNPs = len(g.snps)
  h5file.create_dataset('accessions', data=g.accessions, shape=(NumAcc,))
  h5file.create_dataset('positions', data=g.positions, shape=(NumSNPs,),dtype='i4')
  h5file['positions'].attrs['chrs'] = g.chrs
  h5file['positions'].attrs['chr_regions'] = g.chr_regions
  h5file.create_dataset('snps', shape=(NumSNPs, NumAcc), dtype='int8', compression="gzip", chunks=((NumSNPs, 1)))
  for i in range(NumAcc):
    h5file['snps'][:,i] = numpy.array(g.snps)[:,i]
    if i+1 % 10 == 0:
      logging.info("Written SNP info for %s accessions", i+1)
  h5file['snps'].attrs['data_format'] = g.data_format
  h5file['snps'].attrs['num_snps'] = NumSNPs
  h5file['snps'].attrs['num_accessions'] = NumAcc
  h5file.close()




inOptions = OptionParser()
inOptions.add_option("-c", "--csv_file", dest="csvFile", help="Binary CSV file as input to convert into HDF5", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")

(options, args) = inOptions.parse_args()
logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)


logging.info("Loading genotyping data")
GenotypeData = genotype.load_csv_genotype_data(options.csvFile)
logging.info("Saving binary CSV file into HDF5 file chunked rowwise")
GenotypeData.save_as_hdf5(options.hdf5File)
logging.info("Done writing the HDF5 file chunked row wise")

logging.info("Saving binary CSV file into HDF5 file chunked accession wise")
save_as_hdf5_acc(GenotypeData, options.hdf5accFile)
logging.info("Done writing!!")
