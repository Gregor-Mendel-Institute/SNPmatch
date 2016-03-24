#!/usr/bin/python
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas
# module load vcfnp
import logging
import numpy
import numpy.ma
from pygwas.core import genotype
import vcfnp
import pandas
import scipy
import math

def nCr(n, r):
  f = math.factorial
  return f(n) / (f(r) * f(n-r))

def likeliTest(n, y):
  p = 0.999999  ### Since computing the right likelihood is troubling
  pS = float(y)/n
  a = y * scipy.log(pS/p)
  b = (n - y) * scipy.log((1-pS)/(1-p))
#  c = scipy.log(nCr(n, y))
  return(a+b)

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-i", "--input_file", dest="inFile", help="VCF/BED file for the variants in the sample", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file chunked row-wise", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file chunked column-wise", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-t", "--input_type", dest="input_type", help="Type of the Input given. Possible inputs: 'vcf', 'bed'", type="string")

(options, args) = inOptions.parse_args()

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
num_lines = len(GenotypeData.accessions)

if options.input_type == "vcf":
  logging.info("Reading the VCF file")
  vcf = vcfnp.variants(options.inFile, cache=True).view(numpy.recarray)
  vcfD = vcfnp.calldata_2d(options.inFile, cache=True).view(numpy.recarray)
  DPthres = numpy.mean(vcf.DP[numpy.where(vcf.DP > 0)[0]]) * 4
  snpsREQ = numpy.where((vcfD.is_called[:,0]) & (vcf.QUAL > 30) & (vcf.DP > 0) & (vcf.DP < DPthres))[0]
  snpCHR = numpy.array(numpy.char.replace(vcf.CHROM[snpsREQ], "Chr", "")).astype("int8")
  snpPOS = numpy.array(vcf.POS[snpsREQ])
  snpDP = vcf.DP[snpsREQ]
  snpGT = vcfD.GT[snpsREQ]
#  snpPL = vcfD.PL[snpsREQ, 0]   
#  snpWEI = numpy.copy(snpPL)
#  snpWEI = snpWEI.astype(float)
#  snpWEI = snpWEI/(-10)
#  snpWEI = numpy.exp(snpWEI)
elif options.input_type == "pos":
  logging.info("Reading the BED file")
  targetSNPs = pandas.read_table(options.inFile, header=None, usecols=[0,1,2])
  snpCHR = numpy.array(targetSNPs[0], dtype=int) 
  snpPOS = numpy.array(targetSNPs[1], dtype=int)
  snpGT = numpy.array(targetSNPs[2])
  
  
ScoreList = numpy.zeros(num_lines, dtype="float")
NumInfoSites = numpy.zeros(len(GenotypeData.accessions), dtype="uint32")
NumMatSNPs = 0
chunk_size = 1000

for i in range(1,6):
  perchrTarPos = numpy.where(snpCHR == i)[0]
  perchrtarSNPpos = snpPOS[perchrTarPos]
  logging.info("Loaded %s chromosome positions from the position file", i)
  start = GenotypeData.chr_regions[i-1][0]
  end = GenotypeData.chr_regions[i-1][1]
  chrpositions = GenotypeData.positions[start:end]
  matchedAccInd = numpy.where(numpy.in1d(chrpositions, perchrtarSNPpos))[0] + start
  matchedTarInd = numpy.where(numpy.in1d(perchrtarSNPpos, chrpositions))[0]
  matchedTarGTs = snpGT[perchrTarPos[matchedTarInd],]
  TarGTs = numpy.zeros(len(matchedTarGTs), dtype="int8")
  TarGTs[numpy.where(matchedTarGTs == "1/1")[0]] = 1
  TarGTs1 = numpy.ones(len(matchedTarGTs), dtype="int8")
  TarGTs1[numpy.where(matchedTarGTs == "0/0")[0]] = 0  
  NumMatSNPs = NumMatSNPs + len(matchedAccInd)
  for j in range(0, len(matchedAccInd), chunk_size):
    t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
    samSNPs = numpy.reshape(numpy.repeat(TarGTs[j:j+chunk_size], num_lines), (len(TarGTs[j:j+chunk_size]),num_lines))
    samSNPs1 = numpy.reshape(numpy.repeat(TarGTs1[j:j+chunk_size], num_lines), (len(TarGTs1[j:j+chunk_size]),num_lines))
    ScoreList = ScoreList + (numpy.sum(t1001SNPs == samSNPs, axis=0) + numpy.sum(t1001SNPs == samSNPs1, axis=0))/float(2)
    if(len(TarGTs[j:j+chunk_size]) > 1):
      NumInfoSites = NumInfoSites + len(TarGTs[j:j+chunk_size]) - numpy.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
    elif(len(TarGTs[j:j+chunk_size]) == 1): 
      NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
  logging.info("Done analysing %s positions", NumMatSNPs)

logging.info("Done calculating the scores for each accession")

LikeLiHoods = [likeliTest(NumInfoSites[i], int(ScoreList[i])) for i in range(num_lines)]
LikeLiHoods = numpy.array(LikeLiHoods).astype("float")

TopHit = numpy.amin(LikeLiHoods)
LikeLiHoodRatio = [LikeLiHoods[i]/TopHit for i in range(num_lines)]
LikeLiHoodRatio = numpy.array(LikeLiHoodRatio).astype("float")
TopHitAcc = numpy.where(LikeLiHoodRatio < 3.841)[0]
logging.info("Number of ambiguous accessions: %s", len(TopHitAcc))

outfile = open(options.outFile, 'w')
for i in range(num_lines):
  score = float(ScoreList[i])/NumInfoSites[i]
  outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatio[i], NumMatSNPs, len(snpPOS)))
outfile.close()



