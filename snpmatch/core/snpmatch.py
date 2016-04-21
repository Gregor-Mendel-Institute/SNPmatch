"""
  SNPmatch
"""
import numpy
import numpy.ma
from pygwas.core import genotype
import vcfnp
import pandas
import argparse
import logging
import sys
import os.path

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def likeliTest(n, y):
  p = 0.999999
  if n > 0:
    pS = float(y)/n
    a = y * numpy.log(pS/p)
    b = (n - y) * numpy.log((1-pS)/(1-p))
    return(a+b)
  elif n == y:
    return 1
  else:
    return numpy.nan

def calculate_likelihoods(ScoreList, NumInfoSites):
  num_lines = len(ScoreList)
  LikeLiHoods = [likeliTest(NumInfoSites[i], int(ScoreList[i])) for i in range(num_lines)]
  LikeLiHoods = numpy.array(LikeLiHoods).astype("float")
  TopHit = numpy.amin(LikeLiHoods)
  LikeLiHoodRatios = [LikeLiHoods[i]/TopHit for i in range(num_lines)]
  LikeLiHoodRatios = numpy.array(LikeLiHoodRatios).astype("float")
  return (LikeLiHoods, LikeLiHoodRatios)

def print_out_table(outFile, GenotypeData, ScoreList, NumInfoSites, LikeLiHoods, LikeLiHoodRatios, NumMatSNPs):
  out = open(outFile, 'w')
  for i in range(len(GenotypeData.accessions)):
    score = float(ScoreList[i])/NumInfoSites[i]
    out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatios[i], NumMatSNPs))
  out.close()

def match_bed_to_acc(args):
  logging.info("Reading the position file")
  targetSNPs = pandas.read_table(args['inFile'], header=None, usecols=[0,1,2])
  snpCHR = numpy.array(targetSNPs[0], dtype=int)
  snpPOS = numpy.array(targetSNPs[1], dtype=int)
  snpGT = numpy.array(targetSNPs[2])
  GenotypeData = genotype.load_hdf5_genotype_data(args['hdf5File'])
  GenotypeData_acc = genotype.load_hdf5_genotype_data(args['hdf5accFile'])
  num_lines = len(GenotypeData.accessions)
  ScoreList = numpy.zeros(num_lines, dtype="float")
  NumInfoSites = numpy.zeros(len(GenotypeData.accessions), dtype="uint32")
  NumMatSNPs = 0
  chunk_size = 1000
  for i in numpy.array(GenotypeData.chrs, dtype=int):
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
    TarGTs[numpy.where(matchedTarGTs == "0/1")[0]] = 2
    NumMatSNPs = NumMatSNPs + len(matchedAccInd)
    for j in range(0, len(matchedAccInd), chunk_size):
      t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
      samSNPs = numpy.reshape(numpy.repeat(TarGTs[j:j+chunk_size], num_lines), (len(TarGTs[j:j+chunk_size]),num_lines))
      ScoreList = ScoreList + numpy.sum(t1001SNPs == samSNPs, axis=0) 
      if(len(TarGTs[j:j+chunk_size]) > 1):
        NumInfoSites = NumInfoSites + len(TarGTs[j:j+chunk_size]) - numpy.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
      elif(len(TarGTs[j:j+chunk_size]) == 1):
        NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
    logging.info("Done analysing %s positions", NumMatSNPs)
  (LikeLiHoods, LikeLiHoodRatios) = calculate_likelihoods(ScoreList, NumInfoSites)
  print_out_table(args['outFile'],GenotypeData, ScoreList, NumInfoSites, LikeLiHoods, LikeLiHoodRatios, NumMatSNPs)


def match_vcf_to_acc(args):
  logging.info("Reading the VCF file")
  vcf = vcfnp.variants(args['inFile'], cache=False).view(numpy.recarray)
  vcfD = vcfnp.calldata_2d(args['inFile'], cache=False).view(numpy.recarray)
  DPthres = numpy.mean(vcf.DP[numpy.where(vcf.DP > 0)[0]]) * 4
  snpsREQ = numpy.where((vcfD.is_called[:,0]) & (vcf.QUAL > 30) & (vcf.DP > 0) & (vcf.DP < DPthres))[0]
  snpCHR = numpy.array(numpy.char.replace(vcf.CHROM[snpsREQ], "Chr", "")).astype("int8")
  snpPOS = numpy.array(vcf.POS[snpsREQ]) 
  snpPL = vcfD.PL[snpsREQ, 0]
  snpWEI = numpy.copy(snpPL)
  snpWEI = snpWEI.astype(float)
  snpWEI = snpWEI/(-10)
  snpWEI = numpy.exp(snpWEI)
  GenotypeData = genotype.load_hdf5_genotype_data(args['hdf5File'])
  GenotypeData_acc = genotype.load_hdf5_genotype_data(args['hdf5accFile'])
  num_lines = len(GenotypeData.accessions)
  ScoreList = numpy.zeros(num_lines, dtype="float")
  NumInfoSites = numpy.zeros(len(GenotypeData.accessions), dtype="uint32")
  NumMatSNPs = 0
  chunk_size = 1000
  for i in numpy.array(GenotypeData.chrs, dtype=int):
    perchrTarPos = numpy.where(snpCHR == i)[0]
    perchrtarSNPpos = snpPOS[perchrTarPos]
    logging.info("Loaded %s chromosome positions from the position file", i)
    start = GenotypeData.chr_regions[i-1][0]
    end = GenotypeData.chr_regions[i-1][1]
    chrpositions = GenotypeData.positions[start:end]
    matchedAccInd = numpy.where(numpy.in1d(chrpositions, perchrtarSNPpos))[0] + start
    matchedTarInd = numpy.where(numpy.in1d(perchrtarSNPpos, chrpositions))[0]
    matchedTarWei = snpWEI[perchrTarPos[matchedTarInd],]
    TarGTs0 = numpy.zeros(len(matchedTarInd), dtype="int8")
    TarGTs1 = numpy.ones(len(matchedTarInd), dtype="int8") + 1
    TarGTs2 = numpy.ones(len(matchedTarInd), dtype="int8")
    NumMatSNPs = NumMatSNPs + len(matchedAccInd)
    for j in range(0, len(matchedAccInd), chunk_size):
      t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
      samSNPs0 = numpy.reshape(numpy.repeat(TarGTs0[j:j+chunk_size], num_lines), (len(TarGTs0[j:j+chunk_size]),num_lines))
      samSNPs1 = numpy.reshape(numpy.repeat(TarGTs1[j:j+chunk_size], num_lines), (len(TarGTs1[j:j+chunk_size]),num_lines))
      samSNPs2 = numpy.reshape(numpy.repeat(TarGTs2[j:j+chunk_size], num_lines), (len(TarGTs2[j:j+chunk_size]),num_lines))
      tempScore0 = numpy.sum(numpy.multiply(numpy.array(t1001SNPs == samSNPs0, dtype=int).T, matchedTarWei[j:j+chunk_size,0]).T, axis=0)
      tempScore1 = numpy.sum(numpy.multiply(numpy.array(t1001SNPs == samSNPs1, dtype=int).T, matchedTarWei[j:j+chunk_size,1]).T, axis=0)
      tempScore2 = numpy.sum(numpy.multiply(numpy.array(t1001SNPs == samSNPs2, dtype=int).T, matchedTarWei[j:j+chunk_size,2]).T, axis=0)
      ScoreList = ScoreList + tempScore0 + tempScore1 + tempScore2
      if(len(TarGTs0[j:j+chunk_size]) > 1):
        NumInfoSites = NumInfoSites + len(TarGTs0[j:j+chunk_size]) - numpy.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
      elif(len(TarGTs0[j:j+chunk_size]) == 1):
        NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
    logging.info("Done analysing %s positions", NumMatSNPs)
  (LikeLiHoods, LikeLiHoodRatios) = calculate_likelihoods(ScoreList, NumInfoSites)
  print_out_table(args['outFile'], GenotypeData, ScoreList, NumInfoSites, LikeLiHoods, LikeLiHoodRatios, NumMatSNPs)

  

 



