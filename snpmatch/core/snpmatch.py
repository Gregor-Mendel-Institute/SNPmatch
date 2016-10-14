"""
  SNPmatch
"""
import numpy as np
import numpy.ma
from pygwas.core import genotype
import vcfnp
import pandas
import argparse
import logging
import sys
import os
import StringIO


def setLog(args):
  if args['logDebug']:
    numeric_level = getattr(logging, "DEBUG", None)
  else:
    numeric_level = getattr(logging, "CRITICAL", None)
  logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=numeric_level)

def likeliTest(n, y):
  p = 0.999999
  if n > 0:
    pS = float(y)/n
    a = y * np.log(pS/p)
    b = (n - y) * np.log((1-pS)/(1-p))
    return(a+b)
  elif n == y:
    return 1
  else:
    return np.nan

def calculate_likelihoods(ScoreList, NumInfoSites):
  num_lines = len(ScoreList)
  LikeLiHoods = [likeliTest(NumInfoSites[i], int(ScoreList[i])) for i in range(num_lines)]
  LikeLiHoods = np.array(LikeLiHoods).astype("float")
  TopHit = np.amin(LikeLiHoods)
  LikeLiHoodRatios = [LikeLiHoods[i]/TopHit for i in range(num_lines)]
  LikeLiHoodRatios = np.array(LikeLiHoodRatios).astype("float")
  return (LikeLiHoods, LikeLiHoodRatios)

def print_out_table(outFile, GenotypeData, ScoreList, NumInfoSites, LikeLiHoods, LikeLiHoodRatios, NumMatSNPs, DPmean):
  out = open(outFile, 'w')
  for i in range(len(GenotypeData.accessions)):
    score = float(ScoreList[i])/NumInfoSites[i]
    out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatios[i], NumMatSNPs, DPmean))
  out.close()

def match_bed_to_acc(args):
  setLog(args)
  logging.info("Reading the position file")
  targetSNPs = pandas.read_table(args['inFile'], header=None, usecols=[0,1,2])
  snpCHR = np.array(targetSNPs[0], dtype=int)
  snpPOS = np.array(targetSNPs[1], dtype=int)
  snpGT = np.array(targetSNPs[2])
  GenotypeData = genotype.load_hdf5_genotype_data(args['hdf5File'])
  GenotypeData_acc = genotype.load_hdf5_genotype_data(args['hdf5accFile'])
  num_lines = len(GenotypeData.accessions)
  ScoreList = np.zeros(num_lines, dtype="float")
  NumInfoSites = np.zeros(len(GenotypeData.accessions), dtype="uint32")
  NumMatSNPs = 0
  chunk_size = 1000
  for i in np.array(GenotypeData.chrs, dtype=int):
    perchrTarPos = np.where(snpCHR == i)[0]
    perchrtarSNPpos = snpPOS[perchrTarPos]
    logging.info("Loaded %s chromosome positions from the position file", i)
    start = GenotypeData.chr_regions[i-1][0]
    end = GenotypeData.chr_regions[i-1][1]
    chrpositions = GenotypeData.positions[start:end]
    matchedAccInd = np.where(np.in1d(chrpositions, perchrtarSNPpos))[0] + start
    matchedTarInd = np.where(np.in1d(perchrtarSNPpos, chrpositions))[0]
    matchedTarGTs = snpGT[perchrTarPos[matchedTarInd],]
    TarGTs = np.zeros(len(matchedTarGTs), dtype="int8")
    TarGTs[np.where(matchedTarGTs == "1/1")[0]] = 1
    TarGTs[np.where(matchedTarGTs == "0/1")[0]] = 2
    NumMatSNPs = NumMatSNPs + len(matchedAccInd)
    for j in range(0, len(matchedAccInd), chunk_size):
      t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
      samSNPs = np.reshape(np.repeat(TarGTs[j:j+chunk_size], num_lines), (len(TarGTs[j:j+chunk_size]),num_lines))
      ScoreList = ScoreList + np.sum(t1001SNPs == samSNPs, axis=0) 
      if(len(TarGTs[j:j+chunk_size]) > 1):
        NumInfoSites = NumInfoSites + len(TarGTs[j:j+chunk_size]) - np.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
      elif(len(TarGTs[j:j+chunk_size]) == 1):
        NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
    logging.info("Done analysing %s positions", NumMatSNPs)
  (LikeLiHoods, LikeLiHoodRatios) = calculate_likelihoods(ScoreList, NumInfoSites)
  if args['outFile']:
    print_out_table(args['outFile'],GenotypeData, ScoreList, NumInfoSites, LikeLiHoods, LikeLiHoodRatios, NumMatSNPs, DPmean, num_lines)
  else:
    for i in range(len(GenotypeData.accessions)):
      score = float(ScoreList[i])/NumInfoSites[i]
      sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatios[i], NumMatSNPs, DPmean))


def match_vcf_to_acc(args):
  setLog(args)
  logging.info("Reading the VCF file")
  if args['logDebug']:
    vcf = vcfnp.variants(args['inFile'], cache=True).view(np.recarray)
    vcfD = vcfnp.calldata_2d(args['inFile'], cache=True).view(np.recarray)
  else:
    sys.stderr = StringIO.StringIO()
    vcf = vcfnp.variants(args['inFile'], cache=True).view(np.recarray)
    vcfD = vcfnp.calldata_2d(args['inFile'], cache=True).view(np.recarray)
    sys.stderr = sys.__stderr__
  logging.info("done!")
  DPthres = np.mean(vcf.DP[np.where(vcf.DP > 0)[0]]) * 4
  DPmean = DPthres/4
  snpCHROM =  np.char.replace(vcf.CHROM, "Chr", "")
  snpsREQ = np.where((vcfD.is_called[:,0]) & (vcf.QUAL > 30) & (vcf.DP > 0) & (vcf.DP < DPthres) & (np.char.isdigit(snpCHROM)))[0]
  snpCHR = np.array(snpCHROM[snpsREQ]).astype("int8")
  snpPOS = np.array(vcf.POS[snpsREQ]) 
  snpPL = vcfD.PL[snpsREQ, 0]
  snpWEI = np.copy(snpPL)
  snpWEI = snpWEI.astype(float)
  snpWEI = snpWEI/(-10)
  snpWEI = np.exp(snpWEI)
  logging.info("loading database files")
  GenotypeData = genotype.load_hdf5_genotype_data(args['hdf5File'])
  GenotypeData_acc = genotype.load_hdf5_genotype_data(args['hdf5accFile'])
  logging.info("done!")
  num_lines = len(GenotypeData.accessions)
  ScoreList = np.zeros(num_lines, dtype="float")
  NumInfoSites = np.zeros(len(GenotypeData.accessions), dtype="uint32")
  NumMatSNPs = 0
  chunk_size = 1000
  for i in np.array(GenotypeData.chrs, dtype=int):
    perchrTarPos = np.where(snpCHR == i)[0]
    perchrtarSNPpos = snpPOS[perchrTarPos]
    logging.info("Loaded %s chromosome positions from the position file", i)
    start = GenotypeData.chr_regions[i-1][0]
    end = GenotypeData.chr_regions[i-1][1]
    chrpositions = GenotypeData.positions[start:end]
    matchedAccInd = np.where(np.in1d(chrpositions, perchrtarSNPpos))[0] + start
    matchedTarInd = np.where(np.in1d(perchrtarSNPpos, chrpositions))[0]
    matchedTarWei = snpWEI[perchrTarPos[matchedTarInd],]
    TarGTs0 = np.zeros(len(matchedTarInd), dtype="int8")
    TarGTs1 = np.ones(len(matchedTarInd), dtype="int8") + 1
    TarGTs2 = np.ones(len(matchedTarInd), dtype="int8")
    NumMatSNPs = NumMatSNPs + len(matchedAccInd)
    for j in range(0, len(matchedAccInd), chunk_size):
      t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
      samSNPs0 = np.reshape(np.repeat(TarGTs0[j:j+chunk_size], num_lines), (len(TarGTs0[j:j+chunk_size]),num_lines))
      samSNPs1 = np.reshape(np.repeat(TarGTs1[j:j+chunk_size], num_lines), (len(TarGTs1[j:j+chunk_size]),num_lines))
      samSNPs2 = np.reshape(np.repeat(TarGTs2[j:j+chunk_size], num_lines), (len(TarGTs2[j:j+chunk_size]),num_lines))
      tempScore0 = np.sum(np.multiply(np.array(t1001SNPs == samSNPs0, dtype=int).T, matchedTarWei[j:j+chunk_size,0]).T, axis=0)
      tempScore1 = np.sum(np.multiply(np.array(t1001SNPs == samSNPs1, dtype=int).T, matchedTarWei[j:j+chunk_size,1]).T, axis=0)
      tempScore2 = np.sum(np.multiply(np.array(t1001SNPs == samSNPs2, dtype=int).T, matchedTarWei[j:j+chunk_size,2]).T, axis=0)
      ScoreList = ScoreList + tempScore0 + tempScore1 + tempScore2
      if(len(TarGTs0[j:j+chunk_size]) > 1):
        NumInfoSites = NumInfoSites + len(TarGTs0[j:j+chunk_size]) - np.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
      elif(len(TarGTs0[j:j+chunk_size]) == 1):
        NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
    logging.info("Done analysing %s positions", NumMatSNPs)
  (LikeLiHoods, LikeLiHoodRatios) = calculate_likelihoods(ScoreList, NumInfoSites)
  if args['outFile']:
    print_out_table(args['outFile'],GenotypeData, ScoreList, NumInfoSites, LikeLiHoods, LikeLiHoodRatios, NumMatSNPs, DPmean)
  else:
    for i in range(len(GenotypeData.accessions)):
      score = float(ScoreList[i])/NumInfoSites[i]
      sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatios[i], NumMatSNPs, DPmean))

  

 



