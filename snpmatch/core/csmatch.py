"""
  SNPmatch for crosses (F2s and F3s)
"""
import numpy as np
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

def getBins(g, binLen):
  binLen = int(binLen)
  chrlen = np.array((30427671, 19698289, 23459830, 18585056, 26975502))  # Arabidopsis chromosome length
  ChrIndex = np.zeros(0,dtype="int8")
  LenBins = np.zeros(0, dtype="int16")
  for i in np.array(g.chrs, dtype=int):
    start = g.chr_regions[i-1][0]
    end = g.chr_regions[i-1][1]
    chr_pos = g.positions[start:end]
    for j in range(0, chrlen[i-1], binLen):
      bin_pos = len(np.where((chr_pos >= j) & (chr_pos < j + binLen))[0])
      LenBins = np.append(LenBins, bin_pos)
      ChrIndex = np.append(ChrIndex, i)
  return(ChrIndex, LenBins)

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

def print_out_table(outFile, GenotypeData, ScoreList, NumInfoSites, LikeLiHoods, LikeLiHoodRatios, NumMatSNPs):
  out = open(outFile, 'w')
  for i in range(len(GenotypeData.accessions)):
    score = float(ScoreList[i])/NumInfoSites[i]
    out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatios[i], NumMatSNPs))
  out.close()


def match_bed_to_acc(args):
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
  lr_thres = 3.841
  TotScoreList = np.zeros(num_lines, dtype="uint32")
  TotNumInfoSites = np.zeros(num_lines, dtype="uint32")
  (ChrBins, PosBins) = getBins(GenotypeData, args['binLen'])
  outfile = open(args['scoreFile'], 'w')
  for i in range(len(PosBins)):
    start = np.sum(PosBins[0:i])
    end = start + PosBins[i]
    pos = GenotypeData.positions[start:end]
    perchrTarPos = np.where(snpCHR == ChrBins[i])[0]
    perchrtarSNPpos = snpPOS[perchrTarPos]
    matchedAccInd = np.where(np.in1d(pos, perchrtarSNPpos))[0] + start
    matchedTarInd = np.where(np.in1d(perchrtarSNPpos, pos))[0]
    matchedTarGTs = snpGT[perchrTarPos[matchedTarInd],]
    TarGTs = np.zeros(len(matchedTarGTs), dtype="int8")
    TarGTs[np.where(matchedTarGTs == "1/1")[0]] = 1
    TarGTs[np.where(matchedTarGTs == "0/1")[0]] = 2
    NumMatSNPs = NumMatSNPs + len(matchedAccInd)
    ScoreList = np.zeros(num_lines, dtype="uint32")
    NumInfoSites = np.zeros(num_lines, dtype="uint32")
    for j in range(0, len(matchedAccInd), chunk_size):
      t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
      samSNPs = np.reshape(np.repeat(TarGTs[j:j+chunk_size], num_lines), (len(TarGTs[j:j+chunk_size]),num_lines))
      ScoreList = ScoreList + np.sum(t1001SNPs == samSNPs, axis=0)
      if(len(TarGTs[j:j+chunk_size]) > 1):
        NumInfoSites = NumInfoSites + len(TarGTs[j:j+chunk_size]) - np.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
      elif(len(TarGTs[j:j+chunk_size]) == 1):
        NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
    if i % 10 == 0:
      logging.info("Done analysing %s positions", NumMatSNPs)
    TotScoreList = TotScoreList + ScoreList
    TotNumInfoSites = TotNumInfoSites + NumInfoSites
    likeliScore = [likeliTest(NumInfoSites[k], ScoreList[k]) for k in range(num_lines)]
    likeliScore = np.array(likeliScore)
    if len(likeliScore) > 0:
      TopHit = np.nanmin(likeliScore)
      likeliHoodRatio = [likeliScore[k]/TopHit for k in range(num_lines)]
      likeliHoodRatio = np.array(likeliHoodRatio)
      TopHitAcc = np.where(likeliHoodRatio == 1)[0]
      NumAmb = np.where(likeliHoodRatio < lr_thres)[0]
      if len(TopHitAcc) == 1:
        t = TopHitAcc[0]
        score = float(ScoreList[t])/NumInfoSites[t]
        if len(np.where(likeliHoodRatio >= lr_thres)[0]) > 0:
          nextLikeli = np.nanmin(likeliHoodRatio[np.where(likeliHoodRatio >= lr_thres)[0]])
          outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[t], ScoreList[t], NumInfoSites[t], score, likeliScore[t], nextLikeli, len(NumAmb), i+1))
        else:
          nextLikeli = np.nanmin(likeliHoodRatio[np.where(likeliHoodRatio >= 1)[0]])
          outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[t], ScoreList[t], NumInfoSites[t], score, likeliScore[t], nextLikeli, len(NumAmb), i+1))
      elif len(TopHitAcc) > 1:
        if len(np.where(likeliHoodRatio >= lr_thres)[0]) > 0:
          nextLikeli = np.nanmin(likeliHoodRatio[np.where(likeliHoodRatio >= lr_thres)[0]])
          for k in range(len(TopHitAcc)):
            t = TopHitAcc[k]
            score = float(ScoreList[t])/NumInfoSites[t]
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[t], ScoreList[t], NumInfoSites[t], score, likeliScore[t], nextLikeli, len(NumAmb), i+1))
        else:
          nextLikeli = np.nanmin(likeliHoodRatio[np.where(likeliHoodRatio >= 1)[0]])
          for k in range(len(TopHitAcc)):
            t = TopHitAcc[k]
            score = float(ScoreList[t])/NumInfoSites[t]
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[t], ScoreList[t], NumInfoSites[t], score, likeliScore[t], nextLikeli, len(NumAmb), i+1))
  (LikeLiHoods, LikeLiHoodRatios) = calculate_likelihoods(TotScoreList, TotNumInfoSites)
  print_out_table(args['outFile'],GenotypeData, TotScoreList, TotNumInfoSites, LikeLiHoods, LikeLiHoodRatios, NumMatSNPs)




