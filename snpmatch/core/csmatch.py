"""
  SNPmatch for crosses (F2s and F3s)
"""
import numpy as np
import numpy.ma
import scipy.stats as st
from pygwas.core import genotype
import pandas
import logging
import os
import StringIO
import snpmatch
import json
import itertools

log = logging.getLogger(__name__)

lr_thres = 3.841

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

def writeBinData(out_file, i, GenotypeData, ScoreList, NumInfoSites):
  num_lines = len(GenotypeData.accessions)
  (likeliScore, likeliHoodRatio) = snpmatch.calculate_likelihoods(ScoreList, NumInfoSites)
  if len(likeliScore) > 0:
    NumAmb = np.where(likeliHoodRatio < lr_thres)[0]
    if len(NumAmb) >= 1 and len(NumAmb) < num_lines:
      try:
        nextLikeli = np.nanmin(likeliHoodRatio[np.where(likeliHoodRatio > lr_thres)[0]])
      except:
        nextLikeli = 1
      for k in NumAmb:
        score = float(ScoreList[k])/NumInfoSites[k]
        out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[k], int(ScoreList[k]), NumInfoSites[k], score, likeliScore[k], nextLikeli, len(NumAmb), i+1))

def crossWindower(binLen, snpCHR, snpPOS, snpWEI, DPmean, GenotypeData, outFile):
  NumSNPs = len(snpCHR)
  num_lines = len(GenotypeData.accessions)
  NumMatSNPs = 0
  chunk_size = 1000
  TotScoreList = np.zeros(num_lines, dtype="uint32")
  TotNumInfoSites = np.zeros(num_lines, dtype="uint32")
  (ChrBins, PosBins) = getBins(GenotypeData, binLen)
  out_file = open(outFile, 'w')
  for i in range(len(PosBins)):
    start = np.sum(PosBins[0:i])
    end = start + PosBins[i]
    pos = GenotypeData.positions[start:end]
    perchrTarPos = np.where(snpCHR == ChrBins[i])[0]
    perchrtarSNPpos = snpPOS[perchrTarPos]
    matchedAccInd = np.where(np.in1d(pos, perchrtarSNPpos))[0] + start
    matchedTarInd = np.where(np.in1d(perchrtarSNPpos, pos))[0]
    matchedTarWei = snpWEI[perchrTarPos[matchedTarInd],]
    TarGTs0 = np.zeros(len(matchedTarInd), dtype="int8")
    TarGTs1 = np.ones(len(matchedTarInd), dtype="int8") + 1
    TarGTs2 = np.ones(len(matchedTarInd), dtype="int8")
    NumMatSNPs = NumMatSNPs + len(matchedAccInd)
    ScoreList = np.zeros(num_lines, dtype="uint32")
    NumInfoSites = np.zeros(num_lines, dtype="uint32")
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
        NumInfoSites = NumInfoSites + len(TarGTs0[j:j+chunk_size]) - np.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0)
      elif(len(TarGTs0[j:j+chunk_size]) == 1):
        NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
    TotScoreList = TotScoreList + ScoreList
    TotNumInfoSites = TotNumInfoSites + NumInfoSites
    writeBinData(out_file, i, GenotypeData, ScoreList, NumInfoSites)
    if i % 50 == 0:
      log.info("Done analysing %s positions", NumMatSNPs)
  out_file.close()
  return (TotScoreList, TotNumInfoSites)

def crossInterpreter(outFile, snpmatchJson):
  ## ScoreFile should be one from crossF1genotyper
  ## Output file is from the crossIdentifier
  import operator
  log.info("running cross interpreter!")
  likeLiwind = pandas.read_table(outFile, header=None)
  topHitsDict = json.load(open(snpmatchJson, 'r'))
  topHitsDict['interpreter'] = "csmatch"
  try:
    parents = sorted(topHitsDict['matches'].items(), key=operator.itemgetter(1),reverse=True)[0][0].encode("ascii")
    mother = parents.split("x")[0]
    father = parents.split("x")[1]
    topHitsDict['interpretation']['text'] = "Sample may be and F1 or contaminated!"
    topHitsDict['parents'] = {'mother': [mother,1], 'father': [father,1]}
  except:
    homowind = np.where(likeLiwind[6] == 1)[0]
    clean = np.unique(likeLiwind[0][homowind], return_counts = True)
    homAcc = clean[0][np.argsort(-clean[1])[0:2]]
    homAcccounts = clean[1][np.argsort(-clean[1])[0:2]]
    topHitsDict['interpretation']['text'] = "Sample is be a cross!"
    topHitsDict['parents'] = {'mother': [homAcc[0], homAcccounts[0]], 'father': [homAcc[1], homAcccounts[1]]}
  with open(snpmatchJson, "w") as out_stats:
    out_stats.write(json.dumps(topHitsDict))

def crossIdentifier(binLen, snpCHR, snpPOS, snpWEI, DPmean, GenotypeData, GenotypeData_acc, outFile, scoreFile):
  ## Get tophit accessions
  # sorting based on the final scores
  NumSNPs = len(snpCHR)
  num_lines = len(GenotypeData.accessions)
  (ScoreList, NumInfoSites) = crossWindower(binLen, snpCHR, snpPOS, snpWEI, DPmean, GenotypeData, outFile)
  probScore = [float(ScoreList[i])/NumInfoSites[i] for i in range(num_lines)]
  probScore = np.array(probScore, dtype = float)
  log.info("simulating F1s for top 10 accessions")
  Accessions = np.copy(GenotypeData.accessions)
  TopHitAccs = np.argsort(-probScore)[0:10]
  for (i, j) in itertools.combinations(TopHitAccs, 2):
    p1 = GenotypeData_acc.snps[:,i]
    p2 = GenotypeData_acc.snps[:,j]
    score = 0
    numinfo = 0
    NumMatSNPs = 0
    for c in np.array(GenotypeData.chrs, dtype=int):
      perchrTarPos = np.where(snpCHR == c)[0]
      perchrtarSNPpos = snpPOS[perchrTarPos]
      start = GenotypeData.chr_regions[c-1][0]
      end = GenotypeData.chr_regions[c-1][1]
      chrpositions = GenotypeData.positions[start:end]
      matchedAccInd = np.where(np.in1d(chrpositions, perchrtarSNPpos))[0] + start
      matchedTarInd = np.where(np.in1d(perchrtarSNPpos, chrpositions))[0]
      NumMatSNPs = NumMatSNPs + len(matchedTarInd)
      gtp1 = p1[matchedAccInd]
      gtp2 = p2[matchedAccInd]
      matchedTarWEI = snpWEI[perchrTarPos[matchedTarInd],]
      homalt = np.where((gtp1 == 1) & (gtp2 == 1))[0]
      homref = np.where((gtp1 == 0) & (gtp2 == 0))[0]
      het = np.where((gtp1 != -1) & (gtp2 != -1) & (gtp1 != gtp2))[0]
      score = score + np.sum(matchedTarWEI[homalt, 2]) + np.sum(matchedTarWEI[homref, 0]) + np.sum(matchedTarWEI[het, 1])
      numinfo = numinfo + len(homalt) + len(homref) + len(het)
    ScoreList = np.append(ScoreList, score)
    NumInfoSites = np.append(NumInfoSites, numinfo)
    acc = GenotypeData.accessions[i] + "x" + GenotypeData.accessions[j]
    Accessions = np.append(Accessions, acc)
  log.info("writing output!")
  snpmatch.print_out_table(scoreFile, Accessions, ScoreList, NumInfoSites, NumMatSNPs, DPmean)
  snpmatch.print_topHits(scoreFile + ".matches.json", Accessions, ScoreList, NumInfoSites, float(NumMatSNPs)/NumSNPs, NumMatSNPs)
  crossInterpreter(outFile, scoreFile + ".matches.json")

def potatoCrossIdentifier(args):
  log.info("running SNPmatch parser!")
  (snpCHR, snpPOS, snpGT, snpWEI, DPmean) = snpmatch.parseInput(inFile = args['inFile'], logDebug = args['logDebug'])
  log.info("loading genotype files!")
  GenotypeData = genotype.load_hdf5_genotype_data(args['hdf5File'])
  GenotypeData_acc = genotype.load_hdf5_genotype_data(args['hdf5accFile'])
  log.info("done!")
  log.info("running cross identifier!")
  crossIdentifier(args['binLen'],snpCHR, snpPOS, snpWEI, DPmean, GenotypeData, GenotypeData_acc, args['outFile'], args['scoreFile'])
  log.info("finished!")

def crossGenotyper(args):
  ## Get the VCF file (filtered may be) generated by GATK.
  # inputs:
  # 1) VCF file
  # 2) Parent1 and Parent2
  # 3) SNP matrix (hdf5 file)
  # 4) Bin length, default as 200Kbp
  # 5) Chromosome length
  (snpCHR, snpPOS, snpGT, snpWEI, DPmean) = snpmatch.parseInput(inFile = args['inFile'], logDebug = args['logDebug'])
  parents = args['parents']
  ## need to filter the SNPs present in C and M 
  log.info("loading HDF5 file")
  GenotypeData_acc = genotype.load_hdf5_genotype_data(args['hdf5accFile'])  
  ## die if either parents are not in the dataset
  try:
    indP1 = np.where(GenotypeData_acc.accessions == parents.split("x")[0])[0][0]
    indP2 = np.where(GenotypeData_acc.accessions == parents.split("x")[1])[0][0]
  except:
    snpmatch.die("parents are not in the dataset")
  snpsP1 = GenotypeData_acc.snps[:,indP1]
  snpsP2 = GenotypeData_acc.snps[:,indP2]
  # identifying the segregating SNPs between the accessions
  # only selecting 0 or 1
  segSNPsind = np.where((snpsP1 != snpsP2) & (snpsP1 >= 0) & (snpsP2 >= 0) & (snpsP1 < 2) & (snpsP2 < 2))[0]
  log.info("number of segregating snps between parents: %s", len(segSNPsind))
  (ChrBins, PosBins) = getBins(GenotypeData_acc, args['binLen'])
  log.info("number of bins: %s", len(ChrBins))
  outfile = open(args['outFile'], 'w')
  for i in range(len(PosBins)):
    start = np.sum(PosBins[0:i])
    end = start + PosBins[i]
    # first snp positions which are segregating and are in this window
    reqPOSind = segSNPsind[np.where((segSNPsind < end) & (segSNPsind >= start))[0]]
    reqPOS = GenotypeData_acc.positions[reqPOSind]
    perchrTarPosind = np.where(snpCHR == ChrBins[i])[0]
    perchrTarPos = snpPOS[perchrTarPosind]
    matchedAccInd = reqPOSind[np.where(np.in1d(reqPOS, perchrTarPos))[0]]
    matchedTarInd = perchrTarPosind[np.where(np.in1d(perchrTarPos, reqPOS))[0]]
    matchedTarGTs = snpGT[matchedTarInd]
    try:
      TarGTs = snpmatch.parseGT(matchedTarGTs)
      TarGTs[np.where(TarGTs == 2)[0]] = 4
      genP1 = np.subtract(TarGTs, snpsP1[matchedAccInd])
      genP1no = len(np.where(genP1 == 0)[0])
      if len(genP1) > 0:
        pValP1 = st.binom_test(genP1no, len(genP1), 0.8, alternative = "greater")
        pValP2 = st.binom_test(len(genP1) - genP1no, len(genP1), 0.8, alternative = "greater")
        if pValP1 < 0.05:
          outfile.write("%s\t%s\t%s\t0\t%s\n" % (i+1, genP1no, len(genP1), pValP1))
        elif pValP2 < 0.05:
          outfile.write("%s\t%s\t%s\t1\t%s\n" % (i+1, genP1no, len(genP1), pValP2))
        elif float(genP1no)/len(genP1) >= 0.8 or float(genP1no)/len(genP1) <= 0.2:
          outfile.write("%s\t%s\t%s\tNA\tNA\n" % (i+1, genP1no, len(genP1)))
        else:
          outfile.write("%s\t%s\t%s\t0.5\tNA\n" % (i+1, genP1no, len(genP1)))
      else:
        outfile.write("%s\t%s\t%s\tNA\tNA\n" % (i+1, genP1no, len(genP1)))
    except:
      outfile.write("%s\tNA\tNA\tNA\tNA\n" % (i+1))
    if i % 10 == 0:
      log.info("progress: %s windows", i+10)
  log.info("done!")
  outfile.close()



