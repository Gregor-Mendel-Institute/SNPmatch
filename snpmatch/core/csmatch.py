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
    NumAmb = np.where(likeliHoodRatio < snpmatch.lr_thres)[0]
    if len(NumAmb) >= 1 and len(NumAmb) < num_lines:
      try:
        nextLikeli = np.nanmin(likeliHoodRatio[np.where(likeliHoodRatio > snpmatch.lr_thres)[0]])
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
  return (TotScoreList, TotNumInfoSites, NumMatSNPs)

def getHomoWindows(likeLiwind):
  snp_thres_wind = np.nanmean(likeLiwind[2]) - np.std(likeLiwind[2])
  x_info = np.unique(likeLiwind[7])
  homo_wind = np.zeros(0, dtype = "int")
  for i in x_info:
    eWinds = likeLiwind.iloc[np.where(likeLiwind[7] == i)[0]]
    if np.nanmean(eWinds[3]) > snpmatch.prob_thres and np.nanmean(eWinds[2]) > snp_thres_wind:
      homo_wind = np.append(homo_wind, i)
  return homo_wind

def crossInterpreter(GenotypeData, binLen, outFile, scoreFile):
  ## ScoreFile should be one from crossF1genotyper
  ## Output file is from the crossIdentifier
  cs_thres = 0.9
  log.info("running cross interpreter!")
  num_lines = len(GenotypeData.accessions)
  likeLiwind = pandas.read_table(outFile, header=None)
  ScoreAcc = pandas.read_table(scoreFile, header=None)
  topHitsDict = json.load(open(scoreFile + ".matches.json", 'r'))
  if topHitsDict['interpretation']['case'] == 3:
    homo_wind = getHomoWindows(likeLiwind)
    homo_acc = np.unique(likeLiwind[0][np.where(np.in1d(likeLiwind[7], homo_wind))[0]],return_counts=True)
    matches_dict = [(homo_acc[0][i].astype("string"), homo_acc[1][i]) for i in np.argsort(-homo_acc[1])]
    topHitsDict['matches'] = matches_dict 
    f1matches = ScoreAcc.iloc[~np.in1d(ScoreAcc[0], GenotypeData.accessions)].reset_index()
    topMatch = np.argsort(f1matches[5])[0]  ## top F1 match sorted based on likelihood
    if f1matches[3][topMatch] > cs_thres:
      mother = f1matches[0][topMatch].split("x")[0]
      father = f1matches[0][topMatch].split("x")[1]
      topHitsDict['interpretation']['text'] = "Sample may be a F1! or a contamination!"
      topHitsDict['interpretation']['case'] = 5
      topHitsDict['parents'] = {'mother': [mother,1], 'father': [father,1]}
      topHitsDict['genotype_windows'] = {'chr_bins': None, 'coordinates': {'x': None, 'y': None}}
    else:
      (ChrBins, PosBins) = getBins(GenotypeData, binLen) 
      ## Get exactly the homozygous windows with one count
      clean = np.unique(likeLiwind[0][np.where(likeLiwind[6] == 1)[0]], return_counts = True)
      if len(clean[0]) > 0:  ## Check if there are atlease one homozygous window
        parents = clean[0][np.argsort(-clean[1])[0:2]].astype("string")
        parents_counts = clean[1][np.argsort(-clean[1])[0:2]].astype("int")
        xdict = np.array(np.unique(likeLiwind[7]), dtype="int")
        ydict = np.repeat("NA", len(xdict)).astype("a25")
        if len(parents) == 1:
          topHitsDict['interpretation']['text'] = "Sample may be a F2! but only one parent found!"
          topHitsDict['interpretation']['case'] = 6
          topHitsDict['parents'] = {'mother': [parents[0], parents_counts[0]], 'father': ["NA", "NA"]}
          par1_ind = likeLiwind[7][np.where((likeLiwind[0].astype("string") == parents[0]) & np.in1d(likeLiwind[7], homo_wind))[0]]
          ydict[np.where(np.in1d(xdict,par1_ind))[0]] = parents[0]
        else:
          topHitsDict['interpretation']['text'] = "Sample may be a F2!"
          topHitsDict['interpretation']['case'] = 6
          topHitsDict['parents'] = {'mother': [parents[0], parents_counts[0]], 'father': [parents[1], parents_counts[1]]}
          NumChrs = np.unique(ChrBins, return_counts=True)
          chr_bins = dict(('Chr%s' % NumChrs[0][i], NumChrs[1][i]) for i in range(len(NumChrs[0])))
          par1_ind = np.array(likeLiwind[7][np.where((likeLiwind[0].astype("string") == parents[0]) & np.in1d(likeLiwind[7], homo_wind))[0]])
          par2_ind = np.array(likeLiwind[7][np.where((likeLiwind[0].astype("string") == parents[1]) & np.in1d(likeLiwind[7], homo_wind))[0]])
          ydict[np.where(np.in1d(xdict,par1_ind))[0]] = parents[0]
          ydict[np.where(np.in1d(xdict,par2_ind))[0]] = parents[1]
        xdict = xdict.tolist()
        ydict = ydict.tolist()
        topHitsDict['genotype_windows'] = {'chr_bins': chr_bins, 'coordinates': {'x': xdict, 'y': ydict}}
      else:   ## No homozygous window found!
        topHitsDict['interpretation']['case'] = 7
        topHitsDict['interpretation']['text'] = "Sample may just be contamination!"
        topHitsDict['genotype_windows'] = {'chr_bins': None, 'coordinates': {'x': None, 'y': None}}
        topHitsDict['parents'] = {'mother': [None,0], 'father': [None,1]}
    with open(scoreFile + ".matches.json", "w") as out_stats:
      out_stats.write(json.dumps(topHitsDict))

def crossIdentifier(binLen, snpCHR, snpPOS, snpWEI, DPmean, GenotypeData, GenotypeData_acc, outFile, scoreFile):
  ## Get tophit accessions
  # sorting based on the final scores
  NumSNPs = len(snpCHR)
  num_lines = len(GenotypeData.accessions)
  (ScoreList, NumInfoSites, NumMatSNPs) = crossWindower(binLen, snpCHR, snpPOS, snpWEI, DPmean, GenotypeData, outFile)
  probScore = [float(ScoreList[i])/NumInfoSites[i] for i in range(num_lines)]
  probScore = np.array(probScore, dtype = float)
  snpmatch.print_topHits(scoreFile + ".matches.json", GenotypeData.accessions, ScoreList, NumInfoSites, float(NumMatSNPs)/NumSNPs, NumMatSNPs)
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
  crossInterpreter(GenotypeData, binLen, outFile, scoreFile)

def potatoCrossIdentifier(args):
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



