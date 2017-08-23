"""
  SNPmatch
"""
import numpy as np
import numpy.ma
from pygwas.core import genotype
import vcfnp
import pandas
import logging
import sys
import os
import StringIO
import json

log = logging.getLogger(__name__)
lr_thres = 3.841
snp_thres = 4000
prob_thres = 0.98

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def likeliTest(n, y):
  p = 0.99999999
  if n > 0 and n != y and y > 0:
    pS = float(y)/n
    a = y * np.log(pS/p)
    b = (n - y) * np.log((1-pS)/(1-p))
    return(a+b)
  elif n == y and n > 0:
    return 1
  elif y == 0:
    return likeliTest(n, y + 1)
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

def print_out_table(outFile, AccList, ScoreList, NumInfoSites, NumMatSNPs, DPmean):
  (LikeLiHoods, LikeLiHoodRatios) = calculate_likelihoods(ScoreList, NumInfoSites)
  if outFile:
    out = open(outFile, 'w')
    for i in range(len(AccList)):
      score = float(ScoreList[i])/NumInfoSites[i]
      out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (AccList[i], int(ScoreList[i]), NumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatios[i], NumMatSNPs, DPmean))
    out.close()
  else:
    for i in range(len(AccList)):
      score = float(ScoreList[i])/NumInfoSites[i]
      sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (AccList[i], int(ScoreList[i]), NumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatios[i], NumMatSNPs, DPmean))

def CaseInterpreter(overlap, NumSNPs, topHits, probScore):
  overlap_thres = 0.5
  num_lines = len(probScore)
  case = 10
  if len(topHits) == 1:
    case = 0
    note = "Unique hit"
  elif np.nanmean(probScore[topHits]) > prob_thres:
    case = 2
    note = "Ambiguous sample: Accessions in top hits can be really close"
  elif overlap > overlap_thres:
    case = 3
    note = "Ambiguous sample: Sample might contain mixture of DNA or contamination"
  elif overlap < overlap_thres:
    case = 4
    note = "Ambiguous sample: Overlap of SNPs is very low, sample may not be in database"
  if case > 4:
    case = 1
    note = "Ambiguous sample"
  return (case, note)

def print_topHits(outFile, AccList, ScoreList, NumInfoSites, overlap, NumMatSNPs):
  num_lines = len(ScoreList)
  (LikeLiHoods, LikeLiHoodRatios) = calculate_likelihoods(ScoreList, NumInfoSites)
  topHits = np.where(LikeLiHoodRatios < lr_thres)[0]
  probScore = [float(ScoreList[i])/NumInfoSites[i] for i in range(num_lines)]
  overlapScore = [float(NumInfoSites[i])/NumMatSNPs for i in range(num_lines)]
  probScore = np.array(probScore, dtype = float)
  sorted_order = topHits[np.argsort(-probScore[topHits])]
  (case, note) = CaseInterpreter(overlap, NumMatSNPs, topHits, probScore)
  matches_dict = [(AccList[i], probScore[i], NumInfoSites[i], overlapScore[i]) for i in sorted_order]
  topHitsDict = {'overlap': [overlap, NumMatSNPs], 'matches': matches_dict, 'interpretation':{'case': case, 'text': note}}
  with open(outFile, "w") as out_stats:
    out_stats.write(json.dumps(topHitsDict))

def parseGT(snpGT):
  first = snpGT[0]
  snpBinary = np.zeros(len(snpGT), dtype = "int8")
  if first.find('|') != -1:
    ## GT is phased
    separator = "|"
  elif first.find('/') != -1:
    ## GT is not phased
    separator = "/"
  elif np.char.isdigit(first):
    return np.array(np.copy(snpGT), dtype = "int8")
  else:
    die("unable to parse the format of GT in vcf!")
  hetGT = "0" + separator + "1"
  refGT = "0" + separator + "0"
  altGT = "1" + separator + "1"
  nocall = "." + separator + "."
  snpBinary[np.where(snpGT == altGT)[0]] = 1
  snpBinary[np.where(snpGT == hetGT)[0]] = 2
  snpBinary[np.where(snpGT == nocall)[0]] = -1
  return snpBinary

def parseChrName(targetCHR):
    snpCHROM = np.char.replace(np.core.defchararray.lower(np.array(targetCHR, dtype="string")), "chr", "")
    snpsREQ = np.where(np.char.isdigit(snpCHROM))[0]   ## Filtering positions from mitochondrial and chloroplast
    snpCHR = snpCHROM[snpsREQ]
    return (snpCHR, snpsREQ)

def readBED(inFile, logDebug):
  log.info("reading the position file")
  targetSNPs = pandas.read_table(inFile, header=None, usecols=[0,1,2])
  (snpCHR, snpsREQ) = parseChrName(targetSNPs[0])
  snpPOS = np.array(targetSNPs[1], dtype=int)[snpsREQ]
  snpGT = np.array(targetSNPs[2])[snpsREQ]
  snpBinary = parseGT(snpGT)
  snpWEI = np.ones((len(snpCHR), 3))  ## for homo and het
  snpWEI[np.where(snpBinary != 0),0] = 0
  snpWEI[np.where(snpBinary != 1),2] = 0
  snpWEI[np.where(snpBinary != 2),1] = 0
  return (snpCHR, snpPOS, snpGT, snpWEI)

def readVcf(inFile, logDebug):
  log.info("reading the VCF file")
  if logDebug:
    vcf = vcfnp.variants(inFile, cache=False).view(np.recarray)
    vcfD = vcfnp.calldata_2d(inFile, cache=False).view(np.recarray)
  else:
    sys.stderr = StringIO.StringIO()
    vcf = vcfnp.variants(inFile, cache=False).view(np.recarray)
    vcfD = vcfnp.calldata_2d(inFile, cache=False).view(np.recarray)
    sys.stderr = sys.__stderr__
  DPthres = np.mean(vcfD.DP[np.where(vcfD.DP > 0)[0]]) * 4
  DPmean = DPthres / 4
  (snpCHR, snpsREQ) = parseChrName(vcf.CHROM)
  snpPOS = np.array(vcf.POS[snpsREQ])
  try:
    snpGT = np.array(vcfD.GT[snpsREQ,0])
  except AttributeError:
    die("input VCF file doesnt have required GT field")
  try:
    snpPL = vcfD.PL[snpsREQ, 0]
    snpWEI = np.copy(snpPL)
    snpWEI = snpWEI.astype(float)
    snpWEI = snpWEI/(-10)
    snpWEI = np.exp(snpWEI)
  except AttributeError:
    snpBinary = parseGT(snpGT)
    snpWEI = np.ones((len(snpsREQ), 3))  ## for homo and het
    snpWEI[np.where(snpBinary != 0),0] = 0
    snpWEI[np.where(snpBinary != 1),2] = 0
    snpWEI[np.where(snpBinary != 2),1] = 0
  return (DPmean, snpCHR, snpPOS, snpGT, snpWEI)

def getHeterozygosity(snpGT, outFile='default'):
    snpBinary = parseGT(snpGT)
    numHets = len(np.where(snpBinary == 2)[0])
    if outFile != 'default':
        with open(outFile) as json_out:
            topHitsDict = json.load(json_out)
        topHitsDict['percent_heterozygosity'] = float(numHets)/len(snpGT)
        with open(outFile, "w") as out_stats:
          out_stats.write(json.dumps(topHitsDict))
    return float(numHets)/len(snpGT)

def parseInput(inFile, logDebug, outFile = "parser"):
  if outFile == "parser" or not outFile:
    outFile = inFile + ".snpmatch"
  if os.path.isfile(inFile + ".snpmatch.npz"):
    log.info("snpmatch parser dump found! loading %s", inFile + ".snpmatch.npz")
    snps = np.load(inFile + ".snpmatch.npz")
    (snpCHR, snpPOS, snpGT, snpWEI, DPmean) = (snps['chr'], snps['pos'], snps['gt'], snps['wei'], snps['dp'])
  else:
    _,inType = os.path.splitext(inFile)
    if inType == '.npz':
      log.info("loading snpmatch parser file! %s", inFile)
      snps = np.load(inFile)
      (snpCHR, snpPOS, snpGT, snpWEI, DPmean) = (snps['chr'], snps['pos'], snps['gt'], snps['wei'], snps['dp'])
    else:
      if inType == '.vcf':
        log.info("creating snpmatch parser dump %s", inFile + ".snpmatch.npz")
        (DPmean, snpCHR, snpPOS, snpGT, snpWEI) = readVcf(inFile, logDebug)
      elif inType == '.bed':
        log.info("creating snpmatch parser dump %s", inFile + ".snpmatch.npz")
        (snpCHR, snpPOS, snpGT, snpWEI) = readBED(inFile, logDebug)
        DPmean = "NA"
      else:
        die("input file type %s not supported" % inType)
      np.savez(outFile, chr = snpCHR, pos = snpPOS, gt = snpGT, wei = snpWEI, dp = DPmean)
      NumSNPs = len(snpCHR)
      case = 0
      note = "Sufficient number of SNPs"
      if NumSNPs < snp_thres:
        note = "Attention: low number of SNPs provided"
        case = 1
      snpst = np.unique(snpCHR, return_counts=True)
      snpdict = dict(('Chr%s' % snpst[0][i], snpst[1][i]) for i in range(len(snpst[0])))
      statdict = {"interpretation": {"case": case, "text": note}, "snps": snpdict, "num_of_snps": NumSNPs, "depth": DPmean}
      statdict['percent_heterozygosity'] = getHeterozygosity(snpGT)
      with open(outFile + ".stats.json", "w") as out_stats:
        out_stats.write(json.dumps(statdict))
  log.info("done!")
  return (snpCHR, snpPOS, snpGT, snpWEI, DPmean)

def genotyper(snpCHR, snpPOS, snpGT, snpWEI, DPmean, hdf5File, hdf5accFile, outFile):
  NumSNPs = len(snpCHR)
  log.info("loading database files")
  GenotypeData = genotype.load_hdf5_genotype_data(hdf5File)
  GenotypeData_acc = genotype.load_hdf5_genotype_data(hdf5accFile)
  log.info("done!")
  num_lines = len(GenotypeData.accessions)
  ScoreList = np.zeros(num_lines, dtype="float")
  NumInfoSites = np.zeros(len(GenotypeData.accessions), dtype="uint32")
  NumMatSNPs = 0
  overlappedInds = np.zeros(0, dtype=int)
  chunk_size = 1000
  for i in parseChrName(GenotypeData.chrs):
    perchrTarPos = np.where(snpCHR == i)[0]
    perchrtarSNPpos = snpPOS[perchrTarPos]
    log.info("Analysing chromosome %s positions", i)
    start = GenotypeData.chr_regions[i-1][0]
    end = GenotypeData.chr_regions[i-1][1]
    chrpositions = GenotypeData.positions[start:end]
    matchedAccInd = np.where(np.in1d(chrpositions, perchrtarSNPpos))[0] + start
    matchedTarInd = np.where(np.in1d(perchrtarSNPpos, chrpositions))[0]
    matchedTarWei = snpWEI[perchrTarPos[matchedTarInd],]
    TarGTs0 = np.zeros(len(matchedTarInd), dtype="int8")
    TarGTs1 = np.ones(len(matchedTarInd), dtype="int8") + 1
    TarGTs2 = np.ones(len(matchedTarInd), dtype="int8")
    overlappedInds = np.append(overlappedInds, perchrTarPos[matchedTarInd])
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
    log.info("Done analysing %s positions", NumMatSNPs)
  log.info("writing score file!")
  overlap = float(NumMatSNPs)/NumSNPs
  print_out_table(outFile + '.scores.txt',GenotypeData.accessions, ScoreList, NumInfoSites, NumMatSNPs, DPmean)
  if not outFile:
    outFile = "genotyper"
  print_topHits(outFile + ".matches.json", GenotypeData.accessions, ScoreList, NumInfoSites, overlap, NumMatSNPs)
  getHeterozygosity(snpGT[overlappedInds], outFile + ".matches.json")
  return (ScoreList, NumInfoSites)

def potatoGenotyper(args):
  (snpCHR, snpPOS, snpGT, snpWEI, DPmean) = parseInput(inFile = args['inFile'], logDebug = args['logDebug'])
  log.info("running genotyper!")
  (ScoreList, NumInfoSites) = genotyper(snpCHR, snpPOS, snpGT, snpWEI, DPmean, args['hdf5File'], args['hdf5accFile'], args['outFile'])
  log.info("finished!")

def pairwiseScore(inFile_1, inFile_2, logDebug, outFile):
    (snpCHR1, snpPOS1, snpGT1, snpWEI1, DPmean1) = parseInput(inFile = inFile_1, logDebug = logDebug)
    (snpCHR2, snpPOS2, snpGT2, snpWEI2, DPmean2) = parseInput(inFile = inFile_2, logDebug = logDebug)
    snpmatch_stats = {}
    unique_1, unique_2, common, scores = 0, 0, 0, 0
    chrs = np.union1d(snpCHR1, snpCHR2)
    for i in chrs:
        perchrTarPosInd1 = np.where(snpCHR1 == i)[0]
        perchrTarPosInd2 = np.where(snpCHR2 == i)[0]
        log.info("Analysing chromosome %s positions", i)
        perchrtarSNPpos1 = snpPOS1[perchrTarPosInd1]
        perchrtarSNPpos2 = snpPOS2[perchrTarPosInd2]
        matchedAccInd1 = np.where(np.in1d(perchrtarSNPpos1, perchrtarSNPpos2))[0]
        matchedAccInd2 = np.where(np.in1d(perchrtarSNPpos2, perchrtarSNPpos1))[0]
        unique_1 = unique_1 + len(perchrTarPosInd1) - len(matchedAccInd1)
        unique_2 = unique_2 + len(perchrTarPosInd2) - len(matchedAccInd2)
        common = common + len(matchedAccInd1)
        scores = scores + np.sum(np.array(snpGT1[matchedAccInd1] == snpGT2[matchedAccInd2], dtype = int))
    snpmatch_stats['unique'] = {"%s" % os.path.basename(inFile_1): [float(unique_1)/len(snpCHR1), len(snpCHR1)], "%s" % os.path.basename(inFile_2): [float(unique_2)/len(snpCHR2), len(snpCHR2)]}
    snpmatch_stats['matches'] = [float(scores)/common, common]
    if not outFile:
        outFile = "genotyper"
    log.info("writing output in a file: %s" % outFile + ".matches.json")
    with open(outFile + ".matches.json", "w") as out_stats:
        out_stats.write(json.dumps(snpmatch_stats))
    log.info("finished!")
