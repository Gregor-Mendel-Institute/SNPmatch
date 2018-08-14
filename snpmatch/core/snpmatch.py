"""
  SNPmatch
"""
import numpy as np
import numpy.ma
from pygwas.core import genotype
import logging
import sys
import os
import parsers
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

class GenotyperOutput(object):
    ## class object for main SNPmatch output

    def __init__(self, AccList, ScoreList, NumInfoSites, overlap, NumMatSNPs, DPmean ):
        self.accs = np.array(AccList, dtype="string")
        self.scores = np.array(ScoreList, dtype="int")
        self.ninfo = np.array(NumInfoSites, dtype="int")
        self.overlap = overlap
        self.num_snps = NumMatSNPs
        self.dp = DPmean

    def get_probabilities(self):
        probs = [float(self.scores[i])/self.ninfo[i] for i in range(len(self.accs))]
        self.probabilies = np.array(probs, dtype="float")

    @staticmethod
    def calculate_likelihoods(scores, ninfo, amin = "calc"):
        num_lines = len(scores)
        LikeLiHoods = [likeliTest(ninfo[i], int(scores[i])) for i in range(num_lines)]
        LikeLiHoods = np.array(LikeLiHoods, dtype = "float")
        if amin == "calc":
            TopHit = np.amin(LikeLiHoods)
        else:
            TopHit = float(amin)
        LikeLiHoodRatios = [LikeLiHoods[i]/TopHit for i in range(num_lines)]
        LikeLiHoodRatios = np.array(LikeLiHoodRatios, dtype="float")
        return((LikeLiHoods, LikeLiHoodRatios))

    def get_likelihoods(self, amin = "calc"):
        (self.likelis, self.lrts) = self.calculate_likelihoods(self.scores, self.ninfo, amin)

    def print_out_table(self, outFile):
        self.get_likelihoods()
        self.get_probabilities()
        if outFile:
            out = open(outFile, 'w')
            for i in range(len(self.accs)):
                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.accs[i], self.scores[i], self.ninfo[i], self.probabilies[i], self.likelis[i], self.lrts[i], self.num_snps, self.dp))
            out.close()
        else:
            for i in range(len(self.accs)):
                sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.accs[i], self.scores[i], self.ninfo[i], self.probabilies[i], self.likelis[i], self.lrts[i], self.num_snps, self.dp))

    def print_json_output(self, outFile):
        self.get_likelihoods()
        self.get_probabilities()
        topHits = np.where(self.lrts < lr_thres)[0]
        overlapScore = [float(self.ninfo[i])/self.num_snps for i in range(len(self.accs))]
        sorted_order = topHits[np.argsort(-self.probabilies[topHits])]
        (case, note) = self.case_interpreter(topHits)
        matches_dict = [(self.accs[i], self.probabilies[i], self.ninfo[i], overlapScore[i]) for i in sorted_order]
        topHitsDict = {'overlap': [self.overlap, self.num_snps], 'matches': matches_dict, 'interpretation':{'case': case, 'text': note}}
        with open(outFile, "w") as out_stats:
            out_stats.write(json.dumps(topHitsDict))

    def case_interpreter(self, topHits):
        overlap_thres = 0.5
        case = 1
        note = "Ambiguous sample"
        if len(topHits) == 1:
            case = 0
            note = "Unique hit"
        elif np.nanmean(self.probabilies[topHits]) > prob_thres:
            case = 2
            note = "Ambiguous sample: Accessions in top hits can be really close"
        elif self.overlap > overlap_thres:
            case = 3
            note = "Ambiguous sample: Sample might contain mixture of DNA or contamination"
        elif self.overlap < overlap_thres:
            case = 4
            note = "Ambiguous sample: Overlap of SNPs is very low, sample may not be in database"
        return(case, note)


def getHeterozygosity(snpGT, outFile='default'):
    snpBinary = parsers.parseGT(snpGT)
    numHets = len(np.where(snpBinary == 2)[0])
    if outFile != 'default':
        with open(outFile) as json_out:
            topHitsDict = json.load(json_out)
        topHitsDict['percent_heterozygosity'] = float(numHets)/len(snpGT)
        with open(outFile, "w") as out_stats:
          out_stats.write(json.dumps(topHitsDict))
    return float(numHets)/len(snpGT)

def genotyper(inputs, hdf5File, hdf5accFile, outFile):
    log.info("loading database files")
    GenotypeData = genotype.load_hdf5_genotype_data(hdf5File)
    GenotypeData_acc = genotype.load_hdf5_genotype_data(hdf5accFile)
    log.info("done!")
    inputs.filter_chr_names(GenotypeData)
    num_lines = len(GenotypeData.accessions)
    ScoreList = np.zeros(num_lines, dtype="float")
    NumInfoSites = np.zeros(len(GenotypeData.accessions), dtype="uint32")
    NumMatSNPs = 0
    overlappedInds = np.zeros(0, dtype=int)
    chunk_size = 1000
    for ind,echr in enumerate(inputs.chr_list):
        perchrTarPos = np.where(inputs.chrs_nochr == echr)[0]
        perchrtarSNPpos = inputs.pos[perchrTarPos]
        log.info("Analysing positions in chromosome %s", echr)
        start = GenotypeData.chr_regions[ind][0]
        end = GenotypeData.chr_regions[ind][1]
        chrpositions = GenotypeData.positions[start:end]
        matchedAccInd = np.where(np.in1d(chrpositions, perchrtarSNPpos))[0] + start
        matchedTarInd = np.where(np.in1d(perchrtarSNPpos, chrpositions))[0]
        matchedTarWei = inputs.wei[perchrTarPos[matchedTarInd],]
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
            if(len(TarGTs0[j:j+chunk_size]) >= 1):
                NumInfoSites = NumInfoSites + len(TarGTs0[j:j+chunk_size]) - np.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
        log.info("Done analysing %s positions", NumMatSNPs)
    log.info("writing score file!")
    overlap = float(NumMatSNPs)/len(inputs.filter_inds_ix)
    result = GenotyperOutput(GenotypeData.accessions, ScoreList, NumInfoSites, overlap, NumMatSNPs, inputs.dp)
    result.print_out_table( outFile + '.scores.txt' )
    if not outFile:
        outFile = "genotyper"
    result.print_json_output( outFile + ".matches.json" )
    getHeterozygosity(inputs.gt[overlappedInds], outFile + ".matches.json")
    return(result)

def potatoGenotyper(args):
    inputs = parsers.ParseInputs(inFile = args['inFile'], logDebug = args['logDebug'])
    log.info("running genotyper!")
    result = genotyper(inputs, args['hdf5File'], args['hdf5accFile'], args['outFile'])
    log.info("finished!")

def pairwiseScore(inFile_1, inFile_2, logDebug, outFile):
    inputs_1 = parsers.ParseInputs(inFile = inFile_1, logDebug = logDebug)
    inputs_2 = parsers.ParseInputs(inFile = inFile_2, logDebug = logDebug)
    snpmatch_stats = {}
    unique_1, unique_2, common, scores = 0, 0, 0, 0
    common_chrs = np.intersect1d(np.unique(inputs_1.chrs), np.unique(inputs_2.chrs))
    for i in common_chrs:
        perchrTarPosInd1 = np.where(inputs_1.chrs == i)[0]
        perchrTarPosInd2 = np.where(inputs_2.chrs == i)[0]
        log.info("Analysing chromosome %s positions", i)
        perchrtarSNPpos1 = inputs_1.pos[perchrTarPosInd1]
        perchrtarSNPpos2 = inputs_2.pos[perchrTarPosInd2]
        matchedAccInd1 = np.where(np.in1d(perchrtarSNPpos1, perchrtarSNPpos2))[0]
        matchedAccInd2 = np.where(np.in1d(perchrtarSNPpos2, perchrtarSNPpos1))[0]
        unique_1 = unique_1 + len(perchrTarPosInd1) - len(matchedAccInd1)
        unique_2 = unique_2 + len(perchrTarPosInd2) - len(matchedAccInd2)
        common = common + len(matchedAccInd1)
        scores = scores + np.sum(np.array(inputs_1.gt[matchedAccInd1] == inputs_2.gt[matchedAccInd2], dtype = int))
    snpmatch_stats['unique'] = {"%s" % os.path.basename(inFile_1): [float(unique_1)/len(inputs_1.chrs), len(inputs_1.chrs)], "%s" % os.path.basename(inFile_2): [float(unique_2)/len(inputs_2.chrs), len(inputs_2.chrs)]}
    snpmatch_stats['matches'] = [float(scores)/common, common]
    if not outFile:
        outFile = "genotyper"
    log.info("writing output in a file: %s" % outFile + ".matches.json")
    with open(outFile + ".matches.json", "w") as out_stats:
        out_stats.write(json.dumps(snpmatch_stats))
    log.info("finished!")
