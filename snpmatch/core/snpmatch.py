"""
  SNPmatch
"""
import numpy as np
import scipy as sp
from scipy import stats
import numpy.ma
import logging
import sys
import os
from . import parsers
from . import snp_genotype
import json

log = logging.getLogger(__name__)
lr_thres = 3.841
snp_thres = 4000
prob_thres = 0.98

def die(msg):
    sys.stderr.write('Error: ' + msg + '\n')
    sys.exit(1)

def get_fraction(x, y):
    if y == 0:
        return(np.nan)
    return(float(x)/y)

def likeliTest(n, y):
    ## n == total informative sites
    ## y == number of matched sites
    assert y <= n, "provided y is greater than n"
    p = 0.99999999
    if n == 0:
        return(np.nan)
    pS = float(y)/n
    if y == n:
        return(1)
    if y > 0:
        a = y * np.log(pS/p)
        b = (n - y) * np.log((1-pS)/(1-p))
        return(a+b)
    elif y == 0:
        return(np.nan)

def test_identity(x, n, error_rate = 0.0005, pthres = 0.05):
    ## error rate is divided by 100 since it is considered to be a percentage
    p = 1 - float(error_rate)
    st = stats.binom_test(x, n, p = p, alternative='less')
    if st > pthres:
        return(float(1))
    else:
        return(float(0))

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
        probs = [get_fraction(self.scores[i], self.ninfo[i]) for i in range(len(self.accs))]
        self.probabilies = np.array(probs, dtype="float")

    @staticmethod
    def calculate_likelihoods(scores, ninfo, amin = "calc"):
        num_lines = len(scores)
        nplikeliTest = np.vectorize(likeliTest,otypes=[float])
        LikeLiHoods = nplikeliTest(ninfo, scores)
        if amin == "calc":
            TopHit = np.nanmin(LikeLiHoods)
        else:
            TopHit = float(amin)
        LikeLiHoodRatios = [get_fraction(LikeLiHoods[i], TopHit) for i in range(num_lines)]
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
        overlapScore = [get_fraction(self.ninfo[i], self.num_snps) for i in range(len(self.accs))]
        sorted_order = topHits[np.argsort(-self.probabilies[topHits])]
        (case, note) = self.case_interpreter(topHits)
        matches_dict = [(self.accs[i], self.probabilies[i], self.ninfo[i], overlapScore[i]) for i in sorted_order]
        topHitsDict = {'overlap': [self.overlap, self.num_snps], 'matches': matches_dict, 'interpretation':{'case': case, 'text': note}}
        with open(outFile, "w") as out_stats:
            out_stats.write(json.dumps(topHitsDict, sort_keys=True, indent=4))

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
        topHitsDict['percent_heterozygosity'] = get_fraction(numHets, len(snpGT))
        with open(outFile, "w") as out_stats:
          out_stats.write(json.dumps(topHitsDict, sort_keys=True, indent=4))
    return(get_fraction(numHets, len(snpGT)))

def genotyper(inputs, hdf5File, hdf5accFile, outFile):
    log.info("loading database files")
    g = snp_genotype.Genotype(hdf5File, hdf5accFile)
    log.info("done!")
    inputs.filter_chr_names()
    num_lines = len(g.g.accessions)
    ScoreList = np.zeros(num_lines, dtype="float")
    NumInfoSites = np.zeros(len(g.g.accessions), dtype="uint32")
    chunk_size = 1000
    commonSNPs = g.get_positions_idxs( inputs.chrs, inputs.pos )
    NumMatSNPs = len(commonSNPs[0])
    for j in range(0, NumMatSNPs, chunk_size):
        matchedAccInd = commonSNPs[0][j:j+chunk_size]
        matchedTarInd = commonSNPs[1][j:j+chunk_size]
        matchedTarWei = inputs.wei[matchedTarInd,]
        TarGTs0 = np.zeros(len(matchedTarInd), dtype="int8")
        TarGTs1 = np.ones(len(matchedTarInd), dtype="int8") + 1
        TarGTs2 = np.ones(len(matchedTarInd), dtype="int8")
        t1001SNPs = g.g.snps[matchedAccInd,:]
        samSNPs0 = np.reshape(np.repeat(TarGTs0, num_lines), (len(TarGTs0),num_lines))
        samSNPs1 = np.reshape(np.repeat(TarGTs1, num_lines), (len(TarGTs1),num_lines))
        samSNPs2 = np.reshape(np.repeat(TarGTs2, num_lines), (len(TarGTs2),num_lines))
        tempScore0 = np.sum(np.multiply(np.array(t1001SNPs == samSNPs0, dtype=int).T, matchedTarWei[:,0]).T, axis=0)
        tempScore1 = np.sum(np.multiply(np.array(t1001SNPs == samSNPs1, dtype=int).T, matchedTarWei[:,1]).T, axis=0)
        tempScore2 = np.sum(np.multiply(np.array(t1001SNPs == samSNPs2, dtype=int).T, matchedTarWei[:,2]).T, axis=0)
        ScoreList = ScoreList + tempScore0 + tempScore1 + tempScore2
        NumInfoSites = NumInfoSites + len(TarGTs0) - np.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
        if j % (chunk_size * 50) == 0:
            log.info("Done analysing %s positions", j+chunk_size)
    log.info("writing score file!")
    overlap = get_fraction(NumMatSNPs, len(inputs.pos))
    result = GenotyperOutput(g.g.accessions, ScoreList, NumInfoSites, overlap, NumMatSNPs, inputs.dp)
    result.print_out_table( outFile + '.scores.txt' )
    result.print_json_output( outFile + ".matches.json" )
    getHeterozygosity(inputs.gt[commonSNPs[1]], outFile + ".matches.json")
    return(result)

def potatoGenotyper(args):
    inputs = parsers.ParseInputs(inFile = args['inFile'], logDebug = args['logDebug'])
    log.info("running genotyper!")
    result = genotyper(inputs, args['hdf5File'], args['hdf5accFile'], args['outFile'])
    log.info("finished!")

def pairwiseScore(inFile_1, inFile_2, logDebug, outFile, hdf5File = None):
    snpmatch_stats = {}
    log.info("loading input files")
    inputs_1 = parsers.ParseInputs(inFile = inFile_1, logDebug = logDebug)
    inputs_2 = parsers.ParseInputs(inFile = inFile_2, logDebug = logDebug)
    if hdf5File is not None:
        log.info("loading database file to identify common SNP positions")
        g = snp_genotype.Genotype(hdf5File)
        snpmatch_stats['hdf5'] = hdf5File
        commonSNPs_1 = g.get_positions_idxs( inputs_1.chrs, inputs_1.pos )
        common_inds = snp_genotype.Genotype.get_common_positions( inputs_1.chrs[commonSNPs_1[1]], inputs_1.pos[commonSNPs_1[1]], inputs_2.chrs, inputs_2.pos )
        common_inds = (commonSNPs_1[1][common_inds[0]], common_inds[1])
    else:
        log.info("identify common positions")
        common_inds = snp_genotype.Genotype.get_common_positions( inputs_1.chrs, inputs_1.pos, inputs_2.chrs, inputs_2.pos )
    log.info("done!")
    unique_1 = len(inputs_1.chrs) - len(common_inds[0])
    unique_2 = len(inputs_2.chrs) - len(common_inds[0])
    common = np.zeros(0, dtype=int)
    scores = np.zeros(0, dtype=int)
    inputs_1.filter_chr_names()
    inputs_2.filter_chr_names()
    common_chrs = np.intersect1d(inputs_1.g_chrs_ids, inputs_2.g_chrs_ids)
    for i in common_chrs:
        perchrTarInd = np.where(inputs_1.g_chrs[common_inds[0]] == i)[0]
        log.info("Analysing chromosome %s positions", i)
        t_common = len(perchrTarInd)
        t_scores = np.sum(np.array(inputs_1.gt[common_inds[0][perchrTarInd]] == inputs_2.gt[common_inds[1][perchrTarInd]], dtype = int))
        snpmatch_stats[i] = [get_fraction(t_scores, t_common), t_common]
        common = np.append(common, t_common)
        scores = np.append(scores, t_scores)
    snpmatch_stats['matches'] = [get_fraction(np.sum(scores), np.sum(common)), np.sum(common)]
    snpmatch_stats['unique'] = {"%s" % os.path.basename(inFile_1): [get_fraction(unique_1, len(inputs_1.chrs)), len(inputs_1.chrs)], "%s" % os.path.basename(inFile_2): [get_fraction(unique_2, len(inputs_2.chrs)), len(inputs_2.chrs)] }
    if not outFile:
        outFile = "genotyper"
    log.info("writing output in a file: %s" % outFile + ".matches.json")
    with open(outFile + ".matches.json", "w") as out_stats:
        out_stats.write(json.dumps(snpmatch_stats, sort_keys=True, indent=4))
    log.info("finished!")
