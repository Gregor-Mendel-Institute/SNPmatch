"""
  SNPmatch
"""
import numpy as np
import pandas as pd
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

def get_fraction(x, y, y_min = 0):
    if y <= y_min:
        return(np.nan)
    return(float(x)/y)

np_get_fraction = np.vectorize(get_fraction, excluded = "y_min")
# def np_get_fraction(x, y, y_min = 0):
#     x = np.array(x)
#     y = np.array(y)
#     with np.errstate(divide='ignore', invalid='ignore'):
#         p = np.divide( x, y )
#         p[np.where(y <= y_min)] = np.nan
#     return(p)


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

def np_binom_test(x, n, p, alternative=None):
    """
    Code adapted from statsmodels.proportion.binom_test
    """
    if alternative in ['larger','greater']:
        pval = stats.binom.sf(x-1, n, p)
    elif alternative in ['smaller','less']:
        pval = stats.binom.cdf(x, n, p)
    else:
        np_binom = np.vectorize( stats.binom_test, excluded=['p', 'alternative'])
        pval = np_binom( x, n, p, alternative )
    return(pval)

def np_test_identity(x, n, error_rate = 0.0005, pthres = 0.05):
    st = np_binom_test(n - x, n, p = error_rate, alternative="greater" )
    return(np.array(st >= pthres).astype(int))

def matchGTsAccs(sampleWei, t1001snps, skip_hets_db = False):
    assert sampleWei.shape[0] == t1001snps.shape[0], "please provide same number of positions for both sample and db"
    assert sampleWei.shape[1] == 3, "SNP weights should be a np.array with  shape == n,3"
    ## Initilizing
    if skip_hets_db: 
        t1001snps[t1001snps == 2] = -1
    num_lines = t1001snps.shape[1]
    TarGTs0 = np.zeros(t1001snps.shape, dtype="int8") ## Homo -- ref
    TarGTs1 = np.ones(t1001snps.shape, dtype="int8") + 1  ## Hets 
    TarGTs2 = np.ones(t1001snps.shape, dtype="int8")  ## Homo -- alt
    score = np.zeros( num_lines )
    score = score + np.multiply(np.array(numpy.ma.masked_less(t1001snps, 0) == TarGTs0, dtype=int).T, sampleWei[:,0]).sum(axis = 1)
    score = score + np.multiply(np.array(numpy.ma.masked_less(t1001snps, 0) == TarGTs1, dtype=int).T, sampleWei[:,1]).sum(axis = 1)
    score = score + np.multiply(np.array(numpy.ma.masked_less(t1001snps, 0) == TarGTs2, dtype=int).T, sampleWei[:,2]).sum(axis = 1)
    ninfo = np.repeat(t1001snps.shape[0], num_lines) - np.sum(numpy.ma.masked_less(t1001snps, 0).mask.astype(int ), axis = 0)
    return((score, ninfo))

class GenotyperOutput(object):
    ## class object for main SNPmatch output

    def __init__(self, AccList, ScoreList, NumInfoSites, overlap, NumMatSNPs, DPmean ):
        self.accs = np.array(AccList, dtype="str")
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
        output_table = pd.DataFrame( {
            'accs': self.accs,
            'matches': self.scores,
            'ninfo': self.ninfo,
            'probabilities': self.probabilies,
            'likelihood': self.likelis,
            'lrt': self.lrts,
            'num_snps': self.num_snps,
            'dp': np.nanmean(self.dp)
        } )
        output_table = output_table[ ['accs', 'matches', 'ninfo', 'probabilities', 'likelihood', 'lrt', 'num_snps', 'dp'] ]
        if outFile:
            output_table.to_csv( outFile, header = None, sep = "\t", index = None )
        return( output_table )

    def print_json_output(self, outFile):
        self.get_likelihoods()
        self.get_probabilities()
        topHits = np.where(self.lrts < lr_thres)[0]
        overlapScore = [get_fraction(self.ninfo[i], self.num_snps) for i in range(len(self.accs))]
        sorted_order = topHits[np.argsort(-self.probabilies[topHits])]
        (case, note) = self.case_interpreter(topHits)
        matches_dict = [(str(self.accs[i]), float(self.probabilies[i]), int(self.ninfo[i]), float(overlapScore[i])) for i in sorted_order]
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
            note = "Ambiguous sample: Many input SNP positions are missing in db positions. Maybe sample  not one in database"
        return(case, note)

class Genotyper(object):
    ## class object for main SNPmatch

    def __init__(self, inputs, g, outFile, run_genotyper = True, skip_db_hets = False, chunk_size = 1000):
        assert type(g) is snp_genotype.Genotype, "provide a snp_genotype.Genotype class for genotypes"
        inputs.filter_chr_names()
        self.chunk_size = chunk_size
        self.inputs = inputs
        self.g = g
        self.num_lines = len(self.g.g.accessions)
        self.outFile = outFile
        self._skip_db_hets = skip_db_hets
        if run_genotyper:
            self.result = self.genotyper()
            self.write_genotyper_output( self.result )

    def get_common_positions(self):
        self.commonSNPs = self.g.get_positions_idxs( self.inputs.chrs, self.inputs.pos )

    def filter_tophits(self):
        self.result = self.genotyper()
        self.write_genotyper_output( self.result )
        self.result.get_likelihoods()
        topHits =  np.where(self.result.lrts < lr_thres)[0]
        if len(topHits) == 1:
            log.info("Done! It is a perfect hit")
            return(None)
        log.info("#lines indistinguishable: %s" % len(topHits))
        log.info("refining likelihoods for only indistinguishable lines")
        if len(topHits) > (self.num_lines / 2):
            log.info("too many lines are indistinguishable, skipping refining likelihoods step")
            return(None)
        seg_ix = self.g.identify_segregating_snps( topHits )
        self.result_fine = self.genotyper( filter_pos_ix = seg_ix, mask_acc_ix = np.where(self.result.lrts >= lr_thres)[0]  )
        log.info("writing output: %s" % self.outFile + ".refined.scores.txt")
        self.result_fine.print_out_table( self.outFile + ".refined.scores.txt" )

    def genotyper(self, filter_pos_ix = None, mask_acc_ix = None):
        ScoreList = np.zeros(self.num_lines, dtype="float")
        NumInfoSites = np.zeros(len(self.g.g.accessions), dtype="uint32")
        self.get_common_positions()
        if filter_pos_ix is not None:
            assert type(filter_pos_ix) is np.ndarray, "provide np array for indices to be considered"
            t_ix = np.where(np.in1d(self.commonSNPs[0], filter_pos_ix))[0]
            if t_ix.shape[0] < 100:
                log.info("#positions in segregating sites are are too little: %s" % t_ix.shape[0])
            self.commonSNPs = (self.commonSNPs[0][t_ix], self.commonSNPs[1][t_ix] )
        NumMatSNPs = len(self.commonSNPs[0])
        for j in range(0, NumMatSNPs, self.chunk_size):
            matchedAccInd = self.commonSNPs[0][j:j+self.chunk_size]
            matchedTarInd = self.commonSNPs[1][j:j+self.chunk_size]
            matchedTarWei = self.inputs.wei[matchedTarInd,]
            t1001SNPs = self.g.g.snps[matchedAccInd,:]
            t_s, t_n = matchGTsAccs( matchedTarWei, t1001SNPs, self._skip_db_hets )
            ScoreList = ScoreList + t_s
            NumInfoSites = NumInfoSites + t_n
            if j % (self.chunk_size * 50) == 0:
                log.info("Done analysing %s positions", j+self.chunk_size)
        overlap = get_fraction(NumMatSNPs, len(self.inputs.pos))
        if mask_acc_ix is not None:
            assert type(mask_acc_ix) is np.ndarray, "provide a numpy array of accessions indices to mask"
            mask_acc_to_print = np.setdiff1d(np.arange( self.num_lines ), mask_acc_ix)
            return( GenotyperOutput(self.g.g.accessions[mask_acc_to_print], ScoreList[mask_acc_to_print], NumInfoSites[mask_acc_to_print], overlap, NumMatSNPs, self.inputs.dp) )
        return( GenotyperOutput(self.g.g.accessions, ScoreList, NumInfoSites, overlap, NumMatSNPs, self.inputs.dp) )

    def write_genotyper_output(self, result):
        log.info("writing score file!")
        result.get_likelihoods()
        result.print_out_table( self.outFile + '.scores.txt' )
        result.print_json_output( self.outFile + ".matches.json" )
        getHeterozygosity(self.inputs.gt[self.commonSNPs[1]], self.outFile + ".matches.json")
        return(result)


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


def potatoGenotyper(args):
    inputs = parsers.ParseInputs(inFile = args['inFile'], logDebug = args['logDebug'])
    log.info("loading database files")
    g = snp_genotype.Genotype(args['hdf5File'], args['hdf5accFile'])
    log.info("done!")
    log.info("running genotyper!")
    if args['refine']:
        genotyper = Genotyper(inputs, g, args['outFile'], run_genotyper=False,  skip_db_hets = args['skip_db_hets'])
        genotyper.filter_tophits()
        log.info("finished!")
        return(None)
    genotyper = Genotyper(inputs, g, args['outFile'], run_genotyper=True,  skip_db_hets = args['skip_db_hets'])
    log.info("finished!")

def pairwiseScore(inFile_1, inFile_2, logDebug, outFile = None, hdf5File = None):
    snpmatch_stats = {}
    log.info("loading input files")
    inputs_1 = parsers.ParseInputs(inFile = inFile_1, logDebug = logDebug)
    inputs_2 = parsers.ParseInputs(inFile = inFile_2, logDebug = logDebug)
    if hdf5File is not None:
        log.info("loading database file to identify common SNP positions")
        g = snp_genotype.Genotype(hdf5File, None)
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
    if outFile:
        # outFile = "genotyper"
        log.info("writing output in a file: %s" % outFile + ".matches.json")
        with open(outFile + ".matches.json", "w") as out_stats:
            out_stats.write(json.dumps(snpmatch_stats, sort_keys=True, indent=4))
        log.info("finished!")
    return(snpmatch_stats)
