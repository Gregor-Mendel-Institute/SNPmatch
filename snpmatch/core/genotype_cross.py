"""
  SNPmatch for crosses
"""
import numpy as np
import numpy.ma
import scipy.stats as st
from pygwas.core import genotype
import pandas as pd
import logging
import os
from . import snpmatch
from . import csmatch
import parsers
import json
import itertools

log = logging.getLogger(__name__)

# Arabidopsis chromosome length
chrlen = np.array((30427671, 19698289, 23459830, 18585056, 26975502))
tair_chrs = ['1', '2', '3', '4', '5']
mean_recomb_rates = [3.4, 3.6, 3.5, 3.8, 3.6]  ## cM/Mb ## Salome, P et al. 2011


def getWindowGenotype(matchedNos, totalMarkers, lr_thres, n_marker_thres = 5):
    ## matchedNos == array with matched number of SNPs
    ### Choose lr_thres as 2.706 which is at 0.1 alpha level with 1 degree of freedom
    pval = ''
    geno = 'NA'
    if totalMarkers < n_marker_thres:
        return((geno, 'NA'))
    assert len(matchedNos) == 3
    if np.array_equal(np.array(matchedNos), np.repeat(0, 3)):
        return((geno, 'NA'))
    likes = snpmatch.GenotyperOutput.calculate_likelihoods(matchedNos, np.repeat(totalMarkers, 3).tolist())
    for item in likes[1]:
        if pval == '':
            pval = pval + "%.2f" % item
        else:
            pval = pval + ',' + "%.2f" % item
    if len(np.where( likes[1] == 1 )[0]) > 1: ## It is matching to multiple
        return((1, pval))
    high_match = np.nanargmin(likes[0])
    lr_next = np.nanmin(likes[1][np.nonzero(likes[1]-1)])
    if np.isnan(lr_next):
        lr_next = lr_thres
    if np.nanargmin(likes[0]) == 0 and lr_next >= lr_thres:
        geno = 0
    elif np.nanargmin(likes[0]) == 2 and lr_next >= lr_thres:
        geno = 2
    if high_match == 1:
        geno = 1
    return(geno, pval)

## New class for genotype cross
class GenotypeCross(object):

    def __init__(self, hdf5_acc, parents, binLen, father = None, logDebug=True):
        self.logDebug = logDebug
        self.get_segregating_snps_parents(hdf5_acc, parents, father)
        self.window_size = int(binLen)

    def get_segregating_snps_parents(self, hdf5_acc, parents, father):
        log.info("loading genotype data for parents, and identify segregating SNPs")
        if father is not None:
            log.info("input files: %s and %s" % (parents, father))
            if not os.path.isfile(parents) and os.path.isfile(father):
                die("either of the input files do not exists, please provide VCF/BED file for parent genotype information")
            p1_snps = parsers.ParseInputs(inFile = parents, logDebug = self.logDebug)
            p2_snps = parsers.ParseInputs(inFile = father, logDebug = self.logDebug)
            commonCHRs_ids = np.union1d(p1_snps.chrs, p2_snps.chrs)
            commonSNPsCHR = np.zeros(0, dtype=commonCHRs_ids.dtype)
            commonSNPsPOS = np.zeros(0, dtype=int)
            snpsP1 = np.zeros(0, dtype='int8')
            snpsP2 = np.zeros(0, dtype='int8')
            for i in commonCHRs_ids:
                perchrP1inds = np.where(p1_snps.chrs == i)[0]
                perchrP2inds = np.where(p2_snps.chrs == i)[0]
                perchrPositions = np.union1d(p1_snps.pos[perchrP1inds], p2_snps.pos[perchrP2inds])
                commonSNPsCHR = np.append(commonSNPsCHR, np.repeat(i, len(perchrPositions)))
                commonSNPsPOS = np.append(commonSNPsPOS, perchrPositions)
                perchrsnpsP1_inds = np.where(np.in1d(p1_snps.pos[perchrP1inds], perchrPositions))[0]
                perchrsnpsP2_inds = np.where(np.in1d(p2_snps.pos[perchrP2inds], perchrPositions))[0]
                snpsP1 = np.append(snpsP1, parsers.parseGT(p1_snps.gt[perchrsnpsP1_inds]))
                snpsP2 = np.append(snpsP2, parsers.parseGT(p2_snps.gt[perchrsnpsP2_inds]))
            log.info("done!")
        else:
            ## need to filter the SNPs present in C and M
            log.info("loading HDF5 file")
            g_acc = genotype.load_hdf5_genotype_data(hdf5_acc)
            ## die if either parents are not in the dataset
            assert len(parents.split("x")) == 2, "parents should be provided as '6091x6191'"
            try:
                indP1 = np.where(g_acc.accessions == parents.split("x")[0])[0][0]
                indP2 = np.where(g_acc.accessions == parents.split("x")[1])[0][0]
            except:
                snpmatch.die("parents are not in the dataset")
            snpsP1 = g_acc.snps[:,indP1]
            snpsP2 = g_acc.snps[:,indP2]
            commonSNPsCHR = np.array(g_acc.chromosomes)
            commonSNPsPOS = np.array(g_acc.positions)
            log.info("done!")
        segSNPsind = np.where((snpsP1 != snpsP2) & (snpsP1 >= 0) & (snpsP2 >= 0) & (snpsP1 < 2) & (snpsP2 < 2))[0]
        log.info("number of segregating snps between parents: %s", len(segSNPsind))
        self.commonSNPsCHR = commonSNPsCHR[segSNPsind]
        self.commonSNPsPOS = commonSNPsPOS[segSNPsind]
        self.snpsP1 = snpsP1[segSNPsind]
        self.snpsP2 = snpsP2[segSNPsind]
        log.info("done!")

    def _get_common_positions_ixs(self, accs_ix, sample_ix):
        ## Make sure you have inputs loaded
        reqAccsPOS = self.commonSNPsPOS[accs_ix]
        reqTarPos = self._inputs.pos[sample_ix]
        commonPOS = np.intersect1d( reqAccsPOS, reqTarPos )
        matchedAccInd = np.array(accs_ix)[ np.where( np.in1d(reqAccsPOS, reqTarPos) )[0] ]
        matchedTarInd = np.array(sample_ix)[ np.where( np.in1d(reqTarPos, reqAccsPOS) )[0] ]
        return((matchedAccInd, matchedTarInd))

    @staticmethod
    def get_window_genotype_gts(input_gt, snpsP1_gt, snpsP2_gt, lr_thres):
        # input_gt is only '0/0', '0/1', '1/1'
        # snpsP1_gt and snpsP2_gt is either 0, 1 or 2
        num_snps = len(input_gt)
        assert num_snps == len(snpsP1_gt), "provide same number of SNPs"
        assert num_snps == len(snpsP2_gt), "provide same number of SNPs"
        TarGTBinary = parsers.parseGT(input_gt)
        matP1no = len(np.where(np.equal( TarGTBinary, snpsP1_gt ))[0])
        matP2no = len(np.where(np.equal( TarGTBinary, snpsP2_gt ))[0])
        matHetno = len(np.where(np.equal( TarGTBinary, np.repeat(2, num_snps) ))[0])
        return(getWindowGenotype([matP1no, matHetno, matP2no], num_snps, lr_thres))

    def genotype_each_cross(self, input_file, lr_thres):
        ## Input file
        self._inputs = parsers.ParseInputs(inFile = input_file, logDebug = self.logDebug)
        ## Inputs is the ParseInputs class object
        log.info("running cross genotyper")
        iter_bins_genome = csmatch.get_bins_arrays(self.commonSNPsCHR, self.commonSNPsPOS, self.window_size)
        iter_bins_snps = csmatch.get_bins_arrays(self._inputs.chrs, self._inputs.pos, self.window_size)
        bin_inds = 0
        outfile_str = np.zeros(0, dtype="string")
        for e_b, e_s in itertools.izip(iter_bins_genome, iter_bins_snps):
            # first snp positions which are segregating and are in this window
            bin_str = tair_chrs[e_b[0]] + "\t" + str(e_b[1][0]) + "\t" + str(e_b[1][1])
            matchedAccInd, matchedTarInd = self._get_common_positions_ixs( e_b[2], e_s[2] )
            if len(matchedTarInd) == 0:
                outfile_str = np.append(outfile_str, "%s\t%s\t%s\tNA\tNA" % (bin_str, len(matchedTarInd), len(e_b[2])) )
            else:
                matchedTarGTs = self._inputs.gt[matchedTarInd]
                (geno, pval) = self.get_window_genotype_gts(matchedTarGTs, self.snpsP1[matchedAccInd], self.snpsP2[matchedAccInd], lr_thres)
                outfile_str = np.append(outfile_str, "%s\t%s\t%s\t%s\t%s" % (bin_str, len(matchedTarInd), len(e_b[2]), geno, pval))
            bin_inds += 1
            if bin_inds % 40 == 0:
                log.info("progress: %s windows", bin_inds)
        log.info("done!")
        return(outfile_str)

    @staticmethod
    def get_probabilities_ebin(input_gt, snpsP1_gt, snpsP2_gt, error_prob = 0.01):
        num_snps = len(input_gt)
        ebTarGTs = parsers.parseGT(input_gt)
        ebPolarised = np.zeros(num_snps, dtype=int)
        ebPolarised[np.where(np.equal( ebTarGTs, snpsP1_gt ))[0] ] = 0
        ebPolarised[np.where(np.equal( ebTarGTs, snpsP2_gt ))[0] ] = 2
        ebPolarised[np.where(np.equal( ebTarGTs, np.repeat(2, num_snps) ))[0] ] = 1
        eb_probs = np.repeat(error_prob / 2, num_snps * 3).reshape((-3, 3))
        for i_ix in range(num_snps):
            eb_probs[i_ix,ebPolarised[i_ix]] = 1 - error_prob
        return( eb_probs )


    def genotype_each_cross_hmm(self, input_file, n_marker_thres):
        self._inputs = parsers.ParseInputs(inFile = input_file, logDebug = self.logDebug)
        ## Inputs is the ParseInputs class object
        self._inputs.get_chr_list()
        from hmmlearn import hmm
        ### The below transition probability is for a intercross, adapted from R/qtl
        log.info("running HMM")
        #r_i = float(self.window_size) * 0.38 /1000000
        r_i = 4000 * 3.8 /1000000
        e_p = 0.01
        assert r_i < 1, "Provide either small window size or check the mean recombination rate"
        log.info("recombination fraction between windows: %.2f" % r_i)
        transprob = np.array( [[ (1-r_i)**2, 2*r_i*(1-r_i), r_i**2 ], [r_i*(1-r_i), (1-r_i)**2 + r_i**2, r_i*(1-r_i)], [ r_i**2, 2*r_i*(1-r_i), (1-r_i)**2 ]] )
        log.info("error rate for genotyping: %.2f" % e_p)
        model = hmm.GaussianHMM(n_components=3, covariance_type="full", n_iter = 1000)
        model.startprob_ = np.array([0.25, 0.5, 0.25])
        model.transmat_ = transprob
        log.info("Transision matrix: %s" % transprob)
        model.means_ = np.array( [ [1, 0, 0], [0, 1, 0], [0, 0, 1] ] )
        model.covars_ = np.tile(np.identity(3), (3, 1, 1))
        allSNPGenos = np.zeros(0, dtype=int)
        for ec, eclen in itertools.izip(tair_chrs, chrlen):
            reqPOSind =  np.where(self.commonSNPsCHR == ec)[0]
            reqTARind = np.where( self._inputs.g_chrs == ec )[0]
            ebAccsInds, ebTarInds = self._get_common_positions_ixs(reqPOSind, reqTARind)
            eb_SNPprobs = np.zeros((len(reqPOSind), 3))
            ebin_prob = self.get_probabilities_ebin(self._inputs.gt[ebTarInds], self.snpsP1[ebAccsInds], self.snpsP2[ebAccsInds])
            eb_SNPprobs[ebAccsInds - reqPOSind[0],] = ebin_prob
            allSNPGenos = np.append(allSNPGenos, model.predict( eb_SNPprobs ))
        iter_bins_genome = csmatch.get_bins_arrays(self.commonSNPsCHR, self.commonSNPsPOS, self.window_size)
        iter_bins_snps = csmatch.get_bins_arrays(self._inputs.chrs, self._inputs.pos, self.window_size)
        outfile_str = np.zeros(0, dtype="string")
        for e_b, e_s in itertools.izip(iter_bins_genome, iter_bins_snps):
            matchedAccInd, matchedTarInd = self._get_common_positions_ixs( e_b[2], e_s[2] )
            bin_str = tair_chrs[e_b[0]] + "\t" + str(e_b[1][0]) + "\t" + str(e_b[1][1]) + "\t" + str(len(matchedTarInd)) + "\t" + str(len(e_b[2]))
            t_genos = np.unique(allSNPGenos[e_b[2]])
            t_genos_nums = pd.Series([ len( np.where( allSNPGenos[e_b[2]] == e )[0] )  for e in [0, 1, 2]]).astype(str).str.cat(sep=",")
            if len(e_b[2]) > 0:
                if len(t_genos) == 1:
                    outfile_str = np.append(outfile_str, "%s\t%s\t%s" % (bin_str, t_genos[0], t_genos_nums) )
                else:
                    outfile_str = np.append(outfile_str, "%s\t%s\t%s" % (bin_str, 1, t_genos_nums) )
            else:
                outfile_str = np.append(outfile_str, "%s\t%s\t%s" % (bin_str, 'NA', 'NA' ) )
        return(outfile_str)


    @staticmethod
    def write_output_genotype_cross(outfile_str, output_file ):
        log.info("writing file: %s" % output_file)
        outfile = open(output_file, 'w')
        for ef in outfile_str:
            outfile.write( "%s\n" % ef )
        outfile.close()
        log.info("done!")

    def filter_good_samples(self, snpvcf, good_samples_file):
        if good_samples_file is None:
            return(snpvcf)
        good_samples = np.array(pd.read_table(good_samples_file, header = None), dtype="string")
        good_samples_ix = np.zeros(0, dtype=int)
        for ef_ix,ef in enumerate(snpvcf.columns.values[2:]):
            find_ix = np.where(good_samples == ef.split("_")[0])[0]
            if len(find_ix) == 0:
                find_ix = np.where(good_samples == ef.split("_")[0] + ef.split("_")[1] )[0]
            if len(find_ix) > 0:
                good_samples_ix = np.append(good_samples_ix, ef_ix + 2)
        return(snpvcf.iloc[:, np.append((0,1), good_samples_ix) ])

    @staticmethod
    def print_ids(snpvcf):
        all_samples = pd.Series(snpvcf.columns.values[2:]).str.split("_", n = 2, expand=True)
        if len(np.unique(all_samples.iloc[:,0])) == all_samples.shape[0]:
            return(all_samples.iloc[:,0].str.cat(sep = ','))
        all_samples = all_samples.iloc[:,0].map(str) + all_samples.iloc[:,1].map(str)
        return(all_samples.str.cat(sep = ','))

    def genotype_cross_all_samples(self, sample_file, lr_thres, good_samples_file=None):
        log.info("loading input files!")
        snpvcf = pd.read_table(sample_file)
        snpvcf = self.filter_good_samples(snpvcf, good_samples_file)
        num_samples = snpvcf.shape[1] - 2
        log.info("number of samples printed: %s" % num_samples )
        iter_bins_genome = get_bins_arrays(self.commonSNPsCHR, self.commonSNPsPOS, self.window_size)
        iter_bins_snps = get_bins_arrays(np.array(snpvcf.iloc[:,0]), np.array(snpvcf.iloc[:,1]), self.window_size)
        bin_inds = 0
        outfile_str = np.array(('id,,,' + self.print_ids(snpvcf)), dtype="string")
        outfile_str = np.append(outfile_str, 'pheno,' + ',' + ',0' * num_samples)
        for e_b, e_s in itertools.izip(iter_bins_genome, iter_bins_snps):
            bin_str = tair_chrs[e_b[0]] + ":" + str(e_b[1][0]) + "-" + str(e_b[1][1])
            cm_mid = float(mean_recomb_rates[e_b[0]]) * np.mean(e_b[1]).astype(int) / 1000000
            reqPOS = self.commonSNPsPOS[e_b[2]]
            perchrTarPos = np.array(snpvcf.iloc[e_s[2], 1])
            matchedAccInd = np.array(e_b[2], dtype=int)[ np.where( np.in1d(reqPOS, perchrTarPos) )[0] ]
            matchedTarInd = np.array(e_s[2], dtype=int)[ np.where( np.in1d(perchrTarPos, reqPOS) )[0] ]
            if len(matchedTarInd) == 0:
                outfile_str = np.append(outfile_str, "%s,%s,%s%s" % (bin_str, tair_chrs[e_b[0]],  cm_mid, ',NA' * num_samples ) )
            else:
                geno_samples = ''
                for sample_ix in range(num_samples):
                    (geno, pval) = self.get_window_genotype_gts(np.array(snpvcf.iloc[matchedTarInd, sample_ix + 2]), self.snpsP1[matchedAccInd], self.snpsP2[matchedAccInd], lr_thres)
                    geno_samples = geno_samples + ',' + str(geno)
                outfile_str = np.append(outfile_str, "%s,%s,%s%s" % (bin_str, tair_chrs[e_b[0]],  cm_mid, geno_samples ) )
            bin_inds += 1
            if bin_inds % 40 == 0:
                log.info("progress: %s windows", bin_inds)
        log.info("done!")
        return(outfile_str)


def potatoCrossGenotyper(args):
    ## Get the VCF file (filtered may be) generated by GATK.
    ## inputs:
    # 1) VCF file
    # 2) Parent1 and Parent2
    # 3) SNP matrix (hdf5 file)
    # 4) Bin length, default as 200Kbp
    # 5) Chromosome length
    crossgenotyper = GenotypeCross(args['hdf5accFile'], args['parents'], args['binLen'], args['father'], args['logDebug'])
    if args['all_samples']:
        outfile_str = crossgenotyper.genotype_cross_all_samples( args['inFile'], args['lr_thres'], args['good_samples'] )
    elif args['hmm']:
        outfile_str = crossgenotyper.genotype_each_cross_hmm( args['inFile'], 5 )
    else:
        outfile_str = crossgenotyper.genotype_each_cross( args['inFile'], args['lr_thres'] )
    crossgenotyper.write_output_genotype_cross( outfile_str, args['outFile'] )
