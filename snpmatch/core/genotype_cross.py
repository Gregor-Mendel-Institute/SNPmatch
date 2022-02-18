"""
  SNPmatch for crosses
"""
import numpy as np
import numpy.ma
import scipy.stats as st
import pandas as pd
import logging
import os
from . import snpmatch
from . import snp_genotype
from . import infer
from . import genomes
from . import parsers
import json
import itertools

log = logging.getLogger(__name__)


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

    def __init__(self, g, parents, binLen = 0, father = None, logDebug=True):
        self.logDebug = logDebug
        self.g = g
        self.get_segregating_snps_parents(parents, father)
        self.window_size = int(binLen)

    def get_segregating_snps_parents(self, parents, father):
        log.info("loading genotype data for parents, and identify segregating SNPs")
        if father is not None:
            log.info("input files: %s and %s" % (parents, father))
            if not os.path.isfile(parents) and os.path.isfile(father):
                snpmatch.die("either of the input files do not exists, please provide VCF/BED file for parent genotype information")
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
            ## die if either parents are not in the dataset
            assert len(parents.split("x")) == 2, "parents should be provided as '6091x6191'"
            try:
                indP1 = np.where(self.g.accessions == parents.split("x")[0])[0][0]
                indP2 = np.where(self.g.accessions == parents.split("x")[1])[0][0]
            except:
                snpmatch.die("parents are not in the dataset")
            snpsP1 = self.g.g_acc.snps[:,indP1]
            snpsP2 = self.g.g_acc.snps[:,indP2]
            self.p1_ix = indP1
            self.p2_ix = indP2
            commonSNPsCHR = np.array(self.g.g_acc.chromosomes)
            commonSNPsPOS = np.array(self.g.g_acc.positions)
            log.info("done!")
        ## only considering sites where two parents differ ---
        # 1. HMM can now incorporate parental heterozygous sites. would change the emissions  
        # 2. also excluding sites when it is NA in atleast one parent.
        ## you do not need a lot of markers to construct the map .. but rather perfect sites
        segSNPsind = np.where((snpsP1 != snpsP2) & (snpsP1 >= 0) & (snpsP2 >= 0) )[0]
        # segSNPsind = np.where((snpsP1 != snpsP2) )[0]
        log.info("number of segregating snps between parents: %s", len(segSNPsind))
        self.commonSNPsCHR = commonSNPsCHR[segSNPsind].astype('U')
        self.commonSNPsPOS = commonSNPsPOS[segSNPsind]
        self.snpsP1 = snpsP1[segSNPsind]
        self.snpsP2 = snpsP2[segSNPsind]
        log.info("done!")

    def genotype_cross_hmm(self, input_file, min_na_per_sample = 0.8):
        """
        Function to genotype an F2 individual based on segregating SNPs using a simple HMM model
        """
        snpvcf = parsers.import_vcf_file(
            inFile = input_file, 
            logDebug = self.logDebug, 
            samples_to_load = None, 
            add_fields = ['calldata/DP'] 
        )
        samples_ids = pd.Series(snpvcf['samples']) #.str.replace("_processed_reads_no_clonal.bam", "" )
        # if np.unique(samples_ids.str.split("_", expand = True).iloc[:,0]).shape[0] == samples_ids.shape[0]:
        #     samples_ids = samples_ids.str.split("_", expand = True).iloc[:,0]
        samples_gt = snpvcf['gt']
        samples_gt = pd.DataFrame(samples_gt.astype(str))
        g_chr_names = genome.chrs[pd.Series(self.commonSNPsCHR, dtype = str).apply(genome.get_chr_ind)]
        ## identify positions which are segregating between parents
        segregating_ix = self.g.get_common_positions( self.commonSNPsCHR, self.commonSNPsPOS, snpvcf['chr'], snpvcf['pos'] )
        num_markers = segregating_ix[1].shape[0]
        samples_gt = samples_gt.iloc[segregating_ix[1],:]
        samples_dp = snpvcf['calldata/DP'][segregating_ix[1],:] #
        filter_lowcov_ix = (samples_dp <= 0).sum(axis = 0) / float(num_markers)
        filter_lowcov_ix = np.where( filter_lowcov_ix < min_na_per_sample )[0]
        log.info("filtering %s samples due to very low number of informative markers" % str(samples_ids.shape[0] - filter_lowcov_ix.shape[0] ) )
        samples_gt = samples_gt.iloc[:,filter_lowcov_ix]
        samples_dp = samples_dp[:,filter_lowcov_ix] / 2
        ## Filter SNP positions where you have "missing" reads -- Vetmed et al. 
        #1. bcftools -e 'INFO/RO + sum(INFO/AO) < 0.9 * INFO/DP'
        # reference alleles + alternative alleles < 0.9 of read depth
        # 2. Uneven mean mapping quality of ALT and REF
        # bcftools -e 'MQM/MQMR < 0.9'

        samples_ids = samples_ids.iloc[ filter_lowcov_ix ]
        allSNPGenos_raw = pd.DataFrame( 
            index = pd.Series(self.commonSNPsCHR[segregating_ix[0]]).astype(str) + ":" + pd.Series(self.commonSNPsPOS[segregating_ix[0]]).astype(str),
            columns=samples_gt.columns
        )
        allSNPGenos = pd.DataFrame( index = allSNPGenos_raw.index, columns = allSNPGenos_raw.columns )
        if "recomb_rates" in genome.json.keys():
            mean_recomb_rates = np.mean(np.array(genome.json['recomb_rates']))
        else:
            log.warn("Average recombination rates were missing in genome file. Add rates for each chromosome as an array in genome json file under 'recomb_rates' key. Using default rate of 3.5")
            mean_recomb_rates = 3.5
        
        for ec, eclen in zip(genome.chrs_ids, genome.chrlen):
            reqChrind = np.where( g_chr_names[segregating_ix[0]] == ec )[0]
            for sample_ix in range(samples_gt.shape[1]):
                t_sample_dp = samples_dp[reqChrind,sample_ix]
                
                t_model = infer.IdentifyAncestryF2individual(
                    chromosome_size= eclen/1000000, 
                    snps_p1 = self.snpsP1[segregating_ix[0][reqChrind]],
                    snps_p2 = self.snpsP2[segregating_ix[0][reqChrind]], 
                    recomb_rate = mean_recomb_rates, 
                    base_error = 0.036,
                    sample_depth= t_sample_dp
                )
                t_sample_snps = parsers.parseGT(samples_gt.iloc[reqChrind,sample_ix].values)
                allSNPGenos_raw.iloc[reqChrind,sample_ix] = infer.polarize_snps(t_sample_snps, t_model.params['snps_p1'], t_model.params['snps_p2'] )
                allSNPGenos.iloc[reqChrind,sample_ix] = np.array(t_model.viterbi( t_sample_snps )[0], dtype = int)
        pos_em_cm = pd.Series(allSNPGenos.index, index = allSNPGenos.index).str.replace(":",",").apply(genome.estimated_cM_distance)
        allSNPGenos = allSNPGenos.astype(str).agg(','.join, axis=1)
        to_rqtl_result = pd.Series(allSNPGenos.index, index = allSNPGenos.index )
        to_rqtl_result = to_rqtl_result + "," + pd.Series(allSNPGenos.index, index = allSNPGenos.index ).str.split(":", expand = True).iloc[:,0]
        to_rqtl_result = to_rqtl_result + "," + pos_em_cm.astype(str)
        to_rqtl_result = to_rqtl_result + "," + allSNPGenos
        to_rqtl_result = pd.Series(('pheno,,' + ',0' * len(samples_ids))).append(to_rqtl_result , ignore_index=True)
        to_rqtl_result = pd.Series( ('id,,,' + str(pd.Series(samples_ids).str.cat(sep=",")) ) ).append(to_rqtl_result , ignore_index=True)
        return(to_rqtl_result)
    

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

    def filter_good_samples(self, snpvcf, good_samples_file):
        if good_samples_file is None:
            return(snpvcf)
        good_samples = np.array(pd.read_table(good_samples_file, header = None), dtype=str)
        good_samples_ix = np.zeros(0, dtype=int)
        for ef_ix,ef in enumerate(snpvcf.columns.values[2:]):
            find_ix = np.where(good_samples == ef.split("_")[0])[0]
            if len(find_ix) == 0:
                find_ix = np.where(good_samples == ef.split("_")[0] + ef.split("_")[1] )[0]
            if len(find_ix) > 0:
                good_samples_ix = np.append(good_samples_ix, ef_ix + 2)
        return(snpvcf.iloc[:, np.append((0,1), good_samples_ix) ])

    def genotype_cross(self, input_file, lr_thres, good_samples_file=None):
        log.info("loading input files!")
        snpvcf = parsers.import_vcf_file(inFile = input_file, logDebug = self.logDebug, samples_to_load = None)
        num_samples = snpvcf['samples'].shape[0] 
        log.info("number of samples printed: %s" % num_samples )
        log.warn("Using an average recombination rates of 3. Please change it according or use R/qtl package to generate genetic map.")
        mean_recomb_rates = np.repeat(3, len(genome.chrs_ids))
        iter_bins_genome = genome.get_bins_arrays(self.commonSNPsCHR, self.commonSNPsPOS, self.window_size)
        iter_bins_snps = genome.get_bins_arrays(snpvcf['chr'], snpvcf['pos'], self.window_size)
        bin_inds = 0
        outfile_str = np.array(('id,,,' + pd.Series(snpvcf['samples']).str.cat(sep = ',')  ), dtype=str)
        outfile_str = np.append(outfile_str, 'pheno,' + ',' + ',0' * num_samples)
        for e_b, e_s in zip(iter_bins_genome, iter_bins_snps):
            bin_str = genome.chrs_ids[e_b[0]] + ":" + str(e_b[1][0]) + "-" + str(e_b[1][1])
            cm_mid = genome.estimated_cM_distance( genome.chrs_ids[e_b[0]] + "," + str(int(round(np.mean(e_b[1]))))  )
            reqPOS = self.commonSNPsPOS[e_b[2]]
            perchrTarPos = snpvcf['pos'][e_s[2]]
            matchedAccInd = np.array(e_b[2], dtype=int)[ np.where( np.in1d(reqPOS, perchrTarPos) )[0] ]
            matchedTarInd = np.array(e_s[2], dtype=int)[ np.where( np.in1d(perchrTarPos, reqPOS) )[0] ]
            if len(matchedTarInd) == 0:
                outfile_str = np.append(outfile_str, "%s,%s,%s%s" % (bin_str, genome.chrs_ids[e_b[0]],  cm_mid, ',NA' * num_samples ) )
            else:
                geno_samples = ''
                for sample_ix in range(num_samples):
                    (geno, pval) = self.get_window_genotype_gts(snpvcf['gt'][matchedTarInd,sample_ix], self.snpsP1[matchedAccInd], self.snpsP2[matchedAccInd], lr_thres)
                    geno_samples = geno_samples + ',' + str(geno)
                outfile_str = np.append(outfile_str, "%s,%s,%s%s" % (bin_str, genome.chrs_ids[e_b[0]],  cm_mid, geno_samples ) )
            bin_inds += 1
            if bin_inds % 40 == 0:
                log.info("progress: %s windows", bin_inds)
        log.info("done!")
        return(outfile_str)

    @staticmethod
    def write_output_genotype_cross(outfile_str, output_file ):
        log.info("writing file: %s" % output_file)
        outfile = open(output_file, 'w')
        for ef in outfile_str:
            outfile.write( "%s\n" % ef )
        outfile.close()
        log.info("done!")


def uniq_neighbor(a):
    sorted_a = np.array(a[0], dtype = a.dtype)
    sorted_a_count = np.array([1], dtype = int)
    for ef_ix in range(1, len(a)):
        if a[ef_ix] != a[ef_ix-1]:
            sorted_a = np.append(sorted_a, a[ef_ix])
            sorted_a_count = np.append(sorted_a_count, 1)
        elif a[ef_ix] == a[ef_ix-1]:
            sorted_a_count[-1] += 1
    return((sorted_a, sorted_a_count))

def potatoCrossGenotyper(args):
    ## Get the VCF file (filtered may be) generated by GATK.
    ## inputs:
    # 1) VCF file
    # 2) Parent1 and Parent2
    # 3) SNP matrix (hdf5 file)
    # 4) Bin length, default as 200Kbp
    # 5) Chromosome length
    global genome
    genome = genomes.Genome(args['genome'])
    log.info("loading database files")
    g = snp_genotype.Genotype(args['hdf5File'], args['hdf5accFile'])
    log.info("done!")
    crossgenotyper = GenotypeCross(g, args['parents'], args['binLen'], args['father'], args['logDebug'])
    if args['hmm']:
        outfile_str = crossgenotyper.genotype_cross_hmm( args['inFile'] )        
    else:
        outfile_str = crossgenotyper.genotype_cross( args['inFile'], args['lr_thres'] )
    crossgenotyper.write_output_genotype_cross( outfile_str, args['outFile'] )
