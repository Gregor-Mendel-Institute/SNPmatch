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
        segSNPsind = np.where((snpsP1 != snpsP2) & (snpsP1 >= 0) & (snpsP2 >= 0) & (snpsP1 < 2) & (snpsP2 < 2))[0]
        # segSNPsind = np.where((snpsP1 != snpsP2) )[0]
        log.info("number of segregating snps between parents: %s", len(segSNPsind))
        self.commonSNPsCHR = commonSNPsCHR[segSNPsind].astype('U')
        self.commonSNPsPOS = commonSNPsPOS[segSNPsind]
        self.snpsP1 = snpsP1[segSNPsind]
        self.snpsP2 = snpsP2[segSNPsind]
        log.info("done!")

    def _get_common_positions_ixs(self, accs_pos, samples_pos, accs_ix, sample_ix):
        ## Make sure you have inputs loaded
        reqAccsPOS = accs_pos[accs_ix]
        reqTarPos = samples_pos[sample_ix]
        commonPOS = np.intersect1d( reqAccsPOS, reqTarPos )
        matchedAccInd = np.array(accs_ix)[ np.where( np.in1d(reqAccsPOS, reqTarPos) )[0] ]
        matchedTarInd = np.array(sample_ix)[ np.where( np.in1d(reqTarPos, reqAccsPOS) )[0] ]
        return((matchedAccInd, matchedTarInd))
    
    @staticmethod
    def get_parental_obs(input_gt, snpsP1_gt, snpsP2_gt, polarize = None):
        num_snps = len(input_gt)
        ebTarGTs = parsers.parseGT(input_gt)
        ebPolarised = np.zeros(num_snps, dtype=int)
        ### Here I am having 6 observed states. to consider SNPs where either of the parent shows -1 
        ## ('00', '01', '11', '0', '1', 'NA') 
        ##    0,    1,    2,   3,   4,    5
        ebPolarised[:] = 5
        ebPolarised[np.where(np.equal( ebTarGTs, np.repeat(2, num_snps) ))[0] ] = 1
        snpsP1_gt_mask = numpy.ma.masked_less(numpy.ma.masked_greater(snpsP1_gt, 1), 0)
        snpsP2_gt_mask = numpy.ma.masked_less(numpy.ma.masked_greater(snpsP2_gt, 1), 0)
        ebPolarised[np.where((np.equal( ebTarGTs, snpsP1_gt_mask )) & (~snpsP2_gt_mask) )[0] ] = 0  ## 00
        ebPolarised[np.where((np.equal( ebTarGTs, snpsP1_gt_mask )) & (snpsP2_gt_mask) )[0] ] = 3 ## 0
        ebPolarised[np.where((np.equal( ebTarGTs, snpsP2_gt_mask )) & (~snpsP1_gt_mask) )[0] ] = 2  ## 11
        ebPolarised[np.where((np.equal( ebTarGTs, snpsP2_gt_mask )) & (snpsP1_gt_mask) )[0] ] = 4 ## 1
        return(ebPolarised)
    
    @staticmethod
    def init_hmm(num_markers, chromosome_size = 115, recomb_rate = 3.3):
        """
        Function to initilize a HMM model 
        Input: 
            number of markers
            recomb_rate: assumed rate of recombination (cM per Mb)
            genome_size: size of chromosome in Mb
        """
        from hmmlearn import hmm
        ### The below transition probability is for a intercross, adapted from R/qtl
        log.info("Initialising HMM")
        states = ('aa', 'ab', 'bb')
        observations = ('00', '01', '11', '0', '1', 'NA')
        ## assume A. thaliana genome size of 115 Mb 
        ## First calculate the average spacing between the markers
        ## chromosome_size / num_markers 
        ##_____
        ## given a recombination rate of ~3.3 
        # probability you see a recombinant is 1 / (100 * recomb_rate)
        prob_of_change = ( chromosome_size / num_markers ) * ( 1 / (100 * recomb_rate)  )
        transmission_prob = {
            'aa': {'aa': 1 - prob_of_change, 'ab': prob_of_change/2, 'bb': prob_of_change/2},
            'ab': {'aa': prob_of_change/2, 'ab': 1 - prob_of_change, 'bb': prob_of_change/2},
            'bb': {'aa': prob_of_change/2, 'ab': prob_of_change/2, 'bb': 1 - prob_of_change}
        }
        ## Since I have 6 possible observed states -- I have made an emmission matrix with such a structure.
        emission_prob = {
            'aa': {'00': 0.599, '01': 0.1, '11': 0.001, '0': 0.52, '1': 0, 'NA': 0.33},
            'ab': {'00': 0.4, '01': 0.8, '11': 0.4, '0': 0.48, '1': 0.48, 'NA': 0.34},
            'bb': {'00': 0.001, '01': 0.1, '11': 0.599, '0': 0, '1': 0.52, 'NA': 0.33}
        }
        model = hmm.MultinomialHMM(n_components=3) 
        model.startprob_ = np.array([0.25, 0.5, 0.25])
        model.transmat_ = pd.DataFrame(transmission_prob)
        model.emissionprob_ = pd.DataFrame(emission_prob).T
        return(model)

    def genotype_cross_hmm(self, input_file):
        """
        Function to genotype an F2 individual based on segregating SNPs using a simple HMM model
        """
        snpvcf = parsers.import_vcf_file(inFile = input_file, logDebug = self.logDebug, samples_to_load = None)
        samples_ids = pd.Series(snpvcf['samples']) #.str.replace("_processed_reads_no_clonal.bam", "" )
        # if np.unique(samples_ids.str.split("_", expand = True).iloc[:,0]).shape[0] == samples_ids.shape[0]:
        #     samples_ids = samples_ids.str.split("_", expand = True).iloc[:,0]
        samples_gt = snpvcf['gt']
        samples_gt = pd.DataFrame(samples_gt.astype(str))
        input_chr_names = genome.chrs[pd.Series(snpvcf['chr'], dtype = str).apply(genome.get_chr_ind)]
        g_chr_names = genome.chrs[pd.Series(self.commonSNPsCHR, dtype = str).apply(genome.get_chr_ind)]
        allSNPGenos = pd.DataFrame( - np.ones((0,samples_gt.shape[1]), dtype=int) )
        allSNPGenos_raw = pd.DataFrame( - np.ones((0,samples_gt.shape[1]), dtype=int) )
        for ec, eclen in zip(genome.chrs_ids, genome.chrlen):
            reqPOSind =  np.where(g_chr_names == ec)[0]
            reqTARind = np.where( input_chr_names == ec )[0]
            ebAccsInds, ebTarInds = self._get_common_positions_ixs(self.commonSNPsPOS, snpvcf['pos'], reqPOSind, reqTARind)
            t_model = self.init_hmm( len(ebAccsInds), eclen / 1000000 )
            t_genotypes = pd.DataFrame( - np.ones((len(ebAccsInds),samples_gt.shape[1]), dtype=int), index = ec + ":" + pd.Series(self.commonSNPsPOS[ebAccsInds]).astype(str) )
            t_genotypes_raw = pd.DataFrame( - np.ones((len(ebAccsInds),samples_gt.shape[1]), dtype=int), index = ec + ":" + pd.Series(self.commonSNPsPOS[ebAccsInds]).astype(str) )
            for sample_ix in range(samples_gt.shape[1]):
                ebin_gt_polarized = self.get_parental_obs(samples_gt.iloc[ebTarInds,sample_ix].values, self.snpsP1[ebAccsInds], self.snpsP2[ebAccsInds])
                t_genotypes_raw.iloc[:, sample_ix] = ebin_gt_polarized
                t_genotypes.iloc[:, sample_ix] = t_model.predict(ebin_gt_polarized.reshape((-1, 1)))
            allSNPGenos = allSNPGenos.append(t_genotypes)
            allSNPGenos_raw = allSNPGenos_raw.append(t_genotypes_raw)
            import ipdb; ipdb.set_trace()
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
    sorted_a_count = np.zeros(0, dtype = int)
    sorted_a = np.zeros(0, dtype = a.dtype)
    for ef_ix in range(0, len(a)):
        if a[ef_ix] != a[ef_ix-1]:
            sorted_a = np.append(sorted_a, a[ef_ix])
            sorted_a_count = np.append(sorted_a_count, 1)
        elif np.sum(sorted_a_count) == 0:
            sorted_a = np.append(sorted_a, a[ef_ix])
            sorted_a_count[-1] += 1
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
