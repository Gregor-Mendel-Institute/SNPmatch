import numpy as np
import scipy as sp
import pandas as pd
import logging
from glob import glob
import itertools
import os.path
from hmmlearn import hmm
import numpy.ma



log = logging.getLogger(__name__)

# viterbi code adapted from 
# http://www.adeveloperdiary.com/data-science/machine-learning/implement-viterbi-algorithm-in-hidden-markov-model-using-python-and-r/
def viterbi(init_prob, trans_mat, emission_mat, obs):
    """
    Parameters:
        init_prob: initial probabilities
        trans_mat : transition matrix
        emission_mat : emission probabilty matrix
        obs: observations
    """
    T = obs.shape[0]
    M = trans_mat.shape[0]

    if len(emission_mat.shape) == 2:
        emission_mat = np.tile( emission_mat.T, (T, 1, 1) ).T
 
    omega = np.zeros((T, M))
    omega[0, :] = np.log(init_prob * emission_mat[:, obs[0],0])
    prev = np.zeros((T - 1, M))
    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            probability = omega[t - 1] + np.log(trans_mat[:, j]) + np.log(emission_mat[j, obs[t], t])
            # This is our most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(probability)
            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)
 
    # Path Array
    S = np.zeros(T)
    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])
 
    S[0] = last_state
 
    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1
    
    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)
    return((S, omega))

class IdentifyStrechesofHeterozygosity(object):
    """
    HMM to identify streches of homozygosity in an intercross
    ## currently: observations should only be in segregating SNPs
    Input: 
        number of markers
        avg_sites_segregating:  Average number of given sites that would be segregating between parents
        recomb_rate: assumed rate of recombination (cM per Mb)
        genome_size: size of chromosome in Mb
    """

    def __init__(self,
        num_markers, 
        chromosome_size, 
        sample_depth = 1.5, 
        fraction_homo_parents = 0.99,
        avg_sites_segregating = 0.01, 
        base_error = 0.0001, 
        recomb_rate = 3.3
    ):
        #  
        self.params = {}
        self.params['num_markers'] = num_markers
        if isinstance(sample_depth, (int, float)):
            self.params['sample_depth'] = np.repeat(sample_depth, self.params['num_markers'])
        else:
            self.params['sample_depth'] = sample_depth
        ### fraction of sites parental genomes are homozygous
        self.params['fraction_homo_parents'] = fraction_homo_parents 
        ### fraction of sites which are segregating between parents
        self.params['avg_sites_segregating'] = avg_sites_segregating 
        ### sequencing error rate
        self.params['base_error'] = base_error 
        ## genome size 
        self.params['chromosome_size'] = chromosome_size
        ## recombination rate
        self.params['recomb_rate'] = recomb_rate
        ## number of iterations
        log.info("Initialising HMM")
        ### There are two hidden states
        self.hidden_states = ('HOMO', 'HET')
        ## you observe 
        self.observed_states = ('00||11', '01', 'NA')
        ## First, I am doing it for F2 cross
        self.init_prob = [0.5, 0.5]
        self.transition_prob = self.calc_transition_prob(self.params['num_markers'], self.params['recomb_rate'], self.params['chromosome_size'])
        log.info("transition probability:\n %s" % self.transition_prob )
        log.info("calculating emissions")
        self.emission_prob = self.calc_emissions(self.params['base_error'], self.params['sample_depth'], self.params['fraction_homo_parents'], self.params['avg_sites_segregating'])

    def calc_transition_prob(self, num_markers, recomb_rate, chromosome_size):
        ri = (float(chromosome_size) / num_markers) * recomb_rate / 100
        ## assuming a diploid
        transition_prob = np.array([
                [(1 - ri)**2 + ri**2,   2 * ri * (1 - ri)   ],
                [2 * ri * (1 - ri),     (1 - ri)**2 + ri**2 ]
        ])
        return( pd.DataFrame( transition_prob, index = self.hidden_states, columns = self.hidden_states ) )

    def calc_emissions(self, base_error, sample_depth, fraction_homo_parents, avg_sites_segregating):
        ## if a site is called Homo -- all the reads mapped to it are the same
        ## het if atleast one read is different (default GATK/bcftools callers does that)
        ##_____
        ## Emission probability
        # average_depth
        # base_error 
        # avg_sites_segregating 
        ## Z -- underlying ancestry -- either AA or AB
        ## G = genotype at a given locus -- 00 or 01 -- depends on whether parents are segregating
        ## X = observed states (00, 01)
        emission_prob = np.zeros( (len(self.hidden_states), len(self.observed_states), len(sample_depth)) )
        ## What fraction of sites are heterozygoues in parental genomes. 
        ## fraction_homo_parents
        ## Rows are Z -- AA, AB
        prob_g_given_Z = np.array([
            [fraction_homo_parents, 1 - fraction_homo_parents],
            [1 - avg_sites_segregating, avg_sites_segregating],
        ])
        iter_depth = np.unique( sample_depth ) 
        for ef_depth in iter_depth: 
            req_snp_ix = np.where( sample_depth == ef_depth )[0]
            p_homo_given_gaa = ((1 - base_error) ** ef_depth) + base_error**ef_depth
            p_homo_given_gab = 2 * (0.5**ef_depth)
            t_prob_x_given_g = np.array([
                [p_homo_given_gaa, 1 - p_homo_given_gaa, 1],
                [p_homo_given_gab, 1 - p_homo_given_gab, 1]
            ])
            
            ef_emission = np.dot(prob_g_given_Z, np.abs(t_prob_x_given_g) )
            emission_prob[:,:,req_snp_ix] = np.tile( ef_emission.T, (len(req_snp_ix), 1, 1) ).T
        return( np.array( emission_prob ) )
        
    def viterbi(self, input_snps):
        log.info("initialising HMM")
        obs = self.snp_to_observations( input_snps )
        model = viterbi( self.init_prob, self.transition_prob.values, self.emission_prob, obs )
        return(model)

    @staticmethod
    def snp_to_observations(input_snps):
        t_obs = np.array(input_snps)
        t_obs[t_obs == 1] = 0
        t_obs[t_obs == 2] = 1
        t_obs[t_obs == -1] = 2
        return(t_obs)



class IdentifyAncestryF2individual(object):

    def __init__(self, 
        chromosome_size,
        snps_p1,
        snps_p2, 
        recomb_rate = 3.5,  ### in cM/Mb
        error_p1 = 0.00001,
        error_p2 = 0.00001, 
        base_error = 0.01, 
        sample_depth = 1.5
    ):
        log.info("initialising HMM")
        self.ancestry = ['AA', 'AB', 'BB']
        self.geno_parents = ['00', '01', '11']
        self.observed_states = ['00', '01', '11', 'NA']
        log.info( "observations: %s" % self.observed_states )
        self.params = {}
        assert snps_p1.shape[0] == snps_p2.shape[0], "both the SNP arrays for two parents should be of same size"
        self.params['num_markers'] = snps_p1.shape[0]
        self.params['chromosome_size'] = chromosome_size
        self.params['recomb_rate'] = recomb_rate
        self.params['error_p1'] = error_p1
        self.params['error_p2'] = error_p2
        self.params['snps_p1'] = snps_p1
        self.params['snps_p2'] = snps_p2
        self.params['base_error'] = base_error
        if isinstance(sample_depth, (int, float)):
            self.params['sample_depth'] = np.repeat(sample_depth, self.params['num_markers'])
        else:
            self.params['sample_depth'] = sample_depth
        self.init_prob = [0.25, 0.5, 0.25] ### F2 individual with -- Mendelian segregation
        log.info("init probabilites:\n %s" % pd.Series(self.init_prob, index = self.ancestry) )
        self.transition_prob = self._transition_prob(self.params['chromosome_size'], self.params['num_markers'], self.params['recomb_rate'])
        log.info("transition probability:\n %s" % self.transition_prob )
        log.info("calculating emissions")
        self.emission_prob = self._get_emissions(self.params['error_p1'], self.params['error_p2'], self.params['snps_p1'], self.params['snps_p2'], self.params['base_error'], self.params['sample_depth'])

    def _get_emissions(self, error_p1, error_p2, snps_p1, snps_p2, base_error, sample_depth):
        unique_p1 = np.unique(snps_p1)
        unique_p2 = np.unique(snps_p2)
        emission_prob = np.zeros( (len(self.ancestry), len(self.observed_states), len(snps_p1)) )
        iter_parameter = pd.DataFrame( { "snps_p1": snps_p1, "snps_p2": snps_p2, "dp": sample_depth } ).drop_duplicates()
        for ef_param in iter_parameter.iterrows(): 
            ef_emission = self._calc_emission_given_af( 
                error_p1, 
                error_p2, 
                get_af(ef_param[1]['snps_p1']), 
                get_af(ef_param[1]['snps_p2']), 
                base_error, 
                ef_param[1]['dp'] 
            )
            req_snp_ix = np.where((snps_p1 == ef_param[1]['snps_p1']) & ( snps_p2 == ef_param[1]['snps_p2'] ) & (sample_depth == ef_param[1]['dp']) )[0]
            # for ef_ix in req_snp_ix:
                # emission_prob[:,:,ef_ix] = ef_emission.values
            emission_prob[:,:,req_snp_ix] = np.tile( ef_emission.values.T, (len(req_snp_ix), 1, 1) ).T
        return( emission_prob )
    
    def _calc_emission_given_af(self, error_p1, error_p2, af_p1, af_p2, base_error, avg_depth):
        """
        Calculate emission probability at a given marker with
        error_p* : estimate for sequencing error for parent *
        af_p* : allele frequncy for two alleles (say reference and alternative) for the parent
                for example: 0 for homo (00), 0.5 for het (01) and 1 for homo (11)
        base_error : base error for the sample
        avg_depth: depth for the sample at a given position
        """
        ## Adapted from Andolfatto et al.
        conf_p1 = 1 - error_p1
        conf_p2 = 1 - error_p2
        prob_00_given_aa = (conf_p1**2 * (1 - af_p1)) + (error_p1**2 * af_p1)
        prob_11_given_aa = (conf_p1**2 * af_p1) + (error_p1**2 * (1 - af_p1))

        prob_00_given_bb = (conf_p2**2 * (1 - af_p2)) + (error_p2**2 * af_p2)
        prob_11_given_bb = (conf_p2**2 * af_p2) + (error_p2**2 * (1 - af_p2))

        prob_00_given_ab = (((1-af_p1) * conf_p1) + (af_p1 * error_p1) ) * (((1 - af_p2) * conf_p2) + (af_p2 * error_p2) )
        prob_11_given_ab = ((af_p1 * conf_p1) + ((1-af_p1) * error_p1) ) * ((af_p2 * conf_p2) + ((1 - af_p2) * error_p2))

        prob_g_given_Z = [
            [prob_00_given_aa, 1 - prob_00_given_aa - prob_11_given_aa, prob_11_given_aa],
            [prob_00_given_ab, 1 - prob_11_given_ab - prob_00_given_ab, prob_11_given_ab],
            [prob_00_given_bb, 1 - prob_00_given_bb - prob_11_given_bb, prob_11_given_bb]
        ]
        ## if depth is 2 -- p(00|00) = 0.99**2 , P(11|00) = 0.01**2 
        p_00_given_g00 = (1 - base_error)**avg_depth ## if depth is 10 -- 0.99 ** 10
        p_11_given_g00 = base_error ** avg_depth
        p_01_given_g00 = 1 - p_00_given_g00 - p_11_given_g00
        # qbe = 0.5*base_error
        # qb = 0.5*(1-base_error)
        # p_00_given_g01 = (0.5)**avg_depth 
        #  (avg_depth,0) qb**avg_depth + (avg_depth,1) (qb**(avg_depth-1) * qbe) + avg_depth C2 (qb**(avg_depth-2) * qbe**2)
        p_01_given_g01 = 1 - 2 * (0.5**avg_depth)
        p_00_given_g01 = (1 - p_01_given_g01)/2
        prob_x_given_g = [
            [p_00_given_g00, p_01_given_g00, p_11_given_g00, 1],
            [p_00_given_g01, p_01_given_g01, p_00_given_g01, 1],
            [p_11_given_g00, p_01_given_g00, p_00_given_g00, 1]
        ]
        ## you need to do absolute as with depth of 0, you get probabilites as -1
        prob_x_given_g = np.abs(np.array(prob_x_given_g)) 
        return( pd.DataFrame(np.dot(np.array(prob_g_given_Z), prob_x_given_g), index = self.ancestry, columns = self.observed_states) ) 

    def _transition_prob(self, chromosome_size, num_markers, recomb_rate):
        ri = (float(chromosome_size) / num_markers) * recomb_rate / 100
        ## assuming a diploid
        transition_prob = [
                [(1 - ri)**2,   2 * ri * (1 - ri), ri**2],
                [ri * (1 - ri),     (1 - ri)**2 + ri**2, ri * (1 - ri)],
                [ri**2,   2 * ri * (1 - ri), (1 - ri)**2]
        ]
        return( pd.DataFrame( transition_prob, index = self.ancestry, columns = self.ancestry ) )
    
    def viterbi(self, input_snps):
        obs = self.snp_to_observations( input_snps )
        model = viterbi( self.init_prob, self.transition_prob.values, self.emission_prob, obs )
        return(model)

    @staticmethod
    def snp_to_observations(input_snps):
        # input_snps = parsers.parseGT(input_gt)
        num_snps = len(input_snps)
        ### Here I have 4 observed states
        ## ('00', '01', '11',  'NA') ### removed '0', '1',
        ##    0,    1,    2,   3  ###  4,    5
        t_sample_snps = np.copy(input_snps)
        t_sample_snps[t_sample_snps == -1] = 3
        t_sample_snps[t_sample_snps == 2] = 5
        t_sample_snps[t_sample_snps == 1] = 2
        t_sample_snps[t_sample_snps == 5] = 1
        return(t_sample_snps)

def get_af(snps):
    t_af = np.copy(snps)
    t_af[t_af == 1] = 4
    t_af[t_af == 2] = 1
    t_af[t_af == 4] = 2
    t_af = t_af / 2
    return(t_af)

def polarize_snps(input_snps, snps_p1, snps_p2, polarize_to = None):
    # input_snps = parsers.parseGT(input_gt)
    num_snps = len(input_snps)
    ### Here I have 4 observed states
    ## ('00', '01', '11',  'NA') ### removed '0', '1',
    ##    0,    1,    2,   3  ###  4,    5
    ebPolarised = np.repeat(3, num_snps)
    input_snps_mask = numpy.ma.masked_less(input_snps, 0)
    snpsP1_gt_mask = numpy.ma.masked_less(numpy.ma.masked_greater(snps_p1, 1), 0)
    snpsP2_gt_mask = numpy.ma.masked_less(numpy.ma.masked_greater(snps_p2, 1), 0)
    if polarize_to == "p1":
        ebPolarised[np.where((np.equal( input_snps_mask, snpsP1_gt_mask )) )[0] ] = 0  ## 00
        ebPolarised[np.where((~np.equal( input_snps_mask, snpsP1_gt_mask )) & (input_snps_mask < 2) )[0] ] = 2  ## 11
    elif polarize_to == "p2":
        ebPolarised[np.where((np.equal( input_snps_mask, snpsP2_gt_mask )) )[0] ] = 2  ## 00
        ebPolarised[np.where((~np.equal( input_snps_mask, snpsP2_gt_mask )) & ( input_snps_mask < 2 ) )[0] ] = 0  ## 11
    else:
        ebPolarised[np.where((np.equal( input_snps_mask, snpsP1_gt_mask )) )[0] ] = 0  ## 00
        ebPolarised[np.where((np.equal( input_snps_mask, snpsP2_gt_mask )) )[0] ] = 2  ## 11
    ebPolarised[np.where( np.equal( input_snps_mask, np.repeat(2, num_snps)) & (snpsP1_gt_mask != snpsP2_gt_mask ) )[0] ] = 1
    return(ebPolarised)


def uniq_neighbor(a):
    """
    Function to club nearby positions and give counts
    useful in getting the recombination break points
    input:
        numpy array (1d) for observations.
    """
    sorted_a = np.array(a[0], dtype = a.dtype)
    sorted_a_count = np.array([1], dtype = int)
    for ef_ix in range(1, len(a)):
        if a[ef_ix] != a[ef_ix-1]:
            sorted_a = np.append(sorted_a, a[ef_ix])
            sorted_a_count = np.append(sorted_a_count, 1)
        elif a[ef_ix] == a[ef_ix-1]:
            sorted_a_count[-1] += 1
    return((sorted_a, sorted_a_count))