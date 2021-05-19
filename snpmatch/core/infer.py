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

# code adapted from Stephen Marsland's, Machine Learning An Algorthmic Perspective, Vol. 2
# https://github.com/alexsosn/MarslandMLAlgo/blob/master/Ch16/HMM.py
def viterbi(init_prob, trans_mat, emission_mat, obs):
    """
    Parameters:
        init_prob: initial probabilities
        trans_mat : transition matrix
        emission_mat : emission probabilty matrix
        obs: observations
    """
    nStates = np.shape(init_prob)[0]
    T = np.shape(obs)[0]
    
    # init blank path
    path = np.zeros(T)
    # delta --> highest probability of any path that reaches state i
    delta = np.zeros((nStates, T))
    # phi --> argmax by time step for each state
    phi = np.zeros((nStates, T))
    
    if len(emission_mat.shape) == 2:
        emission_mat = np.tile( emission_mat.T, (T, 1, 1) ).T
    
    # init delta and phi 
    delta[:, 0] = init_prob * emission_mat[:, obs[0], 0]
    phi[:, 0] = 0

    # print('\nStart Walk Forward\n')    
    # the forward algorithm extension
    for t in range(1, T):
        for s in range(nStates):
            delta[s, t] = np.max(delta[:, t-1] * trans_mat[:, s]) * emission_mat[s, obs[t], t] 
            phi[s, t] = np.argmax(delta[:, t-1] * trans_mat[:, s])
            # print('s={s} and t={t}: phi[{s}, {t}] = {phi}'.format(s=s, t=t, phi=phi[s, t]))
    # find optimal path
    # print('-'*50)
    # print('Start Backtrace\n')
    path[T-1] = np.argmax(delta[:, T-1])
    # p('init path\n    t={} path[{}-1]={}\n'.format(T-1, T, path[T-1]))
    for t in range(T-2, -1, -1):
        path[t] = phi[int(path[t+1]), [t+1]]
        # p(' '*4 + 't={t}, path[{t}+1]={path}, [{t}+1]={i}'.format(t=t, path=path[t+1], i=[t+1]))
        # print('path[{}] = {}'.format(t, path[t]))
    return((path, delta, phi))

class IdentifyStrechesofHeterozygosity(object):

    def __init__(self,
        num_markers, 
        chromosome_size, 
        avg_depth = 1.5, 
        fraction_homo_parents = 0.99,
        avg_sites_segregating = 0.01, 
        base_error = 0.0001, 
        recomb_rate = 3.3, 
        n_iter = 1000
    ):
        #  
        self.params = {}
        self.params['num_markers'] = num_markers
        self.params['avg_depth'] = avg_depth
        ### fraction of sites parental genomes are heterozygous
        self.params['delta_het_parents'] = fraction_homo_parents 
        ### fraction of sites which are segregating between parents
        self.params['avg_sites_segregating'] = avg_sites_segregating 
        ### sequencing error rate
        self.params['base_error'] = base_error 
        ## genome size 
        self.params['chromosome_size'] = chromosome_size
        ## recombination rate
        self.params['recomb_rate'] = recomb_rate
        ## number of iterations
        self.params['n_iter'] = n_iter
        self._init_hmm_streches_of_het()


    def _init_hmm_streches_of_het(self):
        """
        HMM to identify streches of homozygosity in an intercross
        ## currently: observations should only be in segregating SNPs
        Input: 
            number of markers
            avg_sites_segregating:  Average number of given sites that would be segregating between parents
            recomb_rate: assumed rate of recombination (cM per Mb)
            genome_size: size of chromosome in Mb
        """
        ### There are two hidden states
        log.info("Initialising HMM")
        states = ('HOMO', 'HET')
        ## you observe 
        observations = ('00||11', '01', 'NA')
        ## if a site is called Homo -- all the reads mapped to it are the same
        ## het if atleast one read is different (default GATK/bcftools callers does that)
        ##_____
        ## What is the transition probabilities? 
        ## First, I am doing it for F2 cross
        self.init_prob = [0.5, 0.5]
        ri = (float(self.params['chromosome_size']) / self.params['num_markers']) * self.params['recomb_rate'] / 100
        ## assuming a diploid
        transition_prob = [
                [(1 - ri)**2 + ri**2,   2 * ri * (1 - ri)   ],
                [2 * ri * (1 - ri),     (1 - ri)**2 + ri**2 ]
        ]
        ## Emission probability
        # average_depth
        # base_error 
        # avg_sites_segregating 
        ## Z -- underlying ancestry -- either AA or AB
        ## G = genotype at a given locus -- aa or ab -- depends on whether parents are segregating
        ## X = observed states (00, 01)
        p_homo_given_gaa = ((1 - self.params['base_error']) ** self.params['avg_depth']) + self.params['base_error']**self.params['avg_depth']
        p_homo_given_gab = 2 * (0.5 ** self.params['avg_depth'])
        ## What fraction of sites are heterozygoues in parental genomes. 
        ## delta_het_parents
        p_homo_given_ZAA = (self.params['delta_het_parents'] * p_homo_given_gaa ) + ((1 - self.params['delta_het_parents']) * p_homo_given_gab )
        p_homo_given_ZAB = ((1 - self.params['avg_sites_segregating']) * p_homo_given_gaa ) + ( self.params['avg_sites_segregating'] * p_homo_given_gab )

        emission_prob = [ 
            [p_homo_given_ZAA,  1 - p_homo_given_ZAA, 0.5],
            [p_homo_given_ZAB, 1 - p_homo_given_ZAB, 0.5]
        ]
        print("transition probability:\n %s" % pd.DataFrame(transition_prob,index = states, columns = states) )
        self.transition_prob = np.array(transition_prob)
        print("emission probability:\n %s" % pd.DataFrame(emission_prob, index = states, columns = observations) )
        self.emission_prob = np.array(emission_prob)
        
    def viterbi(self, obs):
        log.info("initialising HMM")
        model = viterbi( self.init_prob, self.transition_prob, self.emission_prob, obs )
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
        num_markers,
        recomb_rate = 3.5,  ### in cM/Mb
        error_p1 = 0.00001,
        error_p2 = 0.00001, 
        hetfrac_p1 = 0.001,
        hetfrac_p2 = 0.001, 
        base_error = 0.01, 
        avg_depth = 1.5,
        n_iter = 1000
    ):
        self.ancestry = ['AA', 'AB', 'BB']
        self.geno_parents = ['00', '01', '11']
        self.observed_states = ['00', '01', '11', 'NA']
        self.params = {}
        self.params['num_markers'] = num_markers
        self.params['chromosome_size'] = chromosome_size
        self.params['recomb_rate'] = recomb_rate
        self.params['error_p1'] = error_p1
        self.params['error_p2'] = error_p2
        self.params['hetfrac_p1'] = hetfrac_p1
        self.params['hetfrac_p2'] = hetfrac_p2
        self.params['base_error'] = base_error
        self.params['avg_depth'] = avg_depth
        self.params['n_iter'] = n_iter
        self.prob_g_given_z = self._prob_g_given_Z(self.params['error_p1'], self.params['error_p2'], self.params['hetfrac_p1'], self.params['hetfrac_p2'] )
        self.prob_x_given_g = self._prob_x_given_G(self.params['base_error'], self.params['avg_depth'] )
        self.init_prob = [0.25, 0.5, 0.25] ### F2 individual with -- Mendelian segregation
        print("init probabilites:\n %s" % pd.Series(self.init_prob, index = self.ancestry) )
        self.transition_prob = self._transition_prob(self.params['chromosome_size'], self.params['num_markers'], self.params['recomb_rate'])
        print("transition probability:\n %s" % self.transition_prob )
        self.emission_prob =  pd.DataFrame(np.dot(self.prob_g_given_z.values, self.prob_x_given_g.values), index = self.ancestry, columns = self.observed_states)
        print("emission probability:\n %s" % self.emission_prob )

    def _prob_g_given_Z(self, error_p1, error_p2, hetfrac_p1, hetfrac_p2):
        ## Adapted from Andolfatto et al.
        conf_p1 = 1 - error_p1
        conf_p2 = 1 - error_p2
        prob_het_p1 = hetfrac_p1**2
        prob_het_p2 = hetfrac_p2**2
        prob_00_given_aa = (conf_p1**2 * (1 - hetfrac_p1)) + (hetfrac_p1 * conf_p1 * error_p1)
        prob_11_given_aa = (error_p1**2 * (1 - hetfrac_p1)) + (hetfrac_p1 * conf_p1 * error_p1)
        prob_00_given_bb = (error_p2**2 * (1 - hetfrac_p2)) + (hetfrac_p2 * conf_p2 * error_p2)
        prob_11_given_bb = (conf_p2**2 * (1 - hetfrac_p2)) + (hetfrac_p2 * conf_p2 * error_p2)
        prob_00_given_ab = (((1 - hetfrac_p1) * conf_p1) + (hetfrac_p1 * error_p1) ) * (((1 - hetfrac_p2) * error_p2) + (hetfrac_p2 * conf_p2) )
        prob_11_given_ab = (((1 - hetfrac_p2) * conf_p2) + (hetfrac_p2 * error_p2) ) * (((1 - hetfrac_p1) * error_p1) + (hetfrac_p1 * conf_p1) )
        req_prob = [
            [prob_00_given_aa, 1 - prob_00_given_aa - prob_11_given_aa, prob_11_given_aa],
            [prob_00_given_ab, 1 - prob_11_given_ab - prob_00_given_ab, prob_11_given_ab],
            [prob_00_given_bb, 1 - prob_00_given_bb - prob_11_given_bb, prob_11_given_bb]
        ]
        return( pd.DataFrame( req_prob, index = self.ancestry, columns = self.geno_parents ) )
    
    def _prob_x_given_G(self, base_error, avg_depth):
        req_prob = [
            [(1 - base_error)**avg_depth, 1 - base_error**avg_depth - (1 - base_error)**avg_depth, base_error**avg_depth, 0.5],
            [(0.5 * (1 - base_error))**avg_depth, 1 - 2 * ((0.5 * (1 - base_error))**avg_depth),(0.5 * (1 - base_error))**avg_depth,0.5],
            [base_error**avg_depth, 1 - base_error**avg_depth - (1 - base_error)**avg_depth, (1-base_error)**avg_depth, 0.5]
        ]
        return( pd.DataFrame( req_prob, index = self.geno_parents, columns = self.observed_states )  )

    def _transition_prob(self, chromosome_size, num_markers, recomb_rate):
        ri = (float(chromosome_size) / num_markers) * recomb_rate / 100
        ## assuming a diploid
        transition_prob = [
                [(1 - ri)**2,   2 * ri * (1 - ri), ri**2],
                [ri * (1 - ri),     (1 - ri)**2 + ri**2, ri * (1 - ri)],
                [ri**2,   2 * ri * (1 - ri), (1 - ri)**2]
        ]
        return( pd.DataFrame( transition_prob, index = self.ancestry, columns = self.ancestry ) )
    
    def viterbi(self, obs, n_iter = 1000):
        log.info("initialising HMM")
        model = viterbi( self.init_prob, self.transition_prob.values, self.emission_prob.values, obs )
        return(model)

    @staticmethod
    def snp_to_observations(input_snps, snpsP1_gt, snpsP2_gt, polarize = None):
        # input_snps = parsers.parseGT(input_gt)
        num_snps = len(input_snps)
        ### Here I have 4 observed states
        ## ('00', '01', '11',  'NA') ### removed '0', '1',
        ##    0,    1,    2,   3  ###  4,    5
        ebPolarised = np.repeat(3, num_snps)
        input_snps_mask = numpy.ma.masked_less(input_snps, 0)
        snpsP1_gt_mask = numpy.ma.masked_less(numpy.ma.masked_greater(snpsP1_gt, 1), 0)
        snpsP2_gt_mask = numpy.ma.masked_less(numpy.ma.masked_greater(snpsP2_gt, 1), 0)
        if polarize == "p1":
            ebPolarised[np.where((np.equal( input_snps_mask, snpsP1_gt_mask )) )[0] ] = 0  ## 00
            ebPolarised[np.where((~np.equal( input_snps_mask, snpsP1_gt_mask )) & (input_snps_mask < 2) )[0] ] = 2  ## 11
        elif polarize == "p2":
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