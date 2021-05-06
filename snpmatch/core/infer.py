import numpy as np
import scipy as sp
import pandas as pd
import logging
from glob import glob
import itertools
import os.path
from . import parsers



log = logging.getLogger(__name__)

class IdentifyStrechesofHeterozygosity(object):

    def __init__(self,
        sample_gt, 
        chromosome_size, 
        avg_depth = 1.5, 
        delta_het_parents = 0.99,
        avg_sites_segregating = 0.01, 
        base_error = 0.01, 
        recomb_rate = 3.3, 
        n_iter = 1000
    ):
        #  
        self.sample_gt = parsers.parseGT( sample_gt )
        self.params = {}
        self.params['num_markers'] = sample_gt.shape[0]
        self.params['avg_depth'] = avg_depth
        ### fraction of sites parental genomes are heterozygous
        self.params['delta_het_parents'] = delta_het_parents 
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

        self.model = self._init_hmm_streches_of_het()


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
        from hmmlearn import hmm
        ### There are two hidden states
        log.info("Initialising HMM")
        states = ('aa', 'ab')
        ## you observe 
        observations = ('homo', 'het', 'NA')
        ## if a site is called Homo -- all the reads mapped to it are the same
        ## het if atleast one read is different (default GATK/bcftools callers does that)
        ##_____
        ## What is the transition probabilities? 
        ## First, I am doing it for F2 cross
        init_prob = [0.5, 0.5]
        prob_of_change = ( self.params['chromosome_size'] / self.params['num_markers'] ) * ( 1 / (100 * self.params['recomb_rate'])  )
        transition_prob = [
                [1 - prob_of_change,  prob_of_change],
                [prob_of_change,  1 - prob_of_change]
        ]
        ## Emission probability
        # average_depth
        # base_error 
        # avg_sites_segregating 
        ## Z -- underlying ancestry -- either AA or AB
        ## G = genotype at a given locus -- aa or ab -- depends on whether parents are segregating
        ## X = observed states (00, 01)
        p_homo_given_gaa = (1 - self.params['base_error']) ** self.params['avg_depth']
        p_homo_given_gab = 2 * (0.5 ** self.params['avg_depth'])
        ## What fraction of sites are heterozygoues in parental genomes. 
        ## delta_het_parents
        p_homo_given_ZAA = (self.params['delta_het_parents'] * p_homo_given_gaa ) + ((1 - self.params['delta_het_parents']) * p_homo_given_gab )
        p_homo_given_ZAB = ((1 - self.params['avg_sites_segregating']) * p_homo_given_gaa ) + ( self.params['avg_sites_segregating'] * p_homo_given_gab )

        emission_prob = [ 
            [p_homo_given_ZAA,  1 - p_homo_given_ZAA, 0.5],
            [p_homo_given_ZAB, 1 - p_homo_given_ZAB, 0.5]
        ]
        print("transmission prob:\n %s" % pd.DataFrame(transition_prob) )
        print("emission prob:\n %s" % pd.DataFrame(emission_prob) )
        model = hmm.MultinomialHMM(n_components=2, n_iter = self.params['n_iter'], algorithm = "viterbi", init_params='st') 
        model.startprob_ = np.array(init_prob)
        model.transmat_ = pd.DataFrame(transition_prob)
        model.emissionprob_ = pd.DataFrame(emission_prob)
        return(model)

    def _input_obs( self, sample_gt ):
        t_obs = np.array(sample_gt)
        t_obs[t_obs == 1] = 0
        t_obs[t_obs == 2] = 1
        t_obs[t_obs == -1] = 2
        return(t_obs)

    def hmm_decode( self ):
        self._obs = self._input_obs( self.sample_gt )
        t_genotypes = self.model.decode( self._obs.reshape((-1, 1)) )
        return(t_genotypes)
