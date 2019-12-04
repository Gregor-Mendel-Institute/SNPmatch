import numpy as np
from pygwas.core import genotype
import pandas as pd
import logging
import re
from glob import glob
from . import parsers
from . import genomes
import allel
import itertools
import os.path
import numbers

log = logging.getLogger(__name__)

chunk_size = 1000
def load_genotype_files(h5file, hdf5_acc_file=None):
    return(Genotype(h5file, hdf5_acc_file))

## Class object adapted from PyGWAS genotype object
class Genotype(object):

    def __init__(self, hdf5_file, hdf5_acc_file):
        assert hdf5_file is not None or hdf5_acc_file is not None, "Provide atleast one hdf5 genotype file"
        if hdf5_file is None:
            assert os.path.isfile(hdf5_acc_file), "Path to %s seems to be broken" % hdf5_acc_file
            self.g_acc = genotype.load_hdf5_genotype_data(hdf5_acc_file)
            return(None)
        assert os.path.isfile(hdf5_file), "Path to %s seems to be broken" % hdf5_file
        self.g = genotype.load_hdf5_genotype_data(hdf5_file)
        if hdf5_acc_file is None:
            hdf5_acc_file = re.sub('\.hdf5$', '', hdf5_file) + '.acc.hdf5'
            if len(glob(hdf5_acc_file)) > 0:
                self.g_acc = genotype.load_hdf5_genotype_data(hdf5_acc_file)
        else:
            self.g_acc = genotype.load_hdf5_genotype_data(hdf5_acc_file)

    def get_positions_idxs(self, commonSNPsCHR, commonSNPsPOS):
        return(self.get_common_positions( np.array(self.g.chromosomes, dtype="string"), self.g.positions, commonSNPsCHR, commonSNPsPOS ))

    @staticmethod
    def get_common_positions(input_1_chr, input_1_pos, input_2_chr, input_2_pos):
        assert len(input_1_chr) == len(input_1_pos), "Both chromosome and position array provided should be of same length"
        assert len(input_2_chr) == len(input_2_pos), "Both chromosome and position array provided should be of same length"
        common_ins_1 = parsers.ParseInputs("")
        common_ins_1.load_snp_info( snpCHR=input_1_chr, snpPOS=input_1_pos, snpGT="", snpWEI=np.nan, DPmean=0 )
        common_ins_1.filter_chr_names()
        common_ins_2 = parsers.ParseInputs("")
        common_ins_2.load_snp_info(  snpCHR=input_2_chr, snpPOS=input_2_pos, snpGT="", snpWEI=np.nan, DPmean=0 )
        common_ins_2.filter_chr_names()
        common_chr_ids = np.intersect1d(common_ins_1.g_chrs_ids, common_ins_2.g_chrs_ids)
        common_idx_1 = np.zeros(0, dtype=int)
        common_idx_2 = np.zeros(0, dtype=int)
        for i in common_chr_ids:
            perchrTar_ix_1 = np.where(common_ins_1.g_chrs == i)[0]
            perchrTar_ix_2 = np.where(common_ins_2.g_chrs == i)[0]
            perchrTar_pos_1 = np.array(common_ins_1.pos[perchrTar_ix_1], dtype=int)
            perchrTar_pos_2 = np.array(common_ins_2.pos[perchrTar_ix_2], dtype=int)
            common_idx_1 = np.append( common_idx_1, perchrTar_ix_1[np.where( np.in1d(perchrTar_pos_1, perchrTar_pos_2, assume_unique=True) )[0]] )
            common_idx_2 = np.append( common_idx_2, perchrTar_ix_2[np.where( np.in1d(perchrTar_pos_2, perchrTar_pos_1, assume_unique=True) )[0]] )
        return((common_idx_1, common_idx_2))

    def get_matching_accs_ix(self, accs, return_np=False):
        acc_ix = []
        for ea in accs:
            t_ix = np.where(self.g.accessions == ea)[0]
            if len(t_ix) == 0:
                acc_ix.append(None)
            else:
                acc_ix.append(t_ix[0])
        if return_np:
            acc_ix = np.array(acc_ix)[np.where(np.not_equal(acc_ix, None))[0]].astype("int")
        return(acc_ix)

    def load_snps_given_accs_pos(self, out_file, accs_ix = None, pos_ix=None, ):
        ## Returnning snp matrix from a indices.
        if os.path.isfile(out_file + ".npz"):
            return(np.load(out_file + ".npz"))
        if accs_ix is not None and pos_ix is None:
            req_snps = np.zeros((0, len(accs_ix)), dtype = type(self.g.snps[0,0]))
            pos_ix = np.zeros(0, dtype =int)
            for t_ix in range(0, self.g.snps.shape[0], chunk_size):
                t_s = self.g.snps[t_ix:t_ix+chunk_size,][:,accs_ix]
                t_filter = np.where(~np.all(t_s == -1, axis = 1))[0]  ## Here -1 is no information
                req_snps = np.vstack((req_snps, t_s[t_filter,]))
                pos_ix = np.append(pos_ix, t_filter + t_ix)
        elif pos_ix is not None and accs_ix is None:
            req_snps = self.g.snps[pos_ix,:]
            accs_ix = np.arange(self.g.accessions.shape[0])
        elif accs_ix is not None and pos_ix is not None:
            req_snps = self.g.snps[pos_ix,:][:,accs_ix]
        else:
            log.warn("either provide accs_ix or pos_ix!")
            req_snps = np.zeros((0,0), dtype = type(self.g.snps[0,0]))
        np.savez(out_file, req_snps = req_snps, accs_ix=accs_ix, pos_ix=pos_ix)
        return(np.load(out_file + ".npz"))

    def get_polarized_snps(self, polarize_geno = 1, return_float=True):
        if hasattr(self, "req_snps"):
            self.pol_snps = _polarize_snps( np.array(self.req_snps), polarize_geno)
        if  "pol_snps" not in self.__dict__.keys():
            if self.g.snps.shape[0]/chunk_size > 2000: ## It is a huge array
                self.pol_snps = np.zeros((0,self.g.snps.shape[1]), dtype="int8")
                for t_ix in range(0, self.g.snps.shape[0], chunk_size):
                    t_s = _polarize_snps( self.g.snps[t_ix:t_ix+chunk_size,:], polarize_geno )
                    self.pol_snps = np.vstack((self.pol_snps, t_s))
            else:
                self.pol_snps = _polarize_snps( np.array(self.g.snps), polarize_geno)
        if return_float:
            self.pol_snps_fl = self.np_snp_to_pd_df(self.pol_snps, drop_na_all=True)

    def get_af_snps(self, no_accs_missing_info, polarize_geno = 1):
        maf_snps = af_snps_np_array(self.g.snps, no_accs_missing_info, polarize_geno)
        return(maf_snps)

    @staticmethod
    def np_snp_to_pd_df(np_arr, drop_na_all = True):
        np_arr = np.array(np_arr, dtype = float)
        np_arr[np_arr == -1] = np.nan
        np_arr[np_arr == 2] = 0.5
        pd_df = pd.DataFrame(np_arr)
        if drop_na_all:
            return(pd_df.dropna(how = "all"))
        else:
            return(pd_df.dropna(how = "any"))

    def identify_segregating_snps(self, accs_ix ):
        assert type(accs_ix) is np.ndarray, "provide an np array for list of indices to be considered"
        assert len(accs_ix) > 1, "polymorphism happens in more than 1 line"
        if len(accs_ix) > (len(self.g.accessions) / 2):
            return( None )
        if len(accs_ix) < 10:
            t_snps = np.zeros(( self.g_acc.snps.shape[0], len(accs_ix) ))
            for ef in range(len(accs_ix)):
                t_snps[:,ef] = self.g_acc.snps[:,accs_ix[ef]]
            seg_counts = segregting_snps(t_snps)
            div_counts = np.divide(seg_counts[0], seg_counts[1], where = seg_counts[1] != 0 )
            seg_ix = np.setdiff1d(np.where(div_counts  < 1 )[0], np.where(seg_counts[1] == 0)[0])
            return( seg_ix )
        NumSNPs = self.g.positions.shape[0]
        seg_counts = np.zeros(0, dtype=int)
        total_counts = np.zeros(0, dtype=int)
        for j in range(0, NumSNPs, chunk_size):
            t1001SNPs = np.array(self.g.snps[j:j+chunk_size,:][:,accs_ix], dtype=float)
            t1001SNPs = segregting_snps( t1001SNPs )
            seg_counts = np.append(seg_counts, t1001SNPs[0] )
            total_counts = np.append(total_counts, t1001SNPs[1] )
        div_counts = np.divide(seg_counts, total_counts, where = total_counts != 0 )
        seg_ix = np.setdiff1d(np.where(div_counts  < 1 )[0], np.where(total_counts == 0)[0])
        return( seg_ix )

def segregting_snps(t):
    t[t < 0] = np.nan
    t = np.sort(t,axis=1)
    t_r_sum = np.sum( ~np.isnan(t), axis = 1)
    t_sum = np.nansum(t[:,1:] == t[:,:-1], axis=1) + 1
    return((t_sum, t_r_sum))

def _polarize_snps(snps, polarize_geno=1, genotypes=[0, 1]):
    assert len(genotypes) == 2, "assuming it is biallelic"
    t_s =  np.array(snps)
    t_int_ix = np.where(np.sum(t_s == polarize_geno, axis = 1) > float(snps.shape[1])/2)[0]
    t_replace = t_s[t_int_ix, :]
    t_replace[t_replace == genotypes[1]] = 3
    t_replace[t_replace == genotypes[0]] = genotypes[1]
    t_replace[t_replace == 3] = genotypes[0]
    t_s[t_int_ix, :] = t_replace
    return(t_s)

def af_snps_np_array(snps, no_accs_missing_info, polarize_geno=1):
    maf_snps = np.zeros(0, dtype="float")
    for t_ix in range(0, snps.shape[0], chunk_size):
        t_s =  snps[t_ix:t_ix+chunk_size,:]
        t_t = snps.shape[1] - np.sum(t_s == -1, axis = 1)
        t_no_1 = 2 * np.sum(t_s == polarize_geno, axis = 1) + np.sum(t_s == 2, axis = 1)
        t_maf = np.array(t_no_1, dtype=float) / np.multiply(2, t_t)
        t_maf[ np.where(t_t <= no_accs_missing_info)[0] ] = np.nan
        maf_snps = np.append(maf_snps, t_maf )
    return(maf_snps)

def get_sq_diversity_np(snps, acc_ix=None):
    assert type(snps) is pd.core.frame.DataFrame, "please provide pd.Dataframe as input"
    if isinstance(acc_ix, numbers.Integral):
        assert acc_ix < snps.shape[1], "index of a reference to get sq diversity for all the other"
        kin_mat = np.zeros(snps.shape[1], dtype=float)
        for i in range(snps.shape[1]):
            if i == acc_ix:
                kin_mat[acc_ix] = 0
            else:
                t_s = snps.iloc[:,[acc_ix,i]]
                kin_mat[i] = allel.sequence_diversity(range(snps.shape[0]), allel.AlleleCountsArray(np.column_stack((np.sum(t_s == 0, axis =1) * 2, np.sum(t_s == 0.5, axis =1) * 2, np.sum(t_s == 1, axis =1) * 2))))
        return(kin_mat)
    if acc_ix is None:
        acc_ix = np.arange(snps.shape[1])
    assert type(acc_ix) is np.ndarray, "provide an index for samples to get pairwise scores"
    kin_mat = pd.DataFrame(0, index = acc_ix, columns = acc_ix, dtype = float)
    for i,j in itertools.combinations(acc_ix, 2):
        t_k = allel.sequence_diversity(range(snps.shape[0]), allel.AlleleCountsArray(np.column_stack((np.sum(snps.iloc[:,[i,j]] == 0, axis =1) * 2, np.sum(snps.iloc[:,[i,j]] == 0.5, axis =1) * 2, np.sum(snps.iloc[:,[i,j]] == 1, axis =1) * 2))))
        #t_k = np.sum(snps.iloc[:,i] == snps.iloc[:,j])/float(snps.shape[0])
        kin_mat.loc[i,j] = t_k
        kin_mat.loc[j,i] = t_k
    return(kin_mat)

def snpmat_character_to_biallellic(snpmat, polarize = True):
    assert type(snpmat) is pd.DataFrame, "please provide a pd dataframe, ideally a small array"
    snpmat_num = pd.DataFrame(dtype="int", index = snpmat.index, columns = snpmat.columns )
    genotypes=["A", "T", "G", "C"]
    snpmat_num[snpmat == genotypes[0]] = 0
    snpmat_num[snpmat == genotypes[1]] = 1
    snpmat_num[snpmat == genotypes[2]] = 2
    snpmat_num[snpmat == genotypes[3]] = 3
    snpmat_num[~snpmat.isin(genotypes)] = -1
    # snpmat_num[snpmat == "N"] = -1
    for index, row in snpmat_num.iterrows():
        t_r = row.factorize(sort=True)
        t_r[0][t_r[0] == 0] = -1
        t_r[0][t_r[0] == 1] = 0
        t_r[0][t_r[0] == 2] = 1
        # t_r[0][t_r[0] == 3] = 2
        snpmat_num.loc[index] = t_r[0]
    if polarize:
        return(_polarize_snps(snpmat_num, polarize_geno=1, genotypes=[0, 1]))
    return(np.array(snpmat_num))
