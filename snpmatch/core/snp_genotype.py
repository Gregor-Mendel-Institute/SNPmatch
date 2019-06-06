import numpy as np
from pygwas.core import genotype
import pandas as pd
import logging
import re
from glob import glob
from . import parsers
from . import genomes

log = logging.getLogger(__name__)


## Class object adapted from PyGWAS genotype object
class Genotype(object):

    def __init__(self, hdf5_file, hdf5_acc_file=None):
        self.g = genotype.load_hdf5_genotype_data(hdf5_file)
        if hdf5_acc_file is None:
            hdf5_acc_file = re.sub('\.hdf5$', '', hdf5_file) + '.acc.hdf5'
            if len(glob(hdf5_acc_file)) > 0:
                self.g_acc = genotype.load_hdf5_genotype_data(hdf5_acc_file)
        else:
            self.g_acc = genotype.load_hdf5_genotype_data(hdf5_acc_file)

    def get_positions_idxs(self, commonSNPsCHR, commonSNPsPOS):
        common_ins = parsers.ParseInputs("")
        common_ins.load_snp_info( snpCHR=commonSNPsCHR, snpPOS=commonSNPsPOS, snpGT="", snpWEI=np.nan, DPmean=0 )
        common_ins.filter_chr_names()
        g_chrs_ids = np.char.replace(np.core.defchararray.lower(np.array(self.g.chrs, dtype="string")), "chr", "")
        common_chr_ids = np.intersect1d(g_chrs_ids, common_ins.g_chrs_ids)
        acc_idx = np.zeros(0, dtype=int)
        tar_idx = np.zeros(0, dtype=int)
        for i in common_chr_ids:
            perchrTar = np.where(common_ins.g_chrs == i)[0]
            t_g_chr_ix = np.where(g_chrs_ids == i)[0][0]
            perchrTar_POS = np.array(common_ins.pos[perchrTar], dtype=int)
            perchrAcc_POS = self.g.positions[self.g.chr_regions[t_g_chr_ix][0]:self.g.chr_regions[t_g_chr_ix][1]]
            perchrAcc_ix = np.where( np.in1d(perchrAcc_POS, perchrTar_POS, assume_unique=True) )[0] + self.g.chr_regions[t_g_chr_ix][0]
            perchrTar_ix = perchrTar[np.where( np.in1d(perchrTar_POS, perchrAcc_POS, assume_unique=True) )[0]]
            acc_idx = np.append( acc_idx, perchrAcc_ix )
            tar_idx = np.append( tar_idx, perchrTar_ix )
        return((acc_idx, tar_idx))
