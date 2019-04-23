import numpy as np
from pygwas.core import genotype
import pandas as pd
import logging
import re
from glob import glob
from . import parsers

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
        common_ins.load_snp_info( snpCHR=commonSNPsCHR, snpPOS=commonSNPsPOS, snpGT="", snpWEI="", DPmean="" )
        common_ins.filter_chr_names(self.g)
        acc_idx = np.zeros(0, dtype=int)
        tar_idx = np.zeros(0, dtype=int)
        for i in enumerate(common_ins.chr_list):
            perchrTar = np.where(common_ins.chrs_nochr == i[1])[0]
            perchrTar_POS = common_ins.pos[perchrTar]
            perchrAcc_POS = self.g.positions[self.g.chr_regions[i[0]][0]:self.g.chr_regions[i[0]][1]]
            perchrAcc_ix = np.where( np.in1d(perchrAcc_POS, perchrTar_POS) )[0] + self.g.chr_regions[i[0]][0]
            perchrTar_ix = perchrTar[np.where( np.in1d(perchrTar_POS, perchrAcc_POS) )[0]]
            acc_idx = np.append( acc_idx, perchrAcc_ix )
            tar_idx = np.append( tar_idx, perchrTar_ix )
        return((acc_idx, tar_idx))
