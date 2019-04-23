import numpy as np
import numpy.ma
import scipy.stats as st
from pygwas.core import genotype
import pandas as pd
import logging
import os
import re
from glob import glob
import itertools

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
        snp_idx = np.zeros(0, dtype=int)
        for i in enumerate(self.g.chrs):
            perchrTar_POS = commonSNPsPOS[np.where(commonSNPsCHR == i[1])[0]]
            perchrAcc_POS = self.g.positions[self.g.chr_regions[i[0]][0]:self.g.chr_regions[i[0]][1]]
            perchrAcc_ix = np.where( np.in1d(perchrAcc_POS, perchrTar_POS) )[0] + self.g.chr_regions[i[0]][0]
            snp_idx = np.append( snp_idx, perchrAcc_ix )
        return(snp_idx)
