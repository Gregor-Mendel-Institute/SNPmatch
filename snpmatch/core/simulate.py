import math
import random
import logging
import numpy as np
import numpy.ma
import h5py
from pygwas.core import genotype
import scipy
from . import snpmatch
from . import snp_genotype

log = logging.getLogger(__name__)


def simulateSNPs(g, AccID, numSNPs, outFile, err_rate, chunk_size = 1000):
    ## default error rates = 0.001
    try:
        AccToCheck = np.where(g.g.accessions == AccID)[0][0]
    except NameError:
        snpmatch.die("accession is not present in the matrix!")
    tot_snps = np.where(g.g_acc.snps[:,AccToCheck] >= 0)[0]
    log.info("randomly choosing %s SNPs from %s SNPs in accession %s", numSNPs, len(tot_snps), AccID)
    sampleSNPs = np.sort(random.sample(tot_snps, numSNPs))
    log.info("done!")
    ScoreList = np.zeros(len(g.g.accessions), dtype="uint32")
    NumInfoSites = np.zeros(len(g.g.accessions), dtype="uint32")
    NumMatSNPs = 0
    num_lines = len(g.g.accessions)
    for i in range(0, len(sampleSNPs), chunk_size):
        t1001SNPs = g.g.snps[tuple(sampleSNPs[i:i+chunk_size]),]
        sam1kSNPs = t1001SNPs[:, AccToCheck]
        samSNPs = np.reshape(np.repeat(sam1kSNPs, num_lines), (len(sam1kSNPs), num_lines))
        ScoreList = ScoreList + np.sum(t1001SNPs == samSNPs, axis=0)
        NumInfoSites = NumInfoSites + len(sampleSNPs[i:i+chunk_size]) - np.sum(numpy.ma.masked_less(t1001SNPs, 0).mask, axis = 0)
        NumMatSNPs = NumMatSNPs + len(sam1kSNPs)
        if i/chunk_size % 50 == 0:
            log.info("Done analysing %s SNPs", NumMatSNPs)
    logging.info("writing data!")
    ScoreList = np.array([random.uniform(ScoreList[i] - (ScoreList[i] * err_rate), ScoreList[i]) for i in range(num_lines)], dtype=int)
    snpmatch.print_out_table(outFile, g.g.accessions, ScoreList, NumInfoSites, numSNPs, "NA")

def potatoSimulate(args):
    log.info("loading database files")
    g = snp_genotype.Genotype(args['hdf5File'], args['hdf5accFile'] )
    log.info("done!")
    simulateSNPs(g.g, g.g_acc, args['AccID'], args['numSNPs'], args['outFile'], args['err_rate'])
    log.info("finished!")
