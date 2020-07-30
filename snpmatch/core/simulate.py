import logging
import numpy as np
import pandas as pd
from . import snpmatch
from . import parsers
from . import snp_genotype

log = logging.getLogger(__name__)

def simulateSNPs(g, AccID, numSNPs, outFile=None, err_rate=0.001):
    assert type(AccID) is str, "provide Accession ID as a string"
    assert AccID in g.g.accessions, "accession is not present in the matrix!"
    AccToCheck = np.where(g.g.accessions == AccID)[0][0]
    log.info("loading input files")
    acc_snp = g.g_acc.snps[:,AccToCheck]
    informative_snps = np.where(acc_snp >= 0)[0] ## Removing NAs for accession
    input_df = pd.DataFrame(np.column_stack((np.array(g.g.chromosomes)[informative_snps], g.g.positions[informative_snps], acc_snp[informative_snps] )), columns = ["chr", 'pos', 'snp'])
    ## Input -- pandas dataframe with chr, position and genotype
    #assert type(input_df) == pd.core.frame.DataFrame, "please provide a pandas dataframe"
    #assert input_df.shape[1] >= 3, "first three columns are needed in dataframe: chr, pos, snp"
    ## default error rates = 0.001
    log.info("sampling %s positions" % numSNPs)
    sampleSNPs = np.sort(np.random.choice(np.arange(input_df.shape[0]), numSNPs, replace=False))
    input_df = input_df.iloc[sampleSNPs,:]
    log.info("adding in error rate: %s" % err_rate)
    num_to_change = int(err_rate * input_df.shape[0])
    input_df.iloc[np.sort(np.random.choice(np.arange(input_df.shape[0]), num_to_change, replace=False)), 2] = np.random.choice(3, num_to_change)
    input_df.iloc[:, 2] = parsers.snp_binary_to_gt( np.array(input_df.iloc[:,2]) )
    if outFile is not None:
        input_df.to_csv( outFile, sep = "\t", index = None, header = False  )
    return(input_df)

def simulateSNPs_F1(g, parents, numSNPs, outFile = None, err_rate, rm_hets = 1):
    indP1 = np.where(g.g_acc.accessions == parents.split("x")[0])[0][0]
    indP2 = np.where(g.g_acc.accessions == parents.split("x")[1])[0][0]
    log.info("loading files!")
    snpsP1 = g.g_acc.snps[:,indP1]
    snpsP2 = g.g_acc.snps[:,indP2]
    common_ix = np.where((snpsP1 >= 0) & (snpsP2 >= 0) & (snpsP1 < 2) & (snpsP2 < 2))[0]
    segregating_ix = np.where(snpsP1[common_ix] != snpsP2[common_ix] )[0]
    diff_ix = np.setdiff1d( np.arange(len(common_ix)), segregating_ix )
    common_snps = np.zeros(len(common_ix), dtype="int8")
    common_snps[segregating_ix] = 2
    common_snps[diff_ix] = snpsP1[common_ix[diff_ix]]
    input_df = pd.DataFrame( np.column_stack((np.array(g.g_acc.chromosomes)[common_ix], np.array(g.g_acc.positions)[common_ix], common_snps  )), columns = ["chr", 'pos', 'snp'] )
    log.info("sampling %s positions" % numSNPs)
    sampleSNPs = np.sort(np.random.choice(np.arange(input_df.shape[0]), numSNPs, replace=False))
    input_df = input_df.iloc[sampleSNPs,:]
    input_df['snp'] = input_df['snp'].astype(int)
    log.info("adding in error rate: %s" %  err_rate)
    num_to_change = int(err_rate * input_df.shape[0])
    input_df.iloc[np.sort(np.random.choice(np.where(input_df['snp'] != 2)[0], num_to_change, replace=False)), 2] = np.random.choice(2, num_to_change)
    ## Also change hets randomly to homozygous
    het_ix = np.where(input_df['snp'] == 2)[0]
    input_df.iloc[het_ix, 2] = np.random.choice(3, het_ix.shape[0], p=[(1-rm_hets)/2,(1-rm_hets)/2,rm_hets])
    ## Save the file to a bed file
    input_df.iloc[:, 2] = parsers.snp_binary_to_gt( np.array(input_df.iloc[:,2]) )
    if outFile is not None:
        input_df.to_csv( outFile, sep = "\t", index = None, header = False  )
    return(input_df)

def potatoSimulate(args):
    g = snp_genotype.Genotype(args['hdf5File'], args['hdf5accFile'] )
    if args['simF1']:
        simulateSNPs_F1(g, args['AccID'], args['numSNPs'], args['outFile'], args['err_rate'], args['rm_het'])
    else:
        simulateSNPs(g, args['AccID'], args['numSNPs'], args['outFile'], args['err_rate'])
    log.info("finished!")
