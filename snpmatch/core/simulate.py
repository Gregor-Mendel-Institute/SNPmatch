import logging
import numpy as np
import pandas as pd
from . import snpmatch
from . import parsers
from . import snp_genotype

log = logging.getLogger(__name__)

def simulateSNPs(g, input_df, numSNPs, outFile, err_rate):
    ## Input -- pandas dataframe with chr, position and genotype
    assert type(input_df) == pd.core.frame.DataFrame, "please provide a pandas dataframe"
    assert input_df.shape[1] >= 3, "first three columns are needed in dataframe: chr, pos, snp"
    ## default error rates = 0.001
    log.info("sampling %s positions" % numSNPs)
    sampleSNPs = np.sort(np.random.choice(np.arange(input_df.shape[0]), numSNPs, replace=False))
    input_df = input_df.iloc[sampleSNPs,:]
    log.info("adding in error rate: %s" % err_rate)
    num_to_change = int(err_rate * input_df.shape[0])
    input_df.iloc[np.sort(np.random.choice(np.arange(input_df.shape[0]), num_to_change, replace=False)), 2] = np.random.choice(3, num_to_change)
    inputs = parsers.ParseInputs(inFile="")
    inputs.load_snp_info(np.array(input_df.iloc[:,0]), np.array(input_df.iloc[:,1]), parsers.snp_binary_to_gt( np.array(input_df.iloc[:,2]) ), inputs.get_wei_from_GT( parsers.snp_binary_to_gt( np.array(input_df.iloc[:,2]) ) ), "NA")
    snpmatch_result = snpmatch.genotyper(inputs, g, outFile)


def potatoSimulate(args):
    g = snp_genotype.Genotype(args['hdf5File'], args['hdf5accFile'] )
    try:
        AccToCheck = np.where(g.g.accessions == args['AccID'])[0][0]
    except NameError:
        snpmatch.die("accession is not present in the matrix!")
    log.info("loading input files")
    acc_snp = g.g_acc.snps[:,AccToCheck]
    informative_snps = np.where(acc_snp >= 0)[0] ## Removing NAs for accession
    input_df = pd.DataFrame(np.column_stack((np.array(g.g.chromosomes)[informative_snps], g.g.positions[informative_snps], acc_snp[informative_snps] )), columns = ["chr", 'pos', 'snp'])
    log.info("done!")
    simulateSNPs(g, input_df, args['numSNPs'], args['outFile'], args['err_rate'])
    log.info("finished!")
