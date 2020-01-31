import pandas as pd
import numpy as np
import allel
import snpmatch
import logging
import os
import json
import re

log = logging.getLogger(__name__)

def parseGT(snpGT):
    snpBinary = np.zeros(len(snpGT), dtype = "int8")
    if len(snpBinary) == 0:
        return(snpBinary)
    first = snpGT[0]
    if first.find('|') != -1:
        ## GT is phased
        separator = "|"
    elif first.find('/') != -1:
        ## GT is not phased
        separator = "/"
    elif np.char.isdigit(first):
        return np.array(np.copy(snpGT), dtype = "int8")
    else:
        snpmatch.die("unable to parse the format of GT in vcf!")
    hetGT_1 = "0" + separator + "1"
    hetGT_2 = "1" + separator + "0"
    refGT = "0" + separator + "0"
    altGT = "1" + separator + "1"
    nocall = "." + separator + "."
    snpBinary[np.where(snpGT == altGT)[0]] = 1
    snpBinary[np.where((snpGT == hetGT_1) | (snpGT == hetGT_2) )[0]] = 2
    snpBinary[np.where(snpGT == nocall)[0]] = -1
    return(snpBinary)

def snp_binary_to_gt(snpBinary):
    snpBinary = np.array(snpBinary, dtype="int8")
    snpGT = np.zeros(len(snpBinary), dtype="S8")
    snpGT[np.where(snpBinary == -1)[0]] = "./."
    snpGT[np.where(snpBinary == 0)[0]] = "0/0"
    snpGT[np.where(snpBinary == 1)[0]] = "1/1"
    snpGT[np.where(snpBinary == 2)[0]] = "0/1"
    return(snpGT)

class ParseInputs(object):
    ## class object for parsing input files for SNPmatch

    def __init__(self, inFile, logDebug=True, outFile = "parser" ):
        if outFile == "parser" or not outFile:
            outFile = inFile + ".snpmatch"
        if os.path.isfile(inFile + ".snpmatch.npz"):
            log.info("snpmatch parser dump found! loading %s", inFile + ".snpmatch.npz")
            snps = np.load(inFile + ".snpmatch.npz")
            self.load_snp_info(snps['chr'], snps['pos'], snps['gt'], snps['wei'], snps['dp'])
            log.info("done!")
        elif os.path.isfile(inFile):
            _,inType = os.path.splitext(inFile)
            if inType == '.npz':
                log.info("loading snpmatch parser file! %s", inFile)
                snps = np.load(inFile)
                self.load_snp_info(snps['chr'], snps['pos'], snps['gt'], snps['wei'], snps['dp'])
            else:
                log.info('running snpmatch parser!')
                if inType == '.vcf' or len(re.compile(".vcf.gz$").findall(os.path.basename(inFile))) > 0 :
                    (snpCHR, snpPOS, snpGT, snpWEI, DPmean) = self.read_vcf(inFile, logDebug)
                elif inType == '.bed':
                    (snpCHR, snpPOS, snpGT, snpWEI, DPmean) = self.read_bed(inFile, logDebug)
                else:
                    snpmatch.die("input file type %s not supported" % inType)
                self.load_snp_info(snpCHR, snpPOS, snpGT, snpWEI, DPmean)
                self.save_snp_info(outFile)
                self.case_interpret_inputs(outFile + ".stats.json")
            log.info("done!")

    def load_snp_info(self, snpCHR, snpPOS, snpGT, snpWEI, DPmean):
        self.chrs = np.array(snpCHR, dtype="string")
        self.pos = np.array(snpPOS, dtype=int)
        self.gt = np.array(snpGT, dtype="string")
        self.wei = np.array(snpWEI, dtype=float)
        self.dp = DPmean

    def save_snp_info(self, outFile):
        log.info("creating snpmatch parser file: %s", outFile + '.npz')
        np.savez(outFile, chr = self.chrs, pos = self.pos, gt = self.gt, wei = self.wei, dp = self.dp)

    def case_interpret_inputs(self, outFile):
        NumSNPs = len(self.chrs)
        case = 0
        note = "Sufficient number of SNPs"
        if NumSNPs < snpmatch.snp_thres:
            note = "Attention: low number of SNPs provided"
            case = 1
        snpst = np.unique(self.chrs, return_counts=True)
        snpdict = dict(('%s' % snpst[0][i], snpst[1][i]) for i in range(len(snpst[0])))
        statdict = {"interpretation": {"case": case, "text": note}, "snps": snpdict, "num_of_snps": NumSNPs, "depth": self.dp}
        statdict['percent_heterozygosity'] = snpmatch.getHeterozygosity(self.gt)
        with open(outFile , "w") as out_stats:
            out_stats.write(json.dumps(statdict))

    @staticmethod
    def read_bed(inFile, logDebug):
        log.info("reading the position file")
        targetSNPs = pd.read_table(inFile, header=None, usecols=[0,1,2])
        snpCHR = np.array(targetSNPs[0], dtype="string")
        snpPOS = np.array(targetSNPs[1], dtype=int)
        snpGT = np.array(targetSNPs[2])
        snpBinary = parseGT(snpGT)
        snpWEI = np.ones((len(snpCHR), 3))  ## for homo and het
        snpWEI[np.where(snpBinary != 0),0] = 0
        snpWEI[np.where(snpBinary != 1),2] = 0
        snpWEI[np.where(snpBinary != 2),1] = 0
        return((snpCHR, snpPOS, snpGT, snpWEI, "NA"))

    @staticmethod
    def get_wei_from_GT(snpGT):
        snpBinary = parseGT(snpGT)
        snpWEI = np.ones((len(snpGT), 3))  ## for homo and het
        snpWEI[np.where(snpBinary != 0),0] = 0
        snpWEI[np.where(snpBinary != 1),2] = 0
        snpWEI[np.where(snpBinary != 2),1] = 0
        return(snpWEI)

    def read_vcf(self, inFile, logDebug):
        if logDebug:
            vcf = allel.read_vcf(inFile, samples = [0], fields = '*')
        else:
            import StringIO
            import sys
            sys.stderr = StringIO.StringIO()
            vcf = allel.read_vcf(inFile, samples = [0], fields = '*')
            sys.stderr = sys.__stderr__
        try:
            snpGT = allel.GenotypeArray(vcf['calldata/GT']).to_gt()[:, 0]
        except AttributeError:
            snpmatch.die("input VCF file doesnt have required GT field")
        snpsREQ = np.where((snpGT != './.') & (snpGT != '.|.'))[0]
        snpGT = snpGT[snpsREQ]
        if 'calldata/PL' in sorted(vcf.keys()):
            snpWEI = np.copy(vcf['calldata/PL'][snpsREQ, 0]).astype('float')
            missing_pls = np.all(snpWEI == -1, axis = 1)
            snpWEI = snpWEI/(-10)
            snpWEI = np.exp(snpWEI)
            snpWEI[missing_pls,] = self.get_wei_from_GT(snpGT[missing_pls])
        else:
            snpWEI = self.get_wei_from_GT(snpGT)
        snpCHR = np.array(vcf['variants/CHROM'][snpsREQ], dtype="string")
        if 'calldata/DP' in sorted(vcf.keys()):
            DPmean = np.mean(vcf['calldata/DP'][snpsREQ,0])
        else:
            DPmean = "NA"
        snpPOS = np.array(vcf['variants/POS'][snpsREQ])
        return((snpCHR, snpPOS, snpGT, snpWEI, DPmean))

    def filter_chr_names(self):
        ## provide genotypedata (pygwas genotype object)
        self.g_chrs = np.array(pd.Series(self.chrs).str.replace("chr", "", case=False), dtype="string")
        self.g_chrs_ids = np.unique(self.g_chrs)


def potatoParser(inFile, logDebug, outFile = "parser"):
    inputs = ParseInputs(inFile, logDebug, outFile)
    return(inputs.chrs, inputs.pos, inputs.gt, inputs.wei, inputs.dp)
