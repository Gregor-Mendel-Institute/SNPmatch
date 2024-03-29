import pandas as pd
import numpy as np
import allel
from . import snpmatch
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


def parse_sample_names(sample_files, file_sep = "_"):
    sample_files = pd.Series(sample_files).astype(str)
    input_ids = sample_files.apply( os.path.basename ).str.split(file_sep, expand=True)
    if np.unique(input_ids.iloc[:,0]).shape[0] == input_ids.shape[0]:
        input_ids = np.array(input_ids.iloc[:,0])
    elif np.unique(input_ids.iloc[:,0] + file_sep + input_ids.iloc[:,1]).shape[0] == input_ids.shape[0]:
        input_ids = np.array(input_ids.iloc[:,0] + file_sep + input_ids.iloc[:,1])
    else:
        input_ids = input_files.apply( os.path.basename ).str.replace( ".scores.txt", "" ).values
    return(input_ids)
    

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
        self.chrs = np.array(snpCHR, dtype="str")
        self.pos = np.array(snpPOS, dtype=int)
        self.gt = np.array(snpGT, dtype="str")
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
        snpdict = dict(('%s' % snpst[0][i], int(snpst[1][i])) for i in range(len(snpst[0])))
        statdict = {}
        statdict["snps"] = snpdict
        statdict["interpretation"] = {"case": case, "text": note}
        statdict["num_of_snps"] = NumSNPs
        statdict["depth"] = np.nanmean(self.dp)
        statdict['percent_heterozygosity'] = snpmatch.getHeterozygosity(self.gt)
        with open(outFile , "w") as out_stats:
            out_stats.write(json.dumps(statdict))

    @staticmethod
    def read_bed(inFile, logDebug):
        log.info("reading the position file")
        targetSNPs = pd.read_csv(inFile, header=None, sep = None, engine = 'python', usecols=[0,1,2])
        snpCHR = np.array(targetSNPs[0], dtype="str")
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
        snp_inputs = import_vcf_file( inFile, logDebug, samples_to_load = [0] )
        snp_inputs['gt'] = snp_inputs['gt'][:,0]
        snpsREQ = np.where((snp_inputs['gt'] != './.') & (snp_inputs['gt'] != '.|.'))[0]
        snpGT = snp_inputs['gt'][snpsREQ]
        if 'wei' in sorted(snp_inputs.keys()):
            snpWEI = snp_inputs['wei'][snpsREQ, 0]
            missing_pls = np.all(snpWEI == -1, axis = 1)
            snpWEI = snpWEI/(-10)
            snpWEI = np.exp(snpWEI)
            snpWEI[missing_pls,] = self.get_wei_from_GT(snpGT[missing_pls])
        else:
            snpWEI = self.get_wei_from_GT(snpGT)
        snpCHR = snp_inputs['chr'][snpsREQ]
        snpPOS = snp_inputs['pos'][snpsREQ]
        snpDP = snp_inputs['dp'][snpsREQ]
        return((snpCHR, snpPOS, snpGT, snpWEI, snpDP))

    def filter_chr_names(self):
        ## provide genotypedata (pygwas genotype object)
        self.g_chrs = np.array(pd.Series(self.chrs).str.replace("chr", "", case=False), dtype="str")
        _, idx = np.unique(self.g_chrs, return_index=True)
        self.g_chrs_ids = self.g_chrs[np.sort(idx)]

    def save_to_bed(self, outFile):
        # parsers.snp_binary_to_gt( np.array(input_df.iloc[:,2]) )
        input_df = pd.DataFrame(
            np.column_stack((
                self.chrs,
                self.pos, 
                self.gt
            )), 
            columns = ["chr", 'pos', 'gt']
        )
        input_df.to_csv( outFile, sep = "\t", index = None, header = False  )


def import_vcf_file( inFile, logDebug = False, samples_to_load = [0], add_fields = None):
    """
    Function to read a VCF file and load required data. Wrapper for allel.read_vcf
    add_field = array with required names that needed to be loaded
    """
    if logDebug:
        vcf = allel.read_vcf(inFile, samples = samples_to_load, fields = '*')
    else:
        import io
        import sys
        sys.stderr = io.StringIO()
        vcf = allel.read_vcf(inFile, samples = samples_to_load, fields = '*')
        sys.stderr = sys.__stderr__
    snp_inputs = {}
    snp_inputs['samples'] = np.array(vcf['samples']).astype('U')
    if 'calldata/GT' in sorted(vcf.keys()):
        snp_inputs['gt'] = allel.GenotypeArray(vcf['calldata/GT']).to_gt().astype('U')
    else:
        snpmatch.die("input VCF file doesnt have required GT field")
    if 'calldata/PL' in sorted(vcf.keys()):
        snp_inputs['wei'] = np.copy(vcf['calldata/PL']).astype('float')
    snp_inputs['chr'] = np.array(vcf['variants/CHROM'], dtype="str").astype('U')
    snp_inputs['pos'] = np.array(vcf['variants/POS'])
    # if 'calldata/DP' in sorted(vcf.keys()):
    #     snp_inputs['dp'] = vcf['calldata/DP']
    if 'variants/DP' in sorted(vcf.keys()):
        snp_inputs['dp'] = vcf['variants/DP']
    else:
        snp_inputs['dp'] = np.repeat("NA", snp_inputs['pos'].shape[0])
    if add_fields is not None:
        for ef in add_fields:
            if ef in sorted(vcf.keys()):
                snp_inputs[ef] = vcf[ef]
            else:
                log.warn( "Field %s is not present in the VCF file" % ef )
    return(snp_inputs)


def potatoParser(inFile, logDebug, outFile = "parser"):
    inputs = ParseInputs(inFile, logDebug, outFile)
    return(inputs.chrs, inputs.pos, inputs.gt, inputs.wei, inputs.dp)
