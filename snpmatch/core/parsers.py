import pandas as pd
import numpy as np
import allel
import snpmatch
import logging
import os
import json

log = logging.getLogger(__name__)

def parseGT(snpGT):
    first = snpGT[0]
    snpBinary = np.zeros(len(snpGT), dtype = "int8")
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
    hetGT = "0" + separator + "1"
    refGT = "0" + separator + "0"
    altGT = "1" + separator + "1"
    nocall = "." + separator + "."
    snpBinary[np.where(snpGT == altGT)[0]] = 1
    snpBinary[np.where(snpGT == hetGT)[0]] = 2
    snpBinary[np.where(snpGT == nocall)[0]] = -1
    return snpBinary

def parseChrName(targetCHR):
    snpCHROM = np.char.replace(np.core.defchararray.lower(np.array(targetCHR, dtype="string")), "chr", "")
    snpsREQ = np.where(np.char.isdigit(snpCHROM))[0]   ## Filtering positions from mitochondrial and chloroplast
    snpCHR = snpCHROM[snpsREQ]
    return (snpCHR, snpsREQ)

def readBED(inFile, logDebug):
    log.info("reading the position file")
    targetSNPs = pd.read_table(inFile, header=None, usecols=[0,1,2])
    (snpCHR, snpsREQ) = parseChrName(targetSNPs[0])
    snpPOS = np.array(targetSNPs[1], dtype=int)[snpsREQ]
    snpGT = np.array(targetSNPs[2])[snpsREQ]
    snpBinary = parseGT(snpGT)
    snpWEI = np.ones((len(snpCHR), 3))  ## for homo and het
    snpWEI[np.where(snpBinary != 0),0] = 0
    snpWEI[np.where(snpBinary != 1),2] = 0
    snpWEI[np.where(snpBinary != 2),1] = 0
    return (snpCHR, snpPOS, snpGT, snpWEI)

def readVcf(inFile, logDebug):
    log.info("reading the VCF file")
    ## We read only one sample from the VCF file
    if logDebug:
        vcf = allel.read_vcf(inFile, samples = [0], fields = '*')
    else:
        import StringIO
        import sys
        sys.stderr = StringIO.StringIO()
        vcf = allel.read_vcf(inFile, samples = [0], fields = '*')
        #vcf = vcfnp.variants(inFile, cache=False).view(np.recarray)
        #vcfD = vcfnp.calldata_2d(inFile, cache=False).view(np.recarray)
        sys.stderr = sys.__stderr__
    (snpCHR, snpsREQ) = parseChrName(vcf['variants/CHROM'])
    try:
        snpGT = allel.GenotypeArray(vcf['calldata/GT']).to_gt()[snpsREQ, 0]
    except AttributeError:
        snpmatch.die("input VCF file doesnt have required GT field")
    snpsREQ = snpsREQ[np.where(snpGT != './.')[0]]
    snpGT = allel.GenotypeArray(vcf['calldata/GT']).to_gt()[snpsREQ, 0]
    if 'calldata/PL' in sorted(vcf.keys()):
        snpWEI = np.copy(vcf['calldata/PL'][snpsREQ, 0]).astype('float')
        snpWEI = snpWEI/(-10)
        snpWEI = np.exp(snpWEI)

    else:
        snpBinary = parseGT(snpGT)
        snpWEI = np.ones((len(snpsREQ), 3))  ## for homo and het
        snpWEI[np.where(snpBinary != 0),0] = 0
        snpWEI[np.where(snpBinary != 1),2] = 0
        snpWEI[np.where(snpBinary != 2),1] = 0
    snpCHR = snpCHR[snpsREQ]
    DPmean = np.mean(vcf['calldata/DP'][snpsREQ,0])
    snpPOS = np.array(vcf['variants/POS'][snpsREQ])
    return (DPmean, snpCHR, snpPOS, snpGT, snpWEI)

def parseInput(inFile, logDebug, outFile = "parser"):
    if outFile == "parser" or not outFile:
        outFile = inFile + ".snpmatch"
    if os.path.isfile(inFile + ".snpmatch.npz"):
        log.info("snpmatch parser dump found! loading %s", inFile + ".snpmatch.npz")
        snps = np.load(inFile + ".snpmatch.npz")
        (snpCHR, snpPOS, snpGT, snpWEI, DPmean) = (snps['chr'], snps['pos'], snps['gt'], snps['wei'], snps['dp'])
    else:
        _,inType = os.path.splitext(inFile)
        if inType == '.npz':
            log.info("loading snpmatch parser file! %s", inFile)
            snps = np.load(inFile)
            (snpCHR, snpPOS, snpGT, snpWEI, DPmean) = (snps['chr'], snps['pos'], snps['gt'], snps['wei'], snps['dp'])
        else:
            log.info('running snpmatch parser!')
            if inType == '.vcf':
                (DPmean, snpCHR, snpPOS, snpGT, snpWEI) = readVcf(inFile, logDebug)
            elif inType == '.bed':
                (snpCHR, snpPOS, snpGT, snpWEI) = readBED(inFile, logDebug)
                DPmean = "NA"
            else:
                snpmatch.die("input file type %s not supported" % inType)
            log.info("creating snpmatch parser file: %s", outFile + '.npz')
            np.savez(outFile, chr = snpCHR, pos = snpPOS, gt = snpGT, wei = snpWEI, dp = DPmean)
            NumSNPs = len(snpCHR)
            case = 0
            note = "Sufficient number of SNPs"
            if NumSNPs < snpmatch.snp_thres:
                note = "Attention: low number of SNPs provided"
                case = 1
            snpst = np.unique(snpCHR, return_counts=True)
            snpdict = dict(('Chr%s' % snpst[0][i], snpst[1][i]) for i in range(len(snpst[0])))
            statdict = {"interpretation": {"case": case, "text": note}, "snps": snpdict, "num_of_snps": NumSNPs, "depth": DPmean}
            statdict['percent_heterozygosity'] = snpmatch.getHeterozygosity(snpGT)
            with open(outFile + ".stats.json", "w") as out_stats:
                out_stats.write(json.dumps(statdict))
    log.info("done!")
    return (snpCHR, snpPOS, snpGT, snpWEI, DPmean)
