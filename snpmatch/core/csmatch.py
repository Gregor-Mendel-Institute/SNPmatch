"""
  SNPmatch for crosses
"""
import numpy as np
import numpy.ma
import scipy.stats as st
from pygwas.core import genotype
import pandas as pd
import logging
import os
import snpmatch
import parsers
import json
import itertools

log = logging.getLogger(__name__)

# Arabidopsis chromosome length
chrlen = np.array((30427671, 19698289, 23459830, 18585056, 26975502, 154478, 154478))
tair_chrs = ['1', '2', '3', '4', '5']
mean_recomb_rates = [3.4, 3.6, 3.5, 3.8, 3.6]  ## cM/Mb ## Salome, P et al. 2011

def get_bins_echr(real_chrlen, chr_pos, binLen, rel_ix):
    ind = 0
    for t in range(1, real_chrlen, binLen):
        skipped = True
        result = []
        bin_bed = [int(t), int(t) + binLen - 1]
        for epos in chr_pos[ind:]:
            if epos >= bin_bed[0]:
                if epos <= bin_bed[1]:
                    result.append(ind + rel_ix)
                elif epos > bin_bed[1]:
                    skipped = False
                    yield((bin_bed, result))
                    break
                ind = ind + 1
        if skipped:
            yield((bin_bed, result))

def get_bins_genome(g, binLen):
    binLen = int(binLen)
    g_chrs = np.char.replace(np.core.defchararray.lower(np.array(g.chrs[0:len(tair_chrs)], dtype="string")), "chr", "")
    if len(g_chrs) > 7:
        snpmatch.die("Please change the genome sizes in csmatch module")
    if not np.array_equal(g_chrs, np.array(tair_chrs)):
        snpmatch.die("Please change the genome sizes in csmatch module")
    for chr_ix in range(len(g_chrs)):
        start = g.chr_regions[chr_ix][0]
        end = g.chr_regions[chr_ix][1]
        chr_pos = g.positions[start:end]
        echr_bins = get_bins_echr(chrlen[chr_ix], chr_pos, binLen, start)
        for e_bin in echr_bins:
            yield((chr_ix, e_bin[0], e_bin[1]))

def get_bins_arrays(g_chrs, g_snppos, binLen):
    g_chrs = np.char.replace(np.core.defchararray.lower(np.array(g_chrs, dtype="string")), "chr", "")
    g_chrs_uq = np.unique(g_chrs)
    matched_tair = np.where( np.in1d(g_chrs_uq, tair_chrs) )[0]
    if len(matched_tair) == 0:
        snpmatch.die("Please change the genome sizes in csmatch module")
    for chr_ix in range(len(tair_chrs)):
        chr_pos_ix = np.where(g_chrs == tair_chrs[chr_ix])[0]
        if len(chr_pos_ix) > 0:
            echr_bins = get_bins_echr(chrlen[chr_ix], g_snppos[chr_pos_ix], binLen, chr_pos_ix[0])
        else:
            echr_bins = get_bins_echr(chrlen[chr_ix], g_snppos[chr_pos_ix], binLen, 0)
        for e_bin in echr_bins:
            yield((chr_ix, e_bin[0], e_bin[1]))


def writeBinData(out_file, bin_inds, GenotypeData, ScoreList, NumInfoSites):
    num_lines = len(GenotypeData.accessions)
    (likeliScore, likeliHoodRatio) = snpmatch.GenotyperOutput.calculate_likelihoods(ScoreList, NumInfoSites)
    if len(likeliScore) > 0:
        NumAmb = np.where(likeliHoodRatio < snpmatch.lr_thres)[0]
        if len(NumAmb) >= 1 and len(NumAmb) < num_lines:
            try:
                nextLikeli = np.nanmin(likeliHoodRatio[np.where(likeliHoodRatio > snpmatch.lr_thres)[0]])
            except:
                nextLikeli = 1
            for k in NumAmb:
                score = float(ScoreList[k])/NumInfoSites[k]
                out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[k], int(ScoreList[k]), NumInfoSites[k], score, likeliScore[k], nextLikeli, len(NumAmb), bin_inds))

def crossWindower(inputs, GenotypeData, binLen, outFile):
    inputs.filter_chr_names(GenotypeData)
    num_lines = len(GenotypeData.accessions)
    NumMatSNPs = 0
    chunk_size = 1000
    TotScoreList = np.zeros(num_lines, dtype="uint32")
    TotNumInfoSites = np.zeros(num_lines, dtype="uint32")
    iter_bins_genome = get_bins_genome(GenotypeData, binLen)
    iter_bins_snps = get_bins_arrays(inputs.chrs, inputs.pos, binLen)
    out_file = open(outFile, 'w')
    bin_inds = 1
    winds_chrs = np.zeros(0, dtype = GenotypeData.chrs.dtype)
    for e_g, e_s in itertools.izip(iter_bins_genome, iter_bins_snps):
        g_bin_pos = GenotypeData.positions[e_g[2]]
        perchrtarSNPpos = inputs.pos[e_s[2]]
        matchedAccInd = np.array(e_g[2], dtype=int)[np.where(np.in1d(g_bin_pos, perchrtarSNPpos))[0]]
        matchedTarInd = np.array(e_s[2], dtype=int)[np.where(np.in1d(perchrtarSNPpos, g_bin_pos))[0]]
        matchedTarWei = inputs.wei[matchedTarInd,]
        TarGTs0 = np.zeros(len(matchedTarInd), dtype="int8")
        TarGTs1 = np.ones(len(matchedTarInd), dtype="int8") + 1
        TarGTs2 = np.ones(len(matchedTarInd), dtype="int8")
        NumMatSNPs = NumMatSNPs + len(matchedAccInd)
        ScoreList = np.zeros(num_lines, dtype="uint32")
        NumInfoSites = np.zeros(num_lines, dtype="uint32")
        for j in range(0, len(matchedAccInd), chunk_size):
            t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
            samSNPs0 = np.reshape(np.repeat(TarGTs0[j:j+chunk_size], num_lines), (len(TarGTs0[j:j+chunk_size]),num_lines))
            samSNPs1 = np.reshape(np.repeat(TarGTs1[j:j+chunk_size], num_lines), (len(TarGTs1[j:j+chunk_size]),num_lines))
            samSNPs2 = np.reshape(np.repeat(TarGTs2[j:j+chunk_size], num_lines), (len(TarGTs2[j:j+chunk_size]),num_lines))
            tempScore0 = np.sum(np.multiply(np.array(t1001SNPs == samSNPs0, dtype=int).T, matchedTarWei[j:j+chunk_size,0]).T, axis=0)
            tempScore1 = np.sum(np.multiply(np.array(t1001SNPs == samSNPs1, dtype=int).T, matchedTarWei[j:j+chunk_size,1]).T, axis=0)
            tempScore2 = np.sum(np.multiply(np.array(t1001SNPs == samSNPs2, dtype=int).T, matchedTarWei[j:j+chunk_size,2]).T, axis=0)
            ScoreList = ScoreList + tempScore0 + tempScore1 + tempScore2
            if(len(TarGTs0[j:j+chunk_size]) >= 1):
                NumInfoSites = NumInfoSites + len(TarGTs0[j:j+chunk_size]) - np.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0)
        TotScoreList = TotScoreList + ScoreList
        TotNumInfoSites = TotNumInfoSites + NumInfoSites
        writeBinData(out_file, bin_inds, GenotypeData, ScoreList, NumInfoSites)
        winds_chrs = np.append( winds_chrs, inputs.chr_list[e_g[0]] )
        if bin_inds % 50 == 0:
            log.info("Done analysing %s positions", NumMatSNPs)
        bin_inds += 1
    out_file.close()
    overlap = float(NumMatSNPs)/len(inputs.filter_inds_ix)
    result = snpmatch.GenotyperOutput(GenotypeData.accessions, TotScoreList, TotNumInfoSites, overlap, NumMatSNPs, inputs.dp)
    result.winds_chrs = winds_chrs
    return(result)

def getHomoWindows(likeLiwind):
    snp_thres_wind = np.nanmean(likeLiwind[2]) - np.std(likeLiwind[2])
    x_info = np.unique(likeLiwind[7])
    homo_wind = np.zeros(0, dtype = "int")
    for i in x_info:
        eWinds = likeLiwind.iloc[np.where(likeLiwind[7] == i)[0]]
        if np.nanmean(eWinds[3]) > snpmatch.prob_thres and np.nanmean(eWinds[2]) > snp_thres_wind:
            homo_wind = np.append(homo_wind, i)
    return(homo_wind)

def crossInterpreter(snpmatch_result, GenotypeData, binLen, outID):
    ## ScoreFile should be one from crossF1genotyper
    ## Output file is from the crossIdentifier
    cs_thres = 0.9
    outFile = outID + '.windowscore.txt'
    scoreFile = outID + '.scores.txt'
    log.info("running cross interpreter!")
    likeLiwind = pd.read_table(outFile, header=None)
    ScoreAcc = pd.read_table(scoreFile, header=None)
    topHitsDict = json.load(open(scoreFile + ".matches.json", 'r'))
    if topHitsDict['interpretation']['case'] == 3:
        homo_wind = getHomoWindows(likeLiwind)
        homo_acc = np.unique(likeLiwind[0][np.where(np.in1d(likeLiwind[7], homo_wind))[0]],return_counts=True)
        matches_dict = [(homo_acc[0][i].astype("string"), homo_acc[1][i]) for i in np.argsort(-homo_acc[1])]
        topHitsDict['matches'] = matches_dict
        f1matches = ScoreAcc.iloc[~np.in1d(ScoreAcc[0], GenotypeData.accessions)].reset_index()
        topMatch = np.argsort(f1matches[5])[0]  ## top F1 match sorted based on likelihood
        if f1matches[3][topMatch] > cs_thres:
            mother = f1matches[0][topMatch].split("x")[0]
            father = f1matches[0][topMatch].split("x")[1]
            topHitsDict['interpretation']['text'] = "Sample may be a F1! or a contamination!"
            topHitsDict['interpretation']['case'] = 5
            topHitsDict['parents'] = {'mother': [mother,1], 'father': [father,1]}
            topHitsDict['genotype_windows'] = {'chr_bins': None, 'coordinates': {'x': None, 'y': None}}
        else:
            ## Get exactly the homozygous windows with one count
            clean = np.unique(likeLiwind[0][np.where(likeLiwind[6] == 1)[0]], return_counts = True)
            if len(clean[0]) > 0:  ## Check if there are atlease one homozygous window
                parents = clean[0][np.argsort(-clean[1])[0:2]].astype("string")
                parents_counts = clean[1][np.argsort(-clean[1])[0:2]].astype("int")
                xdict = np.array(np.unique(likeLiwind[7]), dtype="int")
                ydict = np.repeat("NA", len(xdict)).astype("a25")
                if len(parents) == 1:
                    topHitsDict['interpretation']['text'] = "Sample may be a F2! but only one parent found!"
                    topHitsDict['interpretation']['case'] = 6
                    topHitsDict['parents'] = {'mother': [parents[0], parents_counts[0]], 'father': ["NA", "NA"]}
                    par1_ind = likeLiwind[7][np.where((likeLiwind[0].astype("string") == parents[0]) & np.in1d(likeLiwind[7], homo_wind))[0]]
                    ydict[np.where(np.in1d(xdict,par1_ind))[0]] = parents[0]
                else:
                    topHitsDict['interpretation']['text'] = "Sample may be a F2!"
                    topHitsDict['interpretation']['case'] = 6
                    topHitsDict['parents'] = {'mother': [parents[0], parents_counts[0]], 'father': [parents[1], parents_counts[1]]}
                    NumChrs = np.unique(snpmatch_result.winds_chrs, return_counts=True)
                    chr_bins = dict(( NumChrs[0][i], NumChrs[1][i]) for i in range(len(NumChrs[0])))
                    par1_ind = np.array(likeLiwind[7][np.where((likeLiwind[0].astype("string") == parents[0]) & np.in1d(likeLiwind[7], homo_wind))[0]])
                    par2_ind = np.array(likeLiwind[7][np.where((likeLiwind[0].astype("string") == parents[1]) & np.in1d(likeLiwind[7], homo_wind))[0]])
                    ydict[np.where(np.in1d(xdict,par1_ind))[0]] = parents[0]
                    ydict[np.where(np.in1d(xdict,par2_ind))[0]] = parents[1]
                xdict = xdict.tolist()
                ydict = ydict.tolist()
                topHitsDict['genotype_windows'] = {'chr_bins': chr_bins, 'coordinates': {'x': xdict, 'y': ydict}}
            else:   ## No homozygous window found!
                topHitsDict['interpretation']['case'] = 7
                topHitsDict['interpretation']['text'] = "Sample may just be contamination!"
                topHitsDict['genotype_windows'] = {'chr_bins': None, 'coordinates': {'x': None, 'y': None}}
                topHitsDict['parents'] = {'mother': [None,0], 'father': [None,1]}
        with open(outID + ".matches.json", "w") as out_stats:
            out_stats.write(json.dumps(topHitsDict))

def crossIdentifier(inputs, GenotypeData, GenotypeData_acc, binLen, outID):
    ## Get tophit accessions
    # sorting based on the final scores
    inputs.filter_chr_names(GenotypeData)
    if not outID:
        outID = "cross.identifier"
    outFile = outID + '.windowscore.txt'
    scoreFile = outID + '.scores.txt'
    snpmatch_result = crossWindower(inputs, GenotypeData, binLen, outFile)
    snpmatch_result.print_json_output( scoreFile + ".matches.json" )
    log.info("simulating F1s for top 10 accessions")
    TopHitAccs = np.argsort(-snpmatch_result.probabilies)[0:10]
    for (i, j) in itertools.combinations(TopHitAccs, 2):
        p1 = GenotypeData_acc.snps[:,i]
        p2 = GenotypeData_acc.snps[:,j]
        score = 0
        numinfo = 0
        NumMatSNPs = 0
        for ind,echr in enumerate(inputs.chr_list):
            perchrTarPos = np.where(inputs.chrs_nochr == echr)[0]
            perchrtarSNPpos = inputs.pos[perchrTarPos]
            start = GenotypeData.chr_regions[ind][0]
            end = GenotypeData.chr_regions[ind][1]
            chrpositions = GenotypeData.positions[start:end]
            matchedAccInd = np.where(np.in1d(chrpositions, perchrtarSNPpos))[0] + start
            matchedTarInd = np.where(np.in1d(perchrtarSNPpos, chrpositions))[0]
            NumMatSNPs = NumMatSNPs + len(matchedTarInd)
            gtp1 = p1[matchedAccInd]
            gtp2 = p2[matchedAccInd]
            matchedTarWEI = inputs.wei[perchrTarPos[matchedTarInd],]
            homalt = np.where((gtp1 == 1) & (gtp2 == 1))[0]
            homref = np.where((gtp1 == 0) & (gtp2 == 0))[0]
            het = np.where((gtp1 != -1) & (gtp2 != -1) & (gtp1 != gtp2))[0]
            score = score + np.sum(matchedTarWEI[homalt, 2]) + np.sum(matchedTarWEI[homref, 0]) + np.sum(matchedTarWEI[het, 1])
            numinfo = numinfo + len(homalt) + len(homref) + len(het)
        snpmatch_result.scores = np.append(snpmatch_result.scores, score)
        snpmatch_result.ninfo = np.append(snpmatch_result.ninfo, numinfo)
        snpmatch_result.accs = np.append( snpmatch_result.accs, GenotypeData.accessions[i] + "x" + GenotypeData.accessions[j] )
    log.info("writing output!")
    snpmatch_result.print_out_table( scoreFile )
    crossInterpreter(snpmatch_result, GenotypeData, binLen, outID)

def potatoCrossIdentifier(args):
    inputs = parsers.ParseInputs(inFile = args['inFile'], logDebug = args['logDebug'])
    log.info("loading genotype files!")
    GenotypeData = genotype.load_hdf5_genotype_data(args['hdf5File'])
    GenotypeData_acc = genotype.load_hdf5_genotype_data(args['hdf5accFile'])
    log.info("done!")
    log.info("running cross identifier!")
    crossIdentifier(inputs, GenotypeData, GenotypeData_acc, args['binLen'], args['outFile'])
    log.info("finished!")
