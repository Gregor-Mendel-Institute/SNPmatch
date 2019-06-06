"""
  SNPmatch for crosses
"""
import numpy as np
import numpy.ma
import pandas as pd
import logging
from . import snpmatch
from . import genomes
from . import snp_genotype
from . import parsers
import json
import itertools

log = logging.getLogger(__name__)
chunk_size = 1000

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
    inputs.filter_chr_names()
    num_lines = len(GenotypeData.accessions)
    NumMatSNPs = 0
    TotScoreList = np.zeros(num_lines, dtype="uint32")
    TotNumInfoSites = np.zeros(num_lines, dtype="uint32")
    TotMatchedTarInds = np.zeros(0, dtype="int")
    iter_bins_genome = genome.get_bins_genome(GenotypeData, binLen)
    iter_bins_snps = genome.get_bins_arrays(inputs.chrs, inputs.pos, binLen)
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
        TotMatchedTarInds = np.append(TotMatchedTarInds, matchedTarInd)
        writeBinData(out_file, bin_inds, GenotypeData, ScoreList, NumInfoSites)
        winds_chrs = np.append( winds_chrs, genome.chrs_ids[e_g[0]] )
        if bin_inds % 50 == 0:
            log.info("Done analysing %s positions", NumMatSNPs)
        bin_inds += 1
    out_file.close()
    overlap = float(NumMatSNPs)/len(inputs.pos)
    result = snpmatch.GenotyperOutput(GenotypeData.accessions, TotScoreList, TotNumInfoSites, overlap, NumMatSNPs, inputs.dp)
    result.matchedTarInd = TotMatchedTarInds
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
    likeLiwind = pd.read_csv(outFile, sep = "\t", header=None)
    ScoreAcc = pd.read_csv(scoreFile, sep = "\t", header=None)
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
            out_stats.write(json.dumps(topHitsDict, sort_keys=True, indent=4))

def crossIdentifier(inputs, g, binLen, outID):
    ## Get tophit accessions
    # sorting based on the final scores
    inputs.filter_chr_names()
    if not outID:
        outID = "cross.identifier"
    outFile = outID + '.windowscore.txt'
    scoreFile = outID + '.scores.txt'
    snpmatch_result = crossWindower(inputs, g.g, binLen, outFile)
    snpmatch_result.print_json_output( scoreFile + ".matches.json" )
    snpmatch.getHeterozygosity(inputs.gt[snpmatch_result.matchedTarInd],  scoreFile + ".matches.json")
    log.info("simulating F1s for top 10 accessions")
    TopHitAccs = np.argsort(-snpmatch_result.probabilies)[0:10]
    commonSNPs = g.get_positions_idxs( inputs.chrs, inputs.pos )
    for (i, j) in itertools.combinations(TopHitAccs, 2):
        p1 = g.g_acc.snps[:,i]
        p2 = g.g_acc.snps[:,j]
        score = 0
        numinfo = 0
        for ind in range(0, len(commonSNPs[0]), chunk_size):
            matchedAccInd = commonSNPs[0][ind:ind+chunk_size]
            matchedTarInd = commonSNPs[1][ind:ind+chunk_size]
            gtp1 = p1[matchedAccInd]
            gtp2 = p2[matchedAccInd]
            matchedTarWEI = inputs.wei[matchedTarInd,]
            homalt = np.where((gtp1 == 1) & (gtp2 == 1))[0]
            homref = np.where((gtp1 == 0) & (gtp2 == 0))[0]
            het = np.where((gtp1 != -1) & (gtp2 != -1) & (gtp1 != gtp2))[0]
            score = score + np.sum(matchedTarWEI[homalt, 2]) + np.sum(matchedTarWEI[homref, 0]) + np.sum(matchedTarWEI[het, 1])
            numinfo = numinfo + len(homalt) + len(homref) + len(het)
        snpmatch_result.scores = np.append(snpmatch_result.scores, score)
        snpmatch_result.ninfo = np.append(snpmatch_result.ninfo, numinfo)
        snpmatch_result.accs = np.append( snpmatch_result.accs, g.g.accessions[i] + "x" + g.g.accessions[j] )
    log.info("writing output!")
    snpmatch_result.print_out_table( scoreFile )
    crossInterpreter(snpmatch_result, g.g, binLen, outID)

def potatoCrossIdentifier(args):
    inputs = parsers.ParseInputs(inFile = args['inFile'], logDebug = args['logDebug'])
    global genome
    genome = genomes.Genome(args['genome'])
    log.info("loading genotype files!")
    g = snp_genotype.Genotype(args['hdf5File'], args['hdf5accFile'])
    log.info("done!")
    log.info("running cross identifier!")
    crossIdentifier(inputs, g, args['binLen'], args['outFile'])
    log.info("finished!")
