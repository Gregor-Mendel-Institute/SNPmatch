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
np_test_identity = np.vectorize(snpmatch.test_identity, excluded=["pthres", "error_rate"])
np_get_fraction = np.vectorize(snpmatch.get_fraction)

class CrossIdentifier(object):
    ## class object for main CSMATCH

    def __init__(self, inputs, g, genome_id, binLen, output_id = "cross.identifier", run_identifier = True, identity_error_rate = 0.02, skip_db_hets = False):
        self.g = g
        assert type(inputs) is parsers.ParseInputs, "provide a parsers class"
        inputs.filter_chr_names()
        self.inputs = inputs
        self.genome = genomes.Genome(genome_id)
        self.binLen = binLen
        self.output_id = output_id
        self.error_rate = identity_error_rate
        self._skip_db_hets = skip_db_hets
        if run_identifier:
            self.cross_identifier()

    def cross_identifier(self):
        window_snpmatch_result = self.window_genotyper(self.output_id + '.windowscore.txt')
        window_snpmatch_result.print_json_output( self.output_id + ".scores.txt.matches.json" )
        snpmatch.getHeterozygosity( self.inputs.gt[window_snpmatch_result.matchedTarInd],  self.output_id + ".scores.txt.matches.json" )
        with open(self.output_id + ".scores.txt.matches.json") as json_out:
            self.cross_identfier_json = json.load(json_out)
        self.result = self.match_insilico_f1s(window_snpmatch_result, self.output_id + '.scores.txt')
        self.cross_interpreter( self.output_id + ".matches.json" )

    @staticmethod
    def get_window_data(bin_inds, AccList, ScoreList, NumInfoSites, error_rate=0.02):
        num_lines = len(AccList)
        (likeliScore, likeliHoodRatio) = snpmatch.GenotyperOutput.calculate_likelihoods(ScoreList, NumInfoSites)
        identity = np_test_identity(x = ScoreList, n = NumInfoSites, error_rate = error_rate)
        NumAmb = np.where(likeliHoodRatio < snpmatch.lr_thres)[0]
        t_window = pd.DataFrame( np.column_stack((AccList, ScoreList, NumInfoSites, np_get_fraction(ScoreList, NumInfoSites), likeliScore, identity )), columns = ["acc", "snps_match", "snps_info", "score", "likelihood", "identical"] )
        t_window["num_amb"] = len(NumAmb)
        t_window['window_index'] = bin_inds
        t_window["acc"] = t_window["acc"].apply(str)
        t_window["snps_match"] = t_window["snps_match"].apply(float).apply(int)
        t_window["snps_info"] = t_window["snps_info"].apply(float).apply(int)
        t_window["identical"] = t_window["identical"].apply(float).apply(int)
        if len(NumAmb) >= 1 and len(NumAmb) < num_lines:
            t_window = t_window.iloc[NumAmb,:]
        else:
            t_window = pd.DataFrame( columns = ["acc", "snps_match", "snps_info", "score", "likelihood", "identical", "num_amb", "window_index"] )
        return(t_window)
        # out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (AccList[k], int(ScoreList[k]), NumInfoSites[k], score, likeliScore[k], nextLikeli, len(NumAmb), bin_inds, identity[k]))

    def window_genotyper(self, out_file, mask_acc_ix = None):
        num_lines = len(self.g.accessions)
        NumMatSNPs = 0
        if mask_acc_ix is not None:
            assert type(mask_acc_ix) is np.ndarray, "please provide numpy array of acc indices to be masked"
            mask_acc_to_print = np.setdiff1d(np.arange( num_lines ), mask_acc_ix)
        else:
            mask_acc_to_print = np.arange( num_lines )
        TotScoreList = np.zeros(num_lines, dtype="uint32")
        TotNumInfoSites = np.zeros(num_lines, dtype="uint32")
        TotMatchedTarInds = np.zeros(0, dtype="int")
        iter_bins_genome = self.genome.get_bins_genome(self.g.g, self.binLen)
        iter_bins_snps = self.genome.get_bins_arrays(self.inputs.chrs, self.inputs.pos, self.binLen)
        self.windows_data = pd.DataFrame( columns = ["acc", "snps_match", "snps_info", "score", "likelihood", "identical", "num_amb", "window_index"] )
        bin_inds = 1
        winds_chrs = np.zeros(0, dtype = self.g.g.chrs.dtype)
        for e_g, e_s in zip(iter_bins_genome, iter_bins_snps):
            g_bin_pos = self.g.g.positions[e_g[2]]
            perchrtarSNPpos = self.inputs.pos[e_s[2]]
            matchedAccInd = np.array(e_g[2], dtype=int)[np.where(np.in1d(g_bin_pos, perchrtarSNPpos))[0]]
            matchedTarInd = np.array(e_s[2], dtype=int)[np.where(np.in1d(perchrtarSNPpos, g_bin_pos))[0]]
            NumMatSNPs = NumMatSNPs + len(matchedAccInd)
            if len(matchedAccInd) > 0:
                ScoreList, NumInfoSites = snpmatch.matchGTsAccs( self.inputs.wei[matchedTarInd,], self.g.g.snps[matchedAccInd,:], self._skip_db_hets )
                TotScoreList = TotScoreList + ScoreList
                TotNumInfoSites = TotNumInfoSites + NumInfoSites
                TotMatchedTarInds = np.append(TotMatchedTarInds, matchedTarInd)
                self.windows_data = self.windows_data.append( self.get_window_data(bin_inds, self.g.accessions[mask_acc_to_print], ScoreList[mask_acc_to_print], NumInfoSites[mask_acc_to_print], self.error_rate), ignore_index=True )
            if bin_inds % 50 == 0:
                log.info("Done analysing %s positions", NumMatSNPs)
            winds_chrs = np.append( winds_chrs, self.genome.chrs_ids[e_g[0]] )
            bin_inds += 1
        overlap = snpmatch.get_fraction(NumMatSNPs, len(self.inputs.pos))
        result = snpmatch.GenotyperOutput(self.g.accessions[mask_acc_to_print], TotScoreList[mask_acc_to_print], TotNumInfoSites[mask_acc_to_print], overlap, NumMatSNPs, self.inputs.dp)
        result.matchedTarInd = TotMatchedTarInds
        result.winds_chrs = winds_chrs
        if out_file is not None:
            self.windows_data.to_csv(out_file, sep = "\t", index = False)
            return(result)
        else:
            return([self.windows_data, result])

    def match_insilico_f1s(self, snpmatch_result, out_file):
        ## Get tophit accessions
        # sorting based on the final scores
        assert type(snpmatch_result) is snpmatch.GenotyperOutput, "Please provide GenotyperOutput class as input"
        if not hasattr(snpmatch_result, 'probabilies'):
            snpmatch_result.get_probabilities()
        log.info("simulating F1s for top 10 accessions")
        TopHitAccs = np.argsort(-snpmatch_result.probabilies)[0:10]
        commonSNPs = self.g.get_positions_idxs( self.inputs.chrs, self.inputs.pos )
        for (i, j) in itertools.combinations(TopHitAccs, 2):
            p1 = self.g.g_acc.snps[:,i]
            p2 = self.g.g_acc.snps[:,j]
            score = 0
            numinfo = 0
            for ind in range(0, len(commonSNPs[0]), chunk_size):
                matchedAccInd = commonSNPs[0][ind:ind+chunk_size]
                matchedTarInd = commonSNPs[1][ind:ind+chunk_size]
                gtp1 = p1[matchedAccInd]
                gtp2 = p2[matchedAccInd]
                matchedTarWEI = self.inputs.wei[matchedTarInd,]
                homalt = np.where((gtp1 == 1) & (gtp2 == 1))[0]
                homref = np.where((gtp1 == 0) & (gtp2 == 0))[0]
                het = np.where((gtp1 != -1) & (gtp2 != -1) & (gtp1 != gtp2))[0]
                score = score + np.sum(matchedTarWEI[homalt, 2]) + np.sum(matchedTarWEI[homref, 0]) + np.sum(matchedTarWEI[het, 1])
                numinfo = numinfo + len(homalt) + len(homref) + len(het)
            snpmatch_result.scores = np.append(snpmatch_result.scores, score)
            snpmatch_result.ninfo = np.append(snpmatch_result.ninfo, numinfo)
            snpmatch_result.accs = np.append( snpmatch_result.accs, self.g.accessions[i] + "x" + self.g.accessions[j] )
        if out_file is not None:
            snpmatch_result.print_out_table( out_file )
        return(snpmatch_result)

    def cross_interpreter(self, out_file):
        assert 'cross_identfier_json' in dir(self), "run cross identifier first!"
        assert 'windows_data' in dir(self), "run window genotyper first!"
        cs_thres = 0.9
        log.info("running cross interpreter!")
        if self.cross_identfier_json['interpretation']['case'] >= 3:
            identical_wind = np.where(self.windows_data.groupby('window_index').max()['identical'] == 1)[0]
            num_winds = np.unique(self.windows_data['window_index']).shape[0]
            self.cross_identfier_json['identical_windows'] = [ snpmatch.get_fraction( identical_wind.shape[0], num_winds ), num_winds ]
            homo_wind = np.intersect1d(self.windows_data['window_index'][np.where(self.windows_data['num_amb'] < 20 )[0]], identical_wind)
            homo_acc = np.unique(self.windows_data.iloc[:,0][np.where(np.in1d(self.windows_data.iloc[:,7], homo_wind))[0]],return_counts=True)
            matches_dict = [(homo_acc[0][i], homo_acc[1][i]) for i in np.argsort(-homo_acc[1])]
            self.cross_identfier_json['matches'] = matches_dict
            topMatch = np.argsort(self.result.likelis)[0]  ## top F1 match sorted based on likelihood
            if topMatch in np.where(~np.in1d(self.result.accs, self.g.accessions))[0]:
                mother = self.result.accs[topMatch].split("x")[0]
                father = self.result.accs[topMatch].split("x")[1]
                self.cross_identfier_json['interpretation']['text'] = "Sample may be a F1! or a contamination!"
                self.cross_identfier_json['interpretation']['case'] = 5
                self.cross_identfier_json['parents'] = {'mother': [mother,1], 'father': [father,1]}
                self.cross_identfier_json['genotype_windows'] = {'chr_bins': None, 'coordinates': {'x': None, 'y': None}}
            else:
                ## Get exactly the homozygous windows with one count
                clean = np.unique(self.windows_data.iloc[:,0][np.where(self.windows_data.iloc[:,6] == 1)[0]], return_counts = True)
                if len(clean[0]) > 0:  ## Check if there are atlease one homozygous window
                    parents = clean[0][np.argsort(-clean[1])[0:2]].astype("string")
                    parents_counts = clean[1][np.argsort(-clean[1])[0:2]].astype("int")
                    xdict = np.array(np.unique(self.windows_data.iloc[:,7]), dtype="int")
                    ydict = np.repeat("NA", len(xdict)).astype("a25")
                    if len(parents) == 1:
                        self.cross_identfier_json['interpretation']['text'] = "Sample may be a F2! but only one parent found!"
                        self.cross_identfier_json['interpretation']['case'] = 6
                        self.cross_identfier_json['parents'] = {'mother': [parents[0], parents_counts[0]], 'father': ["NA", "NA"]}
                        par1_ind = self.windows_data.iloc[:,7][np.where((self.windows_data.iloc[:,0].astype("string") == parents[0]) & np.in1d(self.windows_data.iloc[:,7], homo_wind))[0]]
                        ydict[np.where(np.in1d(xdict,par1_ind))[0]] = parents[0]
                        chr_bins = None
                    else:
                        self.cross_identfier_json['interpretation']['text'] = "Sample may be a F2!"
                        self.cross_identfier_json['interpretation']['case'] = 6
                        self.cross_identfier_json['parents'] = {'mother': [parents[0], parents_counts[0]], 'father': [parents[1], parents_counts[1]]}
                        NumChrs = np.unique(self.result.winds_chrs, return_counts=True)
                        chr_bins = dict(( NumChrs[0][i], NumChrs[1][i]) for i in range(len(NumChrs[0])))
                        par1_ind = np.array(self.windows_data.iloc[:,7][np.where((self.windows_data.iloc[:,0].astype("string") == parents[0]) & np.in1d(self.windows_data.iloc[:,7], homo_wind))[0]])
                        par2_ind = np.array(self.windows_data.iloc[:,7][np.where((self.windows_data.iloc[:,0].astype("string") == parents[1]) & np.in1d(self.windows_data.iloc[:,7], homo_wind))[0]])
                        ydict[np.where(np.in1d(xdict,par1_ind))[0]] = parents[0]
                        ydict[np.where(np.in1d(xdict,par2_ind))[0]] = parents[1]
                    xdict = xdict.tolist()
                    ydict = ydict.tolist()
                    self.cross_identfier_json['genotype_windows'] = {'chr_bins': chr_bins, 'coordinates': {'x': xdict, 'y': ydict}}
                else:   ## No homozygous window found!
                    self.cross_identfier_json['interpretation']['case'] = 7
                    self.cross_identfier_json['interpretation']['text'] = "Sample may just be contamination!"
                    self.cross_identfier_json['genotype_windows'] = {'chr_bins': None, 'coordinates': {'x': None, 'y': None}}
                    self.cross_identfier_json['parents'] = {'mother': [None,0], 'father': [None,1]}
            with open(out_file, "w") as out_stats:
                out_stats.write(json.dumps(self.cross_identfier_json, sort_keys=True, indent=4))


def potatoCrossIdentifier(args):
    inputs = parsers.ParseInputs(inFile = args['inFile'], logDebug = args['logDebug'])
    log.info("loading genotype files!")
    g = snp_genotype.Genotype(args['hdf5File'], args['hdf5accFile'])
    log.info("done!")
    log.info("running cross identifier!")
    ci = CrossIdentifier(inputs, g, args['genome'], args['binLen'], args['outFile'], run_identifier = True, skip_db_hets = args['skip_db_hets'])
    log.info("finished!")
