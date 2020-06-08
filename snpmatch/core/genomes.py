# Loading genome data
import logging
import numpy as np
import pandas as pd
import json
import os.path


log = logging.getLogger(__name__)

def getInd_bin_bed(bin_bed, tair10):
    ## bin_bed = ["Chr1", 1, 1000]
    bin_s = [int(bin_bed[0].replace("Chr", "")) - 1, int(bin_bed[1]), int(bin_bed[2])]
    return(tair10.chr_inds[bin_s[0]] + int(( bin_s[1] + bin_s[2] )/2) )

class Genome(object):
    ## coordinates for ArabidopsisGenome using TAIR 10

    def __init__(self, ref_json):
        ref_genomes_ids = self.get_genome_ids()
        if ref_json in ref_genomes_ids:
            ref_json = os.path.dirname(__file__) + '/../resources/genomes/' + ref_json + '.json'
        assert os.path.exists(ref_json), "Reference json file missing: %s" % ref_json
        with open(ref_json) as ref_genome:
            self.json = json.load(ref_genome)
        self.chrs = np.array(self.json['ref_chrs'], dtype="str")
        self.chrlen = np.array(self.json['ref_chrlen'], dtype = int)
        self.chrs_ids = np.char.replace(np.core.defchararray.lower(self.chrs), "chr", "")

    def get_genome_ids(self):
        from glob import glob
        ref_files = glob(os.path.dirname(__file__) + '/../resources/genomes/*.json')
        ref_ids = []
        for ef in ref_files:
            ref_ids.append( os.path.basename(ef).replace(".json", "") )
        return(ref_ids)

    def get_chr_ind(self, echr):
        real_chrs = np.array( [ ec.replace("Chr", "").replace("chr", "") for ec in self.chrs ] )
        if type(echr) is str or type(echr) is np.string_:
            echr_num = str(echr).replace("Chr", "").replace("chr", "")
            if len(np.where(real_chrs == echr_num)[0]) == 1:
                return(np.where(real_chrs == echr_num)[0][0])
            else:
                return(None)
        echr_num = np.unique( np.array( echr ) )
        ret_echr_ix = np.zeros( len(echr), dtype="int8" )
        for ec in echr_num:
            t_ix = np.where(real_chrs ==  str(ec).replace("Chr", "").replace("chr", "") )[0]
            ret_echr_ix[ np.where(np.array( echr ) == ec)[0] ] = t_ix[0]
        return(ret_echr_ix)

    def estimated_cM_distance(self, snp_position):
        ## snp_position = "Chr1,150000" or "Chr1,1,300000"
        # Data based on
        #Salome, P. A., Bomblies, K., Fitz, J., Laitinen, R. A., Warthmann, N., Yant, L., & Weigel, D. (2011)
        #The recombination landscape in Arabidopsis thaliana F2 populations. Heredity, 108(4), 447-55.
        if "recomb_rates" in self.json.keys():
            mean_recomb_rates = self.json['recomb_rates']
        else:
            log.warn("Average recombination rates were missing in genome file. Add rates for each chromosome as an array in genome json file under 'recomb_rates' key. Using default rate of 3")
            mean_recomb_rates = np.repeat(3, len(genome.chrs_ids))
        assert isinstance(snp_position, basestring), "expected a string!"
        assert len(snp_position.split(",")) >= 2, "input should be 'chr1,1000' or 'chr1,1000,2000'"
        if len(snp_position.split(",")) == 2:
            snp_position = [snp_position.split(",")[0], int(snp_position.split(",")[1])]
        elif len(snp_position.split(",")) == 3:
            snp_position = [snp_position.split(",")[0], (int(snp_position.split(",")[1]) + int(snp_position.split(",")[2])) / 2 ]
        chr_ix = self.get_chr_ind( snp_position[0] )
        return( mean_recomb_rates[chr_ix] * snp_position[1] / 1000000 )


    def get_bins_genome(self, g, binLen):
        binLen = int(binLen)
        g_chrs_ids = np.char.replace(np.core.defchararray.lower(np.array(g.chrs, dtype="str")), "chr", "")
        common_chr_ids = np.intersect1d(g_chrs_ids, self.chrs_ids)
        assert len(g_chrs_ids) <= len(self.chrs_ids), "Please change default --genome option"
        assert len(common_chr_ids) > 0, "Please change default --genome option"
        if len(common_chr_ids) < len(self.chrs_ids):
            log.warn("Some reference contigs are missing in genotype hdf5 file")
        for chr_ix in range(len(self.chrs_ids)):
            t_g_ix = np.where(g_chrs_ids == self.chrs_ids[chr_ix])[0]
            if len(t_g_ix) == 0:
                chr_pos = np.zeros(0, dtype=int)
            else:
                start = g.chr_regions[t_g_ix[0]][0]
                end = g.chr_regions[t_g_ix[0]][1]
                chr_pos = g.positions[start:end]
            echr_bins = get_bins_echr(self.chrlen[chr_ix], chr_pos, binLen, start)
            for e_bin in echr_bins:
                yield((chr_ix, e_bin[0], e_bin[1]))

    def get_bins_arrays(self, g_chrs, g_snppos, binLen):
        g_chrs = np.char.replace(np.core.defchararray.lower(np.array(g_chrs, dtype="str")), "chr", "")
        g_chrs_ids = np.unique(g_chrs)
        common_chr_ids = np.intersect1d(g_chrs_ids, self.chrs_ids)
        assert len(g_chrs_ids) <= len(self.chrs_ids), "Please change default --genome option"
        assert len(common_chr_ids) > 0, "Please change default --genome option"
        if len(common_chr_ids) < len(self.chrs_ids):
            log.warn("Some reference contigs are missing in given SNPs")
        for chr_ix in range(len(self.chrs_ids)):
            chr_pos_ix = np.where(g_chrs == self.chrs_ids[chr_ix])[0]
            if len(chr_pos_ix) > 0:
                echr_bins = get_bins_echr(self.chrlen[chr_ix], g_snppos[chr_pos_ix], binLen, chr_pos_ix[0])
            else:
                echr_bins = get_bins_echr(self.chrlen[chr_ix], g_snppos[chr_pos_ix], binLen, 0)
            for e_bin in echr_bins:
                yield((chr_ix, e_bin[0], e_bin[1]))


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
