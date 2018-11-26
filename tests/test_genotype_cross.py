import pytest
import numpy as np
from snpmatch.core import snpmatch
from snpmatch.core import csmatch
from snpmatch.core import parsers


class TestSNPmatch:

    def test_vcf_parse(self, snps_vcf):
        assert len(snps_vcf.chrs) == 7545
        assert snps_vcf.chrs[0] == 'Chr1'
        assert snps_vcf.gt[0] == '0/0'

    def test_bed_parse(self, snps_bed):
        assert len(snps_bed.chrs) == 10000
        assert snps_bed.chrs[0] == '1'
        assert snps_bed.gt[0] == '0/0'
        assert snps_bed.pos[1] == 51103


    def test_likelihood(self, snp_numbers):
        assert snpmatch.likeliTest(snp_numbers[0], snp_numbers[1]) == 122.8361221819443
        
