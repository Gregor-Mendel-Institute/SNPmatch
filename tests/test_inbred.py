import pytest
import numpy as np
from snpmatch.core import snpmatch
from snpmatch.core import parsers


class TestSNPmatch:

    def test_vcf_parse(self, snps_vcf):
        assert len(snps_vcf.chrs) == 7545
        assert snps_vcf.chrs[0] == 'Chr1'
        assert snps_vcf.gt[0] == '0/0'


    def test_likelihood(self, snp_numbers):
        assert snpmatch.likeliTest(snp_numbers[0], snp_numbers[1]) == 122.8361221819443
