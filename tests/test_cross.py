import pytest
from snpmatch.core import csmatch
import numpy as np

class TestCross:

    def test_number_of_chromomes(self, geno):
        chrs = geno.chrs
        assert len(chrs) == len(csmatch.tair_chrs)
