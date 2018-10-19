import pytest
from os import path
import json

def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true",help="run slow tests")

def pytest_runtest_setup(item):
    if 'slow' in item.keywords and not item.config.getoption("--runslow"):
        pytest.skip("need --runslow option to run")

slow = pytest.mark.slow
resource_path = path.join(path.dirname(__file__), '..', 'sample_files')
num_lines = 10
num_matched = 3

@pytest.fixture
def geno():
    from pygwas.core import genotype
    return(genotype.load_hdf5_genotype_data('%s/all_chromosomes_binary.hdf5' %resource_path))

@pytest.fixture
def geno_acc():
    from pygwas.core import genotype
    return(genotype.load_hdf5_genotype_data('%s/all_chromosomes_binary.acc.hdf5' %resource_path))

@pytest.fixture
def snps_vcf():
    from snpmatch.core import parsers
    return(parsers.ParseInputs(inFile = '%s/701_501.filter.vcf' % resource_path, logDebug = True))

@pytest.fixture
def snps_bed():
    from snpmatch.core import parsers
    return(parsers.ParseInputs(inFile = '%s/701_502.filter.bed' % resource_path, logDebug = True))

@pytest.fixture
def snp_numbers():
    return((num_lines, num_matched))
