from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='SNPmatch',
    version='1.6.1',
    description='A tool to get maximum likely accession in database',
    long_description=long_description,
    url='https://github.com/Gregor-Mendel-Institute/snpmatch',
    author=['Rahul Pisupati'],
    author_email='rahul.bharadwaj.p@gmail.com',
    license='GMI',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.4',
    ],
    keywords='Genotyping Low Coverage sequencing data',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=[
        "scipy >= 0.17.0",
        "numpy >=1.9.0",
        "PyGWAS",
        "vcfnp",
        "pandas"
    ],
    entry_points={
        'console_scripts': [
            'snpmatch=snpmatch:main'
        ],
    },
)

