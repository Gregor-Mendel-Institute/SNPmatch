"""
Class function to read in output csv files from SNPmatch
"""

import pandas as pd
import numpy as np
import os.path
from . import snpmatch


class FollowSNPmatch(object):

    def __init__(self, csv_snpmatch = {}, csv_csmatch = {}):
        #  **kwargs
        '''
        Class function to read in output csv from SNPmatch (intermediate_modified.csv)
        use kwargs to add in multiple csv for different databases (1001g, 250k etc. )
        '''
        self._instances = []
        if csv_snpmatch != {}:
            for req_name in csv_snpmatch.keys():
                req_csv = pd.read_csv(csv_snpmatch[req_name], sep = None, engine = 'python', index_col=0)
                self._instances.append( "snpmatch_" + req_name )
                ## Below we change the column datatype
                for ef in req_csv.columns.intersection( ['TopHitAccession', 'NextHit', 'ThirdHit', 'RefinedTopHit'] ):
                    req_csv[ef] = req_csv[ef].apply(str)
                setattr(self, "snpmatch_" + req_name, req_csv)
                setattr(self, "snpmatch_" + req_name + "_fol", os.path.dirname( csv_snpmatch[req_name] ) )
        if csv_csmatch != {}:
            for req_name in csv_csmatch.keys():
                req_csv = pd.read_csv(csv_csmatch[req_name], sep = None, engine = 'python', index_col=0)
                self._instances.append( "csmatch_" + req_name )
                ## Below we change the column datatype
                for ef in ['TopHit', 'NextHit', 'ObservedParent1','ObservedParent2']:
                    req_csv[ef] = req_csv[ef].apply(str)
                setattr(self, "csmatch_" + req_name, req_csv)
                setattr(self, "csmatch_" + req_name + "_fol", os.path.dirname( csv_csmatch[req_name] ) )
        self._instances = pd.Series(self._instances)
    
    def beauty_print(self, req_name, req_ix = None, req_index = None):
        '''
        Simple function to only print required columns.
        '''
        req_beauty_print = pd.Series([
            "TopHitAccession", 
            "NextHit",
            "ThirdHit",
            "Score",
            "FracScore",
            'identity',
            "SNPsinfoAcc",
            "TopHitsNumber",
            "dist_to_tophit",
            'percent_heterozygosity',
            "IdenticalWindows",
            "RefinedTopHit",
            'RefinedTopHitNumber'
        ])
        common_cols = req_beauty_print[req_beauty_print.isin( self.__getattribute__(req_name).columns )]
        # if req_index:
        #     return(self.__getattribute__(req_name).loc[ req_ix, common_cols ])
        if req_ix is not None:
            return(self.__getattribute__(req_name).loc[ req_ix, common_cols ])
        return(self.__getattribute__(req_name).loc[:, common_cols ])
        


    def get_identity(self, req_name = None, error_rate=0.02):
        '''
        Function to determine whether sample is identical to TopHit
        '''
        if req_name is None:
            for req_name in self._instances[self._instances.str.contains( "snpmatch_" )]:
                self.__getattribute__(req_name)['identity'] = snpmatch.np_test_identity(
                    x = self.__getattribute__(req_name)['Score'] * self.__getattribute__(req_name)['SNPsinfoAcc'], 
                    n = self.__getattribute__(req_name)['SNPsinfoAcc'], 
                    error_rate = error_rate
                )
        else:
            assert req_name in self._instances, "provided %s is not present in instances" % req_name
            self.__getattribute__(req_name)['identity'] = snpmatch.np_test_identity(
                x = self.__getattribute__(req_name)['Score'] * self.__getattribute__(req_name)['SNPsinfoAcc'], 
                n = self.__getattribute__(req_name)['SNPsinfoAcc'], 
                error_rate = error_rate
            )





