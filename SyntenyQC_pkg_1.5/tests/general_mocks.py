# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 11:32:25 2024

@author: u03132tk
"""
#from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


def mock_read_good_gbk(path : str) -> SeqRecord:
    test_features = [SeqFeature(type="CDS",
                                qualifiers = {'translation' : ['a ', 
                                                               'protein1 ', 
                                                               'sequence'],
                                              'pseudo' : None}
                                ),
                     SeqFeature(type="CDS",
                                qualifiers = {'translation' : ['a ', 
                                                               'protein2 ', 
                                                               'sequence']
                                              }
                                ),
                     SeqFeature(type="CDS",
                                qualifiers = {'pseudo' : None,
                                              'translation' : []
                                              }
                                ),
                     SeqFeature(type="CDS",
                                qualifiers = {'pseudo' : None,
                                              'translation' : ['']
                                              }
                                ),
                     SeqFeature(type="CDS",
                                qualifiers = {'pseudo' : None}
                                ),
                     SeqFeature(type="not_CDS",
                                qualifiers = {'translation' : ['not ',
                                                               'a ', 
                                                               'protein2 ', 
                                                               'sequence']}
                                )
                     ]
    return SeqRecord(seq = None,
                     features = test_features)

def mock_makedirs(path : str):
    pass

def mock_get_gbk_files(folder : str) -> list:
    return ["file1.gbk", "file2.gb"]

def mock_listdir(folder : str) -> list:
    '''
    For test simplicity/predictability, replace os.listdir() with a mocked 
    function that returns the value of mock_return when called in 
    test_make_filepath().
    

    Parameters
    ----------
    folder : str
        Folder supplied to os.listdir().

    Returns
    -------
    list
        List of non-existent filenames for testing purposes.
    '''
    return ["accession.gbk", "organism.gbk", "a_directory(1)", 'file3.gb', 'file3.txt', 'a_directory']

