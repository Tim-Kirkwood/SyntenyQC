# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 17:24:45 2024

@author: u03132tk
"""
#import sys
#sys.path.insert(0, 'D:\\')
from SyntenyQC.helpers import get_gbk_files, get_cds_count, check_motif, get_protein_seq
#from Bio.SeqFeature import SeqFeature
#from Bio.SeqRecord import SeqRecord
import pytest
from Bio.SeqFeature import SeqFeature 
from general_mocks import mock_get_gbk_files, mock_read_good_gbk, mock_listdir
#from Bio.SeqFeature import SeqFeature
#from Bio.SeqRecord import SeqRecord

def test_get_gbk_files(monkeypatch):
    def mock_isfile(path):
        if '.' in path:
            defined_filetypes = ['.gbk', '.gb', '.txt']
            return path[path.rindex('.') :] in defined_filetypes
        else:
            #no suffix
            return False    
    def mock_empty_listdir(path):
        return []
    def mock_no_gbk(path):
        return ["accession.txt", "organism.txt", 'file3.txt', 'file3.txt', 'a_directory']
    monkeypatch.setattr('os.listdir', mock_listdir)
    monkeypatch.setattr('os.path.isfile', mock_isfile)
    assert get_gbk_files('a folder') == ["accession.gbk", "organism.gbk", 'file3.gb']
    monkeypatch.setattr('os.listdir', mock_empty_listdir)
    assert get_gbk_files('a folder') == []
    monkeypatch.setattr('os.listdir', mock_no_gbk)
    assert get_gbk_files('a folder') == []

def test_get_cds_count(monkeypatch):
    monkeypatch.setattr('SyntenyQC.helpers.get_gbk_files', mock_get_gbk_files)
    monkeypatch.setattr('SyntenyQC.helpers.read_gbk', mock_read_good_gbk)
    assert get_cds_count('a_folder') == 10

@pytest.mark.parametrize('motif_start,motif_end,expected_msg',  
                         [(-1, 1, 'motif_start -1 < 0\n'),
                          (-2, -1, 'motif_start -2 < 0\nmotif_stop -1 <= 0\n'),
                          (-2, -1, 'motif_start -2 < 0\nmotif_stop -1 <= 0\n'),
                          (1, 1,  'motif_start 1 is >= motif_stop 1\n'),
                          (1, -1, 'motif_stop -1 <= 0\nmotif_start 1 is >= motif_stop -1\n'),
                          (1, 0, 'motif_stop 0 <= 0\nmotif_start 1 is >= motif_stop 0\n'),
                          (1, 2, ''),
                          (0, 2, '')
                          ]
                         )
def test_check_motif(motif_start, motif_end, expected_msg):
    assert check_motif(motif_start, motif_end) == expected_msg

def test_get_protein_seq():
    with pytest.raises(KeyError):
        get_protein_seq(SeqFeature(type="CDS",
                                   qualifiers = {'translation' : [],
                                                 }
                                   )
                        )
    with pytest.raises(KeyError):
        get_protein_seq(SeqFeature(type="CDS",
                                   qualifiers = {'translation' : [''],
                                                 }
                                   )
                        )
    with pytest.raises(KeyError):
        get_protein_seq(SeqFeature(type="CDS",
                                   qualifiers = {}
                                   )
                        )
    assert '' == get_protein_seq(SeqFeature(type="CDS",
                                            qualifiers = {'translation' : [],
                                                          'pseudo' : None}
                                            )
                                 )   
    assert '' == get_protein_seq(SeqFeature(type="CDS",
                                            qualifiers = {'translation' : [''],
                                                          'pseudo' : None}
                                            )
                                 )   
    assert '' == get_protein_seq(SeqFeature(type="CDS",
                                            qualifiers = {'pseudo' : None}
                                            )
                                 )   