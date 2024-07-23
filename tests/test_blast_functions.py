# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 13:12:57 2024

@author: u03132tk
"""
#from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Record import Blast, HSP, Alignment

#import sys
#sys.path.insert(0, 'D:\\')
from SyntenyQC.blast_functions import FastaWriter, get_best_hsp, results_to_hits, hits_to_best_hits, best_hits_to_rbh
from io import StringIO
import pytest
from general_mocks import mock_get_gbk_files, mock_read_good_gbk

#TODO check how BLAST+ handles no hits and hits to some but not all queries - assume they are just not reproted. 
 

def make_alignment(name, test_hsp_parameters):
    hit = Alignment()
    hit.hit_def = name
    hit.hsps = []
    for index, (identites, align_length, score) in enumerate(test_hsp_parameters):
        hsp = HSP()
        hsp.index = index
        hsp.identities, hsp.align_length, hsp.score = identites, align_length, score
        hit.hsps += [hsp]
    return hit

@pytest.mark.parametrize('min_percent_identity,expected_hsp_index',  [(30, 2),
                                                                      (75, 0), 
                                                                      (100, 3)
                                                                      ]
                         )
def test_get_best_hsp(min_percent_identity, expected_hsp_index):
    test_hsp_parameters = [(140, 150, 11), 
                           (140, 150, 11), 
                           (50, 150, 200), 
                           (150, 150, 10)]
    alignment = make_alignment('None', test_hsp_parameters)
    assert get_best_hsp(alignment, min_percent_identity).index == expected_hsp_index





def test_write_fasta(monkeypatch):   
    #https://stackoverflow.com/a/3945057/11357695
    monkeypatch.setattr('SyntenyQC.blast_functions.get_gbk_files', mock_get_gbk_files)
    monkeypatch.setattr('SyntenyQC.blast_functions.read_gbk', mock_read_good_gbk)
    outfile = StringIO()
    FastaWriter.write_fasta(genbank_folder = 'a_genbank_folder', 
                            outfile_handle = outfile)
    outfile.seek(0)
    content = outfile.read()
    assert content == '>file1__0\na protein1 sequence\n'+\
                      '>file1__1\na protein2 sequence\n'+\
                      '>file2__0\na protein1 sequence\n'+\
                      '>file2__1\na protein2 sequence'

def test_write_fasta_fail(monkeypatch):   
    #https://stackoverflow.com/a/3945057/11357695
    
    def mock_read_empty_translation_gbk(path):
        return SeqRecord(seq = None,
                         features = [SeqFeature(type="CDS",
                                                     qualifiers = {'translation' : [],
                                                                   }
                                                     )
                                     ]
                         )
    def mock_read_nonsense_translation(path):
        return SeqRecord(seq = None,
                         features = [SeqFeature(type="CDS",
                                                qualifiers = {'translation' : [''],
                                                              }
                                                )
                                     ]
                         )
    def mock_read_no_translation(path):
        return SeqRecord(seq = None,
                         features = [SeqFeature(type="CDS",
                                                qualifiers = {}
                                                )
                                     ]
                         )
    monkeypatch.setattr('SyntenyQC.blast_functions.get_gbk_files', mock_get_gbk_files)
    
    monkeypatch.setattr('SyntenyQC.blast_functions.read_gbk', mock_read_empty_translation_gbk)
    with pytest.raises(KeyError):       
        FastaWriter.write_fasta(genbank_folder = 'a_genbank_folder', 
                                outfile_handle = StringIO()
                                )
    monkeypatch.setattr('SyntenyQC.blast_functions.read_gbk', mock_read_nonsense_translation)
    with pytest.raises(KeyError):       
        FastaWriter.write_fasta(genbank_folder = 'a_genbank_folder', 
                                outfile_handle = StringIO()
                                )
    monkeypatch.setattr('SyntenyQC.blast_functions.read_gbk', mock_read_no_translation)
    with pytest.raises(KeyError):       
        FastaWriter.write_fasta(genbank_folder = 'a_genbank_folder', 
                                outfile_handle = StringIO()
                                )
        
@pytest.fixture
def expected_reconstructed_dict():
    return {'file1' : {'0' : {'file1' : {'0' : (150, 150, 100),
                                                  '1' : (140, 150, 50)
                                                  },
                                       'file2' : {'0' : (140, 150, 50),
                                                  '1' : (140, 150, 40)
                                                  }
                                       },
                                '1' : {'file1' : {'0' : (140, 150, 50),
                                                   '1' :(150, 150, 100)
                                                   },
                                        'file2' : {'0' : (140, 150, 50),
                                                   '1' : (140, 150, 50)
                                                   }
                                        }
                                },
                     'file2' : {'0' : {'file1' : {'0' : (140, 150, 50),
                                                  '1' :(140, 150, 50)
                                                  },
                                       'file2' : {'0' : (150, 150, 100),
                                                  '1' : (140, 150, 50)
                                                  }
                                       },
                                '1' : {'file1' : {'0' : (140, 150, 50),
                                                  '1' :(140, 150, 50)
                                                  },
                                       'file2' : {'0' : (140, 150, 50),
                                                  '1' : (150, 150, 100)
                                                  }
                                       }
                                }
                     }

#TODO add some NONE hits to check they get ignored
@pytest.fixture
def expected_best_hit_matrix():
    return {'file1' : {'0' : {'file1' : '0',
                              'file2' : '0'},
                       '1' : {'file1' : '1',
                              'file2' : '0'},
                       },

            'file2' : {'0' : {'file1' : '0',
                              'file2' : '0'},
                       '1' : {'file1' : '0',
                              'file2' : '1'}
                              }
            }
@pytest.fixture
def expected_reciprocal_best_hit_matrix():
    return {'file1' : {'0' : {'file1' : '0',
                              'file2' : '0'},
                       '1' : {'file1' : '1'},
                       },

            'file2' : {'0' : {'file1' : '0',
                              'file2' : '0'},
                       '1' : {'file2' : '1'}
                       }
            }
def mock_read_xml(path):
    #deal with missing hits etc
    alignment_parameters = {'file1__0' : {'file1__0' : [(150, 150, 100)],
                                          'file1__1' : [(140, 150, 50)], 
                                          'file2__0' : [(140, 150, 50)],
                                          'file2__1' : [(140, 150, 40)]}, 
                            'file1__1' : {'file1__0' : [(140, 150, 50)],
                                          'file1__1' : [(150, 150, 100)], 
                                          'file2__0' : [(140, 150, 50)],
                                          'file2__1' : [(140, 150, 50)]},
                            'file2__0' : {'file1__0' : [(140, 150, 50)],
                                          'file1__1' : [(140, 150, 50)], 
                                          'file2__0' : [(150, 150, 100)],
                                          'file2__1' : [(140, 150, 50)]},
                            'file2__1' : {'file1__0' : [(140, 150, 50)],
                                          'file1__1' : [(140, 150, 50)], 
                                          'file2__0' : [(140, 150, 50)],
                                          'file2__1' : [(150, 150, 100)]}
                            }
    records = []
    for query, hits in alignment_parameters.items():
        record = Blast()
        record.query = query
        record.alignments = []
        for name, hsp_parameters in hits.items():
            record.alignments += [make_alignment(name, hsp_parameters)]
        records += [record]
    return records

def test_results_to_hits(expected_reconstructed_dict, monkeypatch):
    monkeypatch.setattr('SyntenyQC.blast_functions.read_xml', mock_read_xml)
    hit_matrix =  results_to_hits('a path', min_percent_identity = 0) 
    #no __eq__ method for hsps - so reconstruct with attrs of unterest
    reconstructed_dict = {}
    for query_scaffold, queries in hit_matrix.items():
        reconstructed_dict[query_scaffold] = {}
        for query_index, all_hits in queries.items():
            reconstructed_dict[query_scaffold][query_index] = {}
            for hit_scaffold, hits in all_hits.items():
                reconstructed_dict[query_scaffold][query_index][hit_scaffold] = {}
                for hit_index, hit_hsp in hits.items():
                    hsp_details = (hit_hsp.identities, hit_hsp.align_length, hit_hsp.score)
                    reconstructed_dict[query_scaffold][query_index][hit_scaffold][hit_index] = hsp_details
    assert reconstructed_dict == expected_reconstructed_dict
    

    
    


def test_hits_to_best_hits(expected_best_hit_matrix, monkeypatch):
    monkeypatch.setattr('SyntenyQC.blast_functions.read_xml', mock_read_xml)
    hit_matrix =  results_to_hits('a path', min_percent_identity = 0) 
    best_hit_matrix = hits_to_best_hits(hit_matrix)
    assert best_hit_matrix == expected_best_hit_matrix
    
def test_best_hits_to_rbh(expected_reciprocal_best_hit_matrix, monkeypatch):
    monkeypatch.setattr('SyntenyQC.blast_functions.read_xml', mock_read_xml)
    hit_matrix =  results_to_hits('a path', min_percent_identity = 0) 
    best_hit_matrix = hits_to_best_hits(hit_matrix)
    reciprocal_best_hit_matrix = best_hits_to_rbh(best_hit_matrix)
    assert reciprocal_best_hit_matrix == expected_reciprocal_best_hit_matrix
