# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 13:12:57 2024

@author: u03132tk
"""
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Record import Blast, HSP, Alignment
from SyntenyQC.blast_functions import FastaWriter, get_best_hsp, results_to_hits, hits_to_best_hits, best_hits_to_rbh, make_rbh_matrix 
from io import StringIO
import pytest
from general_mocks import mock_get_gbk_files, mock_read_good_gbk

# =============================================================================
# TEST GAPS
#
# Difficult to test external tools/system-level integrations:
#   TODO - look at testing these at some point (see 
#          https://github.com/gamcil/cblaster/blob/master/tests/test_local.py)
# - run_blast_process
# - makeblastdb_subprocess
# - blastP_subprocess
# - all_vs_all_blast, 
#       
# 
# Thin wrapper of external lib:
# - read_xml
# =============================================================================







# =============================================================================
# helper function for fixtures fake_alignment() and fake_xml_results()
# =============================================================================

def make_alignment(name : str, test_hsp_parameters : list) -> Alignment:
    '''
    Make a fake Alignment with one HSP per tuple in test_hsp_parameters.

    Parameters
    ----------
    name : str
        Hit name.
    test_hsp_parameters : dict
        List of tuples, one tuple per hsp (used to define hsp identities, 
        alignment length and blastp score).

    Returns
    -------
    hit : Alignment
        A biopython BLASTP alignment between query and one hit called 'name', 
        with one or more biopython HSP objects.  Each HSP also has an index 
        attribute, specific to this test, which is used to identify which HSP 
        is returned by get_best_hsp.

    '''
    hit = Alignment()
    hit.hit_def = name
    hit.hsps = []
    for index, (identites, align_length, score) in enumerate(test_hsp_parameters):
        hsp = HSP()
        hsp.index = index
        hsp.identities, hsp.align_length, hsp.score = identites, align_length, score
        hit.hsps += [hsp]
    return hit




# =============================================================================
# TestClasses
# =============================================================================

class TestGetBestHsp:
    
    @pytest.fixture
    def fake_alignment(self) -> Alignment:
        test_hsp_parameters = [(140, 150, 11), #identites, align_length, score
                               (140, 150, 11), 
                               (50, 150, 200), 
                               (150, 150, 10)]
        alignment = make_alignment('doesnt_matter', 
                                   test_hsp_parameters)
        return alignment
    
    @pytest.mark.parametrize('min_percent_identity,expected_hsp_index',  
                             [(30, 2),
                              (75, 0), #takes first hsp that meets criteria
                              (100, 3)
                              ]
                             )
    def test_with_hsp(self,
                      #function args
                      min_percent_identity : int, fake_alignment : Alignment,
                      #test
                      expected_hsp_index : int):
        best_hsp = get_best_hsp(fake_alignment, 
                                min_percent_identity)
        assert best_hsp.index == expected_hsp_index
        
    def test_without_hsp(self, 
                         #function args
                         fake_alignment : Alignment):
        best_hsp = get_best_hsp(fake_alignment, 
                                min_percent_identity = 101)
        assert best_hsp == None



class TestWriteFasta:
    #https://stackoverflow.com/a/3945057/11357695
    
    
# =============================================================================
#     used by all setups
# =============================================================================

    @pytest.fixture 
    def setup_get_gbk_files(self, monkeypatch):
        monkeypatch.setattr('SyntenyQC.blast_functions.get_gbk_files', 
                            mock_get_gbk_files)
        
    
# =============================================================================
#     Setups for different SeqRecords
# =============================================================================

    @pytest.fixture
    def setup_full_run(self, monkeypatch, setup_get_gbk_files):
        monkeypatch.setattr('SyntenyQC.blast_functions.read_gbk', 
                            mock_read_good_gbk)
    
    @pytest.fixture
    def setup_no_cds(self, monkeypatch, setup_get_gbk_files):
        def mock_read_no_cds_gbk(path : str) -> SeqRecord:
            feature = SeqFeature(type="not_CDS",
                                 qualifiers = {'translation' : ['not_a_protein'],
                                               }
                                 )
            return SeqRecord(seq = None,
                             features = [feature]
                             )
        monkeypatch.setattr('SyntenyQC.blast_functions.read_gbk',
                            mock_read_no_cds_gbk)
    
    @pytest.fixture
    def setup_empty_translation(self, monkeypatch, setup_get_gbk_files):
        def mock_read_empty_translation_gbk(path : str) -> SeqRecord:
            feature = SeqFeature(type="CDS",
                                 qualifiers = {'translation' : [],
                                               }
                                 )
            return SeqRecord(seq = None,
                             features = [feature]
                             )
        monkeypatch.setattr('SyntenyQC.blast_functions.read_gbk', 
                            mock_read_empty_translation_gbk)
    
    @pytest.fixture
    def setup_nonsense_translation(self, monkeypatch, setup_get_gbk_files):
        def mock_read_nonsense_translation(path : str) -> SeqRecord:
            feature = SeqFeature(type="CDS",
                                 qualifiers = {'translation' : [''],
                                               }
                                 )
            return SeqRecord(seq = None,
                             features = [feature]
                             )
        monkeypatch.setattr('SyntenyQC.blast_functions.read_gbk', 
                            mock_read_nonsense_translation)
    
    @pytest.fixture
    def setup_no_translation(self, monkeypatch, setup_get_gbk_files):
        def mock_read_no_translation(path : str) -> SeqRecord:
            feature = SeqFeature(type="CDS",
                                 qualifiers = {}
                                 )
            return SeqRecord(seq = None,
                             features = [feature]
                             )
        monkeypatch.setattr('SyntenyQC.blast_functions.read_gbk', 
                            mock_read_no_translation)
        
        
# =============================================================================
#     Tests
# =============================================================================

    @pytest.mark.parametrize('fixture', 
                             ['setup_empty_translation', 
                              'setup_nonsense_translation', 
                              'setup_no_translation']
                             )    
    def test_KeyError(self, fixture, request):
        #cannot parametise fixtures directly (they are input as functions rather than 
        #being called first)
        #suggestion for requests: 
        #   https://github.com/pytest-dev/pytest/issues/349#issuecomment-112203541
        #requests docs: 
        #   https://docs.pytest.org/en/7.1.x/reference/reference.html#pytest.FixtureRequest.getfixturevalue
        request.getfixturevalue(fixture)
        with pytest.raises(KeyError):       
            FastaWriter.write_fasta(genbank_folder = 'a_genbank_folder', 
                                    outfile_handle = StringIO()
                                    )
            
    def test_ValueError(self, setup_no_cds):
        with pytest.raises(ValueError):       
            FastaWriter.write_fasta(genbank_folder = 'a_genbank_folder', 
                                    outfile_handle = StringIO()
                                    )
    
    def test_full(self, setup_full_run):   
        outfile = StringIO() #instead of writing to file
        FastaWriter.write_fasta(genbank_folder = 'a_genbank_folder', 
                                outfile_handle = outfile)
        outfile.seek(0)
        content = outfile.read()
        assert content == '>file1__0\na protein1 sequence\n'\
                          '>file1__1\na protein2 sequence\n'\
                          '>file2__0\na protein1 sequence\n'\
                          '>file2__1\na protein2 sequence'

        


class TestMakeRbhMatrix:
    
# =============================================================================
#     test fixtures
# =============================================================================
    @pytest.fixture
    def expected_reconstructed_dict(self) -> dict:
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
    
    @pytest.fixture
    def expected_best_hit_matrix(self) -> dict:
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
    def expected_reciprocal_best_hit_matrix(self) -> dict:
        return {'file1' : {'0' : {'file1' : '0',
                                  'file2' : '0'},
                           '1' : {'file1' : '1'},
                           },

                'file2' : {'0' : {'file1' : '0',
                                  'file2' : '0'},
                           '1' : {'file2' : '1'}
                           }
                }
    
    @pytest.fixture 
    def fake_xml_results(self) -> list:
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
                record.alignments += [make_alignment(name, 
                                                     hsp_parameters)]
            records += [record]
        return records
    
# =============================================================================
#     setup fixture
# =============================================================================
    @pytest.fixture
    def setup_fake_xml(self, monkeypatch, fake_xml_results):
        def mock_read_xml(path : str) -> list:
            return fake_xml_results
        monkeypatch.setattr('SyntenyQC.blast_functions.read_xml', 
                            mock_read_xml)
        
        
# =============================================================================
#     Helper function    
# =============================================================================
    def reconstruct_dict(self, hit_matrix : dict) -> dict:
        '''
        SyntenyQC.blast_functions.results_to_hits() returns a dict mapping 
        querys -> hits -> hsps.  Bio.Blast.Record.HSP has no __eq__ method so I 
        can't comapre the HSPs generated vs expected.  Instead, I pull 
        attributes of interest for each HSP (identities, alignment length and 
        blastp score) and check they match what I expect.

        Parameters
        ----------
        hit_matrix : dict
            Output of SyntenyQC.blast_functions.results_to_hits().

        Returns
        -------
        dict
            Dictionary of hits, of the format:
                
                {query_file_A : 
                    {query_protein_1 : 
                        {hit_file_A : 
                            {hit_protein_1A : (identities, align_length, score), #from best HSP
                             hit_protein_2A : (identities, align_length, score), #from best HSP
                             ...
                             },
                         hit_file_B : 
                             {hit_protein_1B : (identities, align_length, score), #from best HSP
                              ...
                             }
                        },
                    query_protein_2 : {..................................},
                    ...
                    },
                query_file_B : {.........................................},
                ...
                }

        '''
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
        return reconstructed_dict
            
# =============================================================================
#     tests
# =============================================================================
    def test_unit(self,
                  #test comparison
                  expected_reconstructed_dict : dict,  
                  expected_best_hit_matrix : dict, 
                  expected_reciprocal_best_hit_matrix : dict, 
                  #test setup
                  fake_xml_results : list):
        #fake_xml_results = read_xml('xml_path')
        hit_matrix =  results_to_hits(fake_xml_results, 
                                      min_percent_identity = 0) 
        #no __eq__ method for hsps - so reconstruct with attrs of unterest
        assert self.reconstruct_dict(hit_matrix) == expected_reconstructed_dict
        
        best_hit_matrix = hits_to_best_hits(hit_matrix)
        assert best_hit_matrix == expected_best_hit_matrix
        
        reciprocal_best_hit_matrix = best_hits_to_rbh(best_hit_matrix)
        assert reciprocal_best_hit_matrix == expected_reciprocal_best_hit_matrix
    
    def test_run(self, 
                 #test setup
                 setup_fake_xml, 
                 #test comparison
                 expected_reciprocal_best_hit_matrix : dict):
        rbh_matrix =  make_rbh_matrix(xml_path = 'xml_path', 
                                      min_percent_identity = 0) 
        assert rbh_matrix == expected_reciprocal_best_hit_matrix