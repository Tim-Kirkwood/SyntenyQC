# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 17:24:45 2024

@author: u03132tk
"""
from SyntenyQC.helpers import get_gbk_files, get_cds_count, check_motif, \
                              get_protein_seq, initalise_log, log_params, \
                              ingest_binary_csv
import pytest
from Bio.SeqFeature import SeqFeature 
from Bio.SeqRecord import SeqRecord
from general_mocks import mock_read_good_gbk, mock_listdir
import logging
import pandas as pd

#not tested - read_gbk (thin wrapper so no point) 
#not tested - TestGetLogHandlers (not sure how to mock the file handler)

# =============================================================================
# File/log handling fixtures - monkeypatch setup 
# =============================================================================



@pytest.fixture
def initialised_log(log_setup) -> logging.Logger:
    logger = initalise_log(log_file = 'dosent matter', 
                           logger_name = 'test')
    return logger

@pytest.fixture 
def setup_is_file(monkeypatch):
    def mock_isfile(path : str) -> bool:
        if '.' in path:
            #filetypes from mock_listdir
            defined_filetypes = ['.gbk', '.gb', '.txt']
            suffix = path[path.rindex('.') :]
            #in case test changes in future
            assert suffix in defined_filetypes
            return True
        else:
            #no suffix - dir
            return False   
    monkeypatch.setattr('os.path.isfile', 
                        mock_isfile)
    
@pytest.fixture 
def setup_folder_with_gbks(monkeypatch, setup_is_file):
    monkeypatch.setattr('os.listdir', 
                        mock_listdir)
    
@pytest.fixture   
def setup_folder_without_files(monkeypatch, setup_is_file):
    #setup_is_file should not be needed but is included for uniform test 
    #conditions 
    def mock_empty_listdir(path : str) -> list:
        return []
    monkeypatch.setattr('os.listdir', 
                        mock_empty_listdir)
    
@pytest.fixture 
def setup_folder_without_gbks(monkeypatch, setup_is_file):
    def mock_no_gbk(path : str) -> list:
        return ["accession.txt", "organism.txt", 'file3.txt', 'file3.txt', 'a_directory']
    monkeypatch.setattr('os.listdir', 
                        mock_no_gbk)

# =============================================================================
# SeqFeature/SeqRecord fixtures - objects
# =============================================================================

@pytest.fixture 
def a_cds() -> SeqFeature:
    return SeqFeature(type="CDS",
               qualifiers = {'translation' : ['a_protetin_seq'],
                             }
               )

@pytest.fixture 
def not_a_cds() -> SeqFeature:
    return SeqFeature(type="not_a_CDS",
               qualifiers = {'translation' : ['not_a_protetin_seq'],
                             }
               )

@pytest.fixture 
def cds_no_translation() -> SeqFeature:
    return SeqFeature(type="CDS",
               qualifiers = {'translation' : [],
                             }
               )
    
@pytest.fixture 
def cds_empty_translation() -> SeqFeature:
    return SeqFeature(type="CDS",
               qualifiers = {'translation' : [''],
                             }
               )

@pytest.fixture 
def cds_no_seq() -> SeqFeature:
    return SeqFeature(type="CDS",
               qualifiers = {}
               )

@pytest.fixture 
def cds_pseudo_no_translation() -> SeqFeature:
    return SeqFeature(type="CDS",
               qualifiers = {'translation' : [],
                             'pseudo' : None}
               )

@pytest.fixture 
def cds_pseudo_empty_translation() -> SeqFeature:
    return SeqFeature(type="CDS",
               qualifiers = {'translation' : [''],
                             'pseudo' : None}
               )
                             
@pytest.fixture 
def cds_pseudo_no_seq() -> SeqFeature:
    return SeqFeature(type="CDS",
               qualifiers = {'pseudo' : None}
               )     


# =============================================================================
# tests
# =============================================================================


class TestGetGbkFiles:

    def test_normal(self, setup_folder_with_gbks):
        assert get_gbk_files('a folder') == ["accession.gbk", "organism.gbk", 'file3.gb']
    
    def test_empty_folder(self, setup_folder_without_files):
        assert get_gbk_files('a folder') == []
    
    def test_no_gbks(self, setup_folder_without_gbks):
        assert get_gbk_files('a folder') == []

class TestGetCdsCount:
    
    # =============================================================================
    # SeqFeature/SeqRecord fixtures - monkeypatch setup
    # =============================================================================
    @pytest.fixture
    def setup_read_gbk_without_cds(self, monkeypatch, setup_folder_with_gbks, 
                                   not_a_cds : SeqFeature):
        def mock_read_bad_gbk(path : str) -> SeqRecord:
            return SeqRecord(seq = None,
                             features = [not_a_cds])
        monkeypatch.setattr('SyntenyQC.helpers.read_gbk', 
                            mock_read_bad_gbk)

    @pytest.fixture
    def setup_read_gbk_no_translation(self, monkeypatch, setup_folder_with_gbks, 
                                      cds_no_translation : SeqFeature):
        def mock_read_bad_gbk(path : str) -> SeqRecord:
            return SeqRecord(seq = None,
                             features = [cds_no_translation])
        monkeypatch.setattr('SyntenyQC.helpers.read_gbk', 
                            mock_read_bad_gbk)

    @pytest.fixture
    def setup_read_gbk_empty_translation(self, monkeypatch, setup_folder_with_gbks, 
                                         cds_empty_translation : SeqFeature):
        def mock_read_bad_gbk(path : str) -> SeqRecord:
            return SeqRecord(seq = None,
                             features = [cds_empty_translation])
        monkeypatch.setattr('SyntenyQC.helpers.read_gbk', 
                            mock_read_bad_gbk)
        
    @pytest.fixture
    def setup_read_gbk_no_seq(self, monkeypatch, setup_folder_with_gbks, 
                              cds_no_seq : SeqFeature):
        def mock_read_bad_gbk(path : str) -> SeqRecord:
            return SeqRecord(seq = None,
                             features = [cds_no_seq])
        monkeypatch.setattr('SyntenyQC.helpers.read_gbk', 
                            mock_read_bad_gbk)
        
    @pytest.fixture
    def setup_read_gbk(self, monkeypatch, setup_folder_with_gbks):
        monkeypatch.setattr('SyntenyQC.helpers.read_gbk', 
                            mock_read_good_gbk)
        
    # =============================================================================
    #     tests - normal
    # =============================================================================
    def test_empty_folder(self, setup_folder_without_files):
        assert get_cds_count('a_folder') == 0
    def test_no_gbk(self, setup_folder_without_gbks):
        assert get_cds_count('a_folder') == 0
    def test_no_cds(self, setup_read_gbk_without_cds):
        assert get_cds_count('a_folder') == 0
    def test_gbk(self, setup_read_gbk):
        assert get_cds_count('a_folder') == 6
    
    # =============================================================================
    #     tests - exceptions
    # =============================================================================
    def test_no_translation(self, setup_read_gbk_no_translation):
        with pytest.raises(KeyError):
            get_cds_count('a_folder')
    def test_empty_translation(self, setup_read_gbk_empty_translation):
        with pytest.raises(KeyError):
            get_cds_count('a_folder')
    def test_no_seq(self, setup_read_gbk_no_seq):
        with pytest.raises(KeyError):
            get_cds_count('a_folder')
        

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
def test_check_motif(motif_start : int, motif_end : int, expected_msg : str):
    assert check_motif(motif_start, motif_end) == expected_msg

class TestGetProteinSeq:
    # =============================================================================
    #     normal run
    # =============================================================================
    def test_cds(self, a_cds : SeqFeature):
        assert get_protein_seq(a_cds) == 'a_protetin_seq'
    
    # =============================================================================
    #     pseudo proteins
    # =============================================================================
    def test_cds_pseudo_no_translation(self, cds_pseudo_no_translation : SeqFeature):
        assert get_protein_seq(cds_pseudo_no_translation) == ''
    def test_cds_pseudo_empty_translation(self, cds_pseudo_empty_translation : SeqFeature):
        assert get_protein_seq(cds_pseudo_empty_translation) == ''
    def test_cds_pseudo_no_seq(self, cds_pseudo_no_seq : SeqFeature):
        assert get_protein_seq(cds_pseudo_no_seq) == ''
    
    # =============================================================================
    #     bad format features
    # =============================================================================
    def test_not_a_cds(self, not_a_cds : SeqFeature):
        with pytest.raises(ValueError):
            get_protein_seq(not_a_cds)
    def test_cds_no_translation(self, cds_no_translation : SeqFeature):
        with pytest.raises(KeyError):
            get_protein_seq(cds_no_translation)
    def test_cds_empty_translation(self, cds_empty_translation : SeqFeature):
        with pytest.raises(KeyError):
            get_protein_seq(cds_empty_translation)
    def test_cds_no_seq(self, cds_no_seq : SeqFeature):
        with pytest.raises(KeyError):
            get_protein_seq(cds_no_seq)
    
            

    
    
class TestLogParams:
    
    @pytest.fixture
    def collect_vars(self) -> dict:
        return {"binary_path" : "bp",
                "write_genomes" : "wg",
                "strict_span" : "sp",
                "email" : "e",
                "neighbourhood_size" : "n",
                "filenames" : "f",
                "results_dir" : "rd"}
    
    @pytest.fixture
    def bad_collect_vars(self) -> dict:
        return {"binary_path" : "bp",
                "write_genomes" : "wg",
                "strict_span" : "sp",
                "email" : "e",
                "neighbourhood_size" : "n",
                "filenames" : "f",
                #"results_dir" : "rd"
                }
    
    @pytest.fixture
    def sieve_vars(self) -> dict:
        return {"input_genbank_dir" : "igd",
                "similarity_filter" : "sf",
                "output_genbank_dir" : "ogd",
                "e_value" : "ev",
                "results_dir" : "rd",
                "output_vis_dir" : "ovd",
                "min_percent_identity" : "mpi",
                "output_blast_dir" : "obd",
                "min_edge_view" : "mev"}
    
    @pytest.fixture
    def bad_sieve_vars(self) -> dict:
        return {"input_genbank_dir" : "igd",
                "similarity_filter" : "sf",
                "output_genbank_dir" : "ogd",
                "e_value" : "ev",
                "results_dir" : "rd",
                "output_vis_dir" : "ovd",
                "min_percent_identity" : "mpi",
                "output_blast_dir" : "obd",
                #"min_edge_view" : "mev",
                }
        
    def test_collect(self, initialised_log : logging.Logger, collect_vars : dict, 
                     caplog):
        with caplog.at_level(logging.INFO):
            log_params(collect_vars, 
                       command = 'collect',
                       logger_name = 'test')
        message = '---PARAMETERS---\n'\
                      'Command: collect\n'\
                          'binary_path: bp\n'\
                              'strict_span: sp\n'\
                                  'neighbourhood_size: n\n'\
                                      'write_genomes: wg\n'\
                                          'email: e\n'\
                                              'filenames: f\n'\
                                                  'results_dir: rd\n\n\n'
        
        assert caplog.messages == [message]
    
    def test_sieve(self, initialised_log : logging.Logger, sieve_vars : dict, 
                   caplog):
        with caplog.at_level(logging.INFO):
            log_params(sieve_vars, 
                       command = 'sieve',
                       logger_name = 'test')
        message = '---PARAMETERS---\n'\
                    'Command: sieve\n'\
                        'input_genbank_dir: igd\n'\
                            'e_value: ev\n'\
                                'min_percent_identity: mpi\n'\
                                    'similarity_filter: sf\n'\
                                        'results_dir: rd\n'\
                                            'output_blast_dir: obd\n'\
                                                'output_genbank_dir: ogd\n'\
                                                    'output_vis_dir: ovd\n'\
                                                        'min_edge_view: mev\n\n\n'
        assert caplog.messages == [message]
        
        
    
    def test_bad_collect(self, initialised_log : logging.Logger, bad_collect_vars : dict):
        with pytest.raises(KeyError):
            log_params(bad_collect_vars, 
                       command = 'collect',
                       logger_name = 'test')
    
    def test_bad_sieve(self, initialised_log : logging.Logger, bad_sieve_vars : dict):
        with pytest.raises(KeyError):
            log_params(bad_sieve_vars, 
                       command = 'sieve',
                       logger_name = 'test')
    def test_reverse(self, initialised_log : logging.Logger, collect_vars : dict):
        with pytest.raises(KeyError):
            log_params(collect_vars, 
                       command = 'sieve',
                       logger_name = 'test')
        
    def test_wrong_command(self, initialised_log : logging.Logger, collect_vars : dict):
        with pytest.raises(KeyError):
            log_params(collect_vars, 
                       command = 'not_recognised',
                       logger_name = 'test')
        
        

         

def test_initalise_log(log_setup):
    #TODO should check for stream handler but not sure which attribute it is 
    #stored in 
    logger = initalise_log('path', 
                            'name')
    assert logger.name == 'name'
    
class TestIngestBinaryCsv:
    @pytest.fixture 
    def setup_good_df(self, monkeypatch):
        def mock_good_df(path : str, sep : str) -> pd.DataFrame:
            data = [['orgaism', 'accession', 1, 1, 1, 1]]
            return pd.DataFrame(data, 
                                columns = ['col']*6
                                )
        monkeypatch.setattr('pandas.read_csv', 
                            mock_good_df)

    @pytest.fixture 
    def setup_bad_df(self, monkeypatch):
        def mock_bad_df(path : str, sep : str) -> pd.DataFrame:
            data = [['orgaism', 'accession', 1, 1, 1, ]]
            return pd.DataFrame(data, 
                                columns = ['col']*5
                                )
        monkeypatch.setattr('pandas.read_csv', 
                                mock_bad_df)
        
    def test_good_csv(self, setup_good_df):
        df_out = ingest_binary_csv('binary_path', 
                                 'logger_name')
        expected_out = [['orgaism', 'accession', 1, 1, 1, 1]]
        assert df_out.values.tolist() == expected_out
        
    def test_bad_csv(self, initialised_log : logging.Logger, setup_bad_df, caplog):
        with caplog.at_level(logging.INFO):
            with pytest.raises(ValueError):
                ingest_binary_csv('binary_path', 'logger_name')
        message = "Unexpected binary file format\nexpected columns - "\
                    "['Organism', 'Scaffold', 'Start', 'End', 'Score', one or "\
                         "more Query Gene Names].\n"\
                             "Actual columns: ['col', 'col', 'col', 'col', 'col']"
        assert caplog.messages == [message]