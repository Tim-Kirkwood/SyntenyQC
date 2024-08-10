# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 13:42:51 2024

@author: u03132tk
"""
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from SyntenyQC.neighbourhood import Neighbourhood, get_record, make_filepath, write_results, log_neighbourhood_details
from urllib.request import HTTPError
from http.client import IncompleteRead
import pytest
import logging
import os 
from general_mocks import mock_listdir
import string


# =============================================================================
# Fixtures
# =============================================================================
    
@pytest.fixture 
def neighbourhood_http_error(monkeypatch) -> Neighbourhood:
    @staticmethod
    def mock_scrape_HTTPError(number_of_attempts : int, accession : str):
        raise HTTPError(None, None, None, None, None)
    monkeypatch.setattr("SyntenyQC.neighbourhood.Neighbourhood.scrape_genome", 
                        mock_scrape_HTTPError)
    return Neighbourhood('accession', 
                         motif_start = 1, 
                         motif_stop = 2, 
                         neighbourhood_size = 3, 
                         strict_span = False)

@pytest.fixture 
def neighbourhood_incomplete_read(monkeypatch) -> Neighbourhood:
    @staticmethod
    def mock_scrape_IncompleteRead(number_of_attempts : int, accession : str):
        raise IncompleteRead(None)
    monkeypatch.setattr("SyntenyQC.neighbourhood.Neighbourhood.scrape_genome", 
                        mock_scrape_IncompleteRead)
    return Neighbourhood('accession', 
                         motif_start = 1, 
                         motif_stop = 2, 
                         neighbourhood_size = 3, 
                         strict_span = False)

@pytest.fixture 
def neighbourhood_value_error(monkeypatch) -> Neighbourhood:
    @staticmethod
    def mock_scrape_ValueError(number_of_attempts : int, accession : str):
        raise ValueError
    monkeypatch.setattr("SyntenyQC.neighbourhood.Neighbourhood.scrape_genome", 
                        mock_scrape_ValueError)
    return Neighbourhood('accession', 
                            motif_start = 1, 
                            motif_stop = 2, 
                            neighbourhood_size = 3, 
                            strict_span = False)

@pytest.fixture 
def setup_scrape_alphabet(monkeypatch):
    @staticmethod
    def mock_scrape_seq(number_of_attempts : int, accession : str) -> SeqRecord:
        return SeqRecord(seq = Seq(string.ascii_lowercase),
                         features = [],
                         id = 'pseudo_accession',
                         annotations = {'organism' : 'pseudo_organism/with\\slashes'}
                         )
    monkeypatch.setattr("SyntenyQC.neighbourhood.Neighbourhood.scrape_genome", 
                        mock_scrape_seq)

@pytest.fixture 
def neighbourhood_motif_to_long(setup_scrape_alphabet) -> Neighbourhood:
    return Neighbourhood('pseudo_accession', 
                         motif_start = 1, 
                         motif_stop = 60, 
                         neighbourhood_size = 25, 
                         strict_span = False)

@pytest.fixture 
def neighbourhood_overlapping_termini(setup_scrape_alphabet) -> Neighbourhood:
    return Neighbourhood('pseudo_accession', 
                            motif_start = 1, 
                            motif_stop = 25, 
                            neighbourhood_size = 60, 
                            strict_span = True)

@pytest.fixture 
def pseudo_neighbourhood_sanitised(setup_scrape_alphabet) -> Neighbourhood:
    return Neighbourhood('pseudo_accession', 
                         motif_start = 1, 
                         motif_stop = 25, 
                         neighbourhood_size = 60, 
                         strict_span = False)


@pytest.fixture 
def pseudo_neighbourhood(setup_scrape_alphabet) -> Neighbourhood:
    return Neighbourhood('pseudo_accession', 
                        motif_start = 10, 
                        motif_stop = 15, 
                        neighbourhood_size = 10, 
                        strict_span = True)

@pytest.fixture(scope = 'module')
def local_record() -> SeqRecord:
    test_path = f'{os.path.dirname(os.path.abspath(__file__))}\\NZ_CP042324.gb'
    with open(test_path, 'r') as test_file:
        test_record = SeqIO.read(test_file, 'genbank')
    return test_record






# =============================================================================
# Unit test Neighbourhood methods
# =============================================================================

class TestScrapeGenome:
    
    @pytest.fixture
    def email_address(self) -> str:
        path = f'{os.path.dirname(os.path.abspath(__file__))}\\email.txt'
        with open(path, 'r') as file:
            email = file.read() 
        return email    
    
    def test_normal(self,
                    #function args
                    local_record : SeqRecord, email_address : str,
                    #print email - only need to print once, so only include in this function
                    capsys):
        with capsys.disabled():
            print(f"EMAIL: {email_address}")
        updated_email = email_address != 'an_email@domain.com' 
        ok_email = '@' in email_address
        if ok_email and updated_email:
            Entrez.email = email_address
            scraped_data = Neighbourhood.scrape_genome(5, 
                                                       local_record.id)
            assert local_record.seq == scraped_data.seq
            assert local_record.annotations == scraped_data.annotations
        else:
            if not updated_email:
                assert False, f'please update email from an_email@domain.com in {os.path.dirname(os.path.abspath(__file__))}\\email.txt'
            else:
                assert False, f'{email_address} must have @'
        
    def test_exceptions(self, email_address : str):
        updated_email = email_address != 'an_email@domain.com' 
        ok_email = '@' in email_address
        if ok_email and updated_email:
            Entrez.email = email_address
            with pytest.raises(HTTPError):
                Neighbourhood.scrape_genome(5, 
                                            'wrong_accession')
            with pytest.raises(ValueError):
                Neighbourhood.scrape_genome(5, 
                                            '')
            #TODO work out how to replicate IncompleteRead - mocking seqio.read() to 
            #raise an error and then confirming the error was raised is a bit pointless
        else:
            if not updated_email:
                assert False, f'please update email from an_email@domain.com in {os.path.dirname(os.path.abspath(__file__))}\\email.txt'
            else:
                assert False, f'{email_address} must have @'

class TestDefineNeighbourhood:
    
    @pytest.mark.parametrize("motif_start,motif_stop,neighbourhood_size, expected_value", 
                             [(5, 10, 10, (3, 12)), 
                              (-10, -5, 10, (-12, -3))]
                             )
    def test_normal(self, motif_start : int, motif_stop : int, 
                    neighbourhood_size : int, expected_value : tuple):
        assert expected_value == Neighbourhood.define_neighbourhood(motif_start, 
                                                                    motif_stop,
                                                                    neighbourhood_size 
                                                                    )
    
    def test_exceptions(self):
        with pytest.raises(ValueError):
            Neighbourhood.define_neighbourhood(motif_start = 5, 
                                               motif_stop = 10, 
                                               neighbourhood_size = 2 
                                               )
    
class TestSanitiseNeighbourhood:
    
    @pytest.mark.parametrize("neighbourhood_start,neighbourhood_stop,genome_length,strict_span,expected_values", 
                             [(-5, 5, 10, False, (0, 5)), 
                              (0, 15, 10, False, (0, 10)),
                              (3, 4, 10, True, (3, 4)),
                              (3, 4, 10, False, (3, 4))]
                             ) 
    def test_normal(self, neighbourhood_start : int, neighbourhood_stop : int, 
                    genome_length : int, strict_span : bool, expected_values : tuple):
        assert expected_values == Neighbourhood.sanitise_neighbourhood(neighbourhood_start, 
                                                                    neighbourhood_stop, 
                                                                    genome_length, 
                                                                    strict_span)

    @pytest.mark.parametrize("neighbourhood_start,neighbourhood_stop", 
                             [(-5, 5), 
                              (0, 15),
                              ]
                             ) 
    def test_exceptions(self, neighbourhood_start : int, neighbourhood_stop : int):
        with pytest.raises(ValueError):
            Neighbourhood.sanitise_neighbourhood(neighbourhood_start, 
                                                 neighbourhood_stop, 
                                                 genome_length = 10, 
                                                 strict_span = True)

@pytest.mark.parametrize('neighbourhood_start,neighbourhood_stop,first_cds,last_cds',
                         #NOTE - first and last cds were determined based on 
                         #manual inspection of record in snapgene
                         [(1843790, 1860786, 'FQ762_RS08775', 'FQ762_RS08855'),
                          (1843792, 1860784, 'FQ762_RS08780', 'FQ762_RS08850')
                          ]
                         )
def test_get_neighbourhood(neighbourhood_start : int, neighbourhood_stop : int, 
                           first_cds : str, last_cds : str, local_record : SeqRecord):
    cut_neighbourhood = Neighbourhood.get_neighbourhood(neighbourhood_start, 
                                                        neighbourhood_stop, 
                                                        record = local_record)
    cut_cds = [f for f in cut_neighbourhood.features if f.type == 'CDS']
    assert ''.join(cut_cds[0].qualifiers['locus_tag']) == first_cds
    assert ''.join(cut_cds[-1].qualifiers['locus_tag']) == last_cds
    assert cut_neighbourhood.annotations == local_record.annotations
    assert cut_neighbourhood.seq == local_record.seq[neighbourhood_start : neighbourhood_stop]
        
    
    
        
# =============================================================================
# Neighbourhood integration tests
# =============================================================================

class TestNeighbourhood:
    
    def test_scrape_http_error(self, neighbourhood_http_error : Neighbourhood):
        assert neighbourhood_http_error.scrape_error == 'HTTP error'
        assert neighbourhood_http_error.overlapping_termini == False
        assert neighbourhood_http_error.motif_to_long == False
    
    def test_scrape_incomplete_read(self, neighbourhood_incomplete_read : Neighbourhood):
        assert neighbourhood_incomplete_read.scrape_error == 'IncompleteRead'
        assert neighbourhood_incomplete_read.overlapping_termini == False
        assert neighbourhood_incomplete_read.motif_to_long == False
    
    def test_scrape_value_error(self, neighbourhood_value_error : Neighbourhood):
        assert neighbourhood_value_error.scrape_error == 'ValueError'
        assert neighbourhood_value_error.overlapping_termini == False
        assert neighbourhood_value_error.motif_to_long == False
    
    def test_define_neighbourhod_motif_to_long(self, neighbourhood_motif_to_long : Neighbourhood):
        assert neighbourhood_motif_to_long.scrape_error == None
        assert neighbourhood_motif_to_long.overlapping_termini == False
        assert neighbourhood_motif_to_long.motif_to_long == True

    def test_sanitise_neighbourhod_fail(self, neighbourhood_overlapping_termini : Neighbourhood):
        assert neighbourhood_overlapping_termini.scrape_error == None
        assert neighbourhood_overlapping_termini.overlapping_termini == True
        assert neighbourhood_overlapping_termini.motif_to_long == False
    
    def test_sanitise_neighbourhod(self, pseudo_neighbourhood_sanitised : Neighbourhood):
        assert pseudo_neighbourhood_sanitised.scrape_error == None
        assert pseudo_neighbourhood_sanitised.overlapping_termini == False
        assert pseudo_neighbourhood_sanitised.motif_to_long == False
        assert pseudo_neighbourhood_sanitised.neighbourhood_start == 0
        assert pseudo_neighbourhood_sanitised.neighbourhood_stop == 26
    
    def test_get_neighbourhood(self, pseudo_neighbourhood : Neighbourhood):
        assert pseudo_neighbourhood.scrape_error == None
        assert pseudo_neighbourhood.overlapping_termini == False
        assert pseudo_neighbourhood.motif_to_long == False
        assert pseudo_neighbourhood.neighbourhood_start == 8
        assert pseudo_neighbourhood.neighbourhood_stop == 17
        assert pseudo_neighbourhood.genome.seq == string.ascii_lowercase
        assert pseudo_neighbourhood.neighbourhood.seq == string.ascii_lowercase[8 : 17]




# =============================================================================
# Function tests
# =============================================================================

class TestGetRecord:
    
    def test_genome(self, pseudo_neighbourhood : Neighbourhood):
        record = get_record (pseudo_neighbourhood, 
                             scale = 'genome')
        assert record.seq == pseudo_neighbourhood.genome.seq
        
    def test_neighbourhood(self, pseudo_neighbourhood : Neighbourhood):
        record = get_record (pseudo_neighbourhood, 
                             scale = 'neighbourhood')
        sub_start = pseudo_neighbourhood.neighbourhood_start
        sub_end = pseudo_neighbourhood.neighbourhood_stop
        sub_record = pseudo_neighbourhood.genome[sub_start : sub_end]
        assert record.seq == sub_record.seq
        
    def test_get_record_exceptions(self, pseudo_neighbourhood : Neighbourhood):
        with pytest.raises(ValueError):
            get_record (pseudo_neighbourhood, 
                        'wrong_scale')

@pytest.mark.parametrize('name_type,type_map,folder,expected_outpath',  
                         [('accession', {'accession' : 'accession', 
                                         'organism' : 'organism'},
                           'a_folder','a_folder\\accession_(1).gbk'
                           ),
                          ('organism', {'accession' : 'accession', 
                                        'organism' : 'organism'},
                           'a_folder', 'a_folder\\organism_(1).gbk'
                           ),
                          #some organism names can have slashes
                          ('accession', {'accession' : 'accession/accession', 
                                         'organism' : 'organism\\organism'},
                           'a_folder', 'a_folder\\accession_accession.gbk'
                           ),
                          ('organism', {'accession' : 'accession/accession', 
                                        'organism' : 'organism\\organism'},
                           'a_folder', 'a_folder\\organism_organism.gbk'
                           )
                          ]
                         )
def test_make_filepath(name_type : str, type_map : dict, folder : str, 
                       expected_outpath : str, monkeypatch):
    monkeypatch.setattr(os, "listdir", mock_listdir)
    assert make_filepath(name_type, type_map, folder) == expected_outpath
    
class TestWriteResults:
        
    @pytest.fixture
    def setup_and_track(self, monkeypatch) -> dict:
        
        def mock_listdir_track_input(folder : str) -> list:
            nonlocal monkey_patch_params
            monkey_patch_params['listdir'] = folder
            return ["accession.gbk", "organism.gbk"]
        
        def mock_makedirs_track_input(folder : str, exist_ok : bool):
            nonlocal monkey_patch_params
            monkey_patch_params['makedirs'] = folder
        
        def mock_write_genbank_file_track_input(record : SeqRecord, path : str):
            nonlocal monkey_patch_params
            monkey_patch_params['write_genbank_file'] = {'record' : record, 
                                                         'path' : path}
            
        monkey_patch_params = {}
        monkeypatch.setattr("os.listdir", 
                            mock_listdir_track_input)
        monkeypatch.setattr("os.makedirs", 
                            mock_makedirs_track_input)
        monkeypatch.setattr("SyntenyQC.neighbourhood.write_genbank_file", 
                            mock_write_genbank_file_track_input)
        return monkey_patch_params
        
    @pytest.mark.parametrize('filenames,expected_results_folder,expected_filepath',  
                             [('accession','a_folder\\genome', 
                               'a_folder\\genome\\pseudo_accession.gbk'),
                              ('organism', 'a_folder\\genome',
                               'a_folder\\genome\\pseudo_organism_with_slashes.gbk'),
                              ]
                             )
    def test_genome(self,
                    #function args
                    filenames : str, expected_results_folder : str,
                    expected_filepath : str, pseudo_neighbourhood : Neighbourhood, 
                    #test setup
                    log_setup, setup_and_track, caplog):
        with caplog.at_level(logging.INFO):
            results = write_results(results_folder = 'a_folder', 
                                    neighbourhood = pseudo_neighbourhood, 
                                    filenames = filenames, 
                                    scale = 'genome',
                                    logger_name = 'collect')
        assert expected_results_folder == setup_and_track['makedirs']
        assert expected_filepath == setup_and_track['write_genbank_file']['path']
        assert pseudo_neighbourhood.genome.seq == setup_and_track['write_genbank_file']['record'].seq
        assert caplog.messages == [f'written {pseudo_neighbourhood.accession} GENOME to {results}']
        caplog.clear()

    @pytest.mark.parametrize('filenames,expected_results_folder,expected_filepath',  
                             [('accession', 'a_folder\\neighbourhood',
                               'a_folder\\neighbourhood\\pseudo_accession.gbk'),
                              ('organism', 'a_folder\\neighbourhood',
                               'a_folder\\neighbourhood\\pseudo_organism_with_slashes.gbk')]
                             )
    def test_neighbourhood(self, 
                           #function args
                           filenames : str, expected_results_folder : str, 
                           expected_filepath : str, pseudo_neighbourhood : Neighbourhood,
                           #test setup
                           monkeypatch, log_setup, setup_and_track, caplog):
        with caplog.at_level(logging.INFO):
            results = write_results(results_folder = 'a_folder', 
                                    neighbourhood = pseudo_neighbourhood, 
                                    filenames = filenames, 
                                    scale = 'neighbourhood',
                                    logger_name = 'collect')
        assert expected_results_folder == setup_and_track['makedirs']
        assert expected_filepath == setup_and_track['write_genbank_file']['path']
        assert pseudo_neighbourhood.neighbourhood.seq == setup_and_track['write_genbank_file']['record'].seq
        assert caplog.messages == [f'written {pseudo_neighbourhood.accession} NEIGHBOURHOOD to {results}']
        caplog.clear()

class TestLogNeighbourhoodDetails:
    
    def test_log_http_error(self,
                            #function arg
                            neighbourhood_http_error : Neighbourhood, 
                            #test setup
                            caplog, log_setup):
        with caplog.at_level(logging.INFO):
            log_neighbourhood_details(neighbourhood_http_error, 
                                      logger_name = 'collect')
        assert caplog.messages == ['scrape_fail - HTTP error - accession']
        
    def test_log_incomplete_read(self, 
                                 #function arg
                                 neighbourhood_incomplete_read : Neighbourhood, 
                                 #test setup
                                 caplog, log_setup):
        with caplog.at_level(logging.INFO):
            log_neighbourhood_details(neighbourhood_incomplete_read, 
                                      logger_name = 'collect')
        assert caplog.messages == ['scrape_fail - IncompleteRead - accession']
    
    def test_log_value_error(self, 
                             #function arg
                             neighbourhood_value_error : Neighbourhood, 
                             #test setup
                             caplog, log_setup):
        with caplog.at_level(logging.INFO):
            log_neighbourhood_details(neighbourhood_value_error, 
                                      logger_name = 'collect')
        assert caplog.messages == ['scrape_fail - ValueError - accession']
    
    def test_log_motif_to_long(self, 
                               #function arg
                               neighbourhood_motif_to_long : Neighbourhood, 
                               #test setup
                               caplog, log_setup):
        with caplog.at_level(logging.INFO):
            log_neighbourhood_details(neighbourhood_motif_to_long, 
                                      logger_name = 'collect')
        assert caplog.messages == ['motif is longer than specified neighbourhood - pseudo_accession']
        
    def test_log_overlapping_termini(self, 
                                     #function arg
                                     neighbourhood_overlapping_termini : Neighbourhood, 
                                     #test setup
                                     caplog, log_setup):
        with caplog.at_level(logging.INFO):
            log_neighbourhood_details(neighbourhood_overlapping_termini, 
                                      logger_name = 'collect')
        assert caplog.messages == ['overlapping_termini - pseudo_accession']
    
    def test_log_neighbourhood(self, 
                               #function arg
                               pseudo_neighbourhood : Neighbourhood, 
                               #test setup
                               caplog, log_setup):
        with caplog.at_level(logging.INFO):
            log_neighbourhood_details(pseudo_neighbourhood, 
                                      logger_name = 'collect')
        assert caplog.messages == ['neighbourhood = 8 -> 17']