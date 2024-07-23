# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 13:42:51 2024

@author: u03132tk
"""
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#import sys
#add D to path.  Note, make sure you try running this when the package is 
#published - need to check innit files are all required
#sys.path.insert(0, 'D:\\')
from SyntenyQC.neighbourhood import Neighbourhood, get_record, make_filepath, write_results
from urllib.request import HTTPError
from http.client import IncompleteRead
import pytest
#import mock
import os 
from general_mocks import mock_listdir
import string
 
#write_genbank_file - skip, fileio is hard to test (look into mocking at a later date, but still difficult to do more than shallow tests of what params are called)
#write_results - skip

#TODO write integration test for Neighbourhood - test error handling is as expected outside of static functions
    
@pytest.fixture(scope = 'module')
def a_neighbourhood():
    '''
    A test Neighbourhood
    '''
    return Neighbourhood(accession = 'NZ_CP042324', 
                              motif_start = 20000, 
                              motif_stop = 30000, 
                              neighbourhood_size = 20000, 
                              strict_span = True)

@pytest.mark.parametrize('scale',  ['genome', 'neighbourhood'])
def test_get_record_genome(scale : str, a_neighbourhood : Neighbourhood):
    '''
    Test that get_record(genome) returns a record of the correct scale (either 
    a genome or a specific neighbourhood within the genome).

    Parameters
    ----------
    scale : str
        Scale of record to extract ("neighbourhood" or "genome")
    a_neighbourhood : Neighbourhood
        A test Neighbourhood.
    '''
    assert scale in ['genome', 'neighbourhood']
    record = get_record (a_neighbourhood, scale)
    if scale == 'genome':
        assert record.seq == a_neighbourhood.genome.seq
    else:
        sub_start = a_neighbourhood.neighbourhood_start
        sub_end = a_neighbourhood.neighbourhood_stop
        sub_record = a_neighbourhood.genome[sub_start : sub_end]
        assert len(sub_record.seq) == sub_end - sub_start
        assert record.seq == sub_record.seq
    



def test_get_record_exceptions(a_neighbourhood):
    '''
    Test that get_record(genome) raises a ValueError if scale is not "neighbourhood" or "genome".

    Parameters
    ----------
    scale : str
        Scale of record to extract ("neighbourhood" or "genome")
    a_neighbourhood : Neighbourhood
        A test Neighbourhood.
    '''
    with pytest.raises(ValueError):
        get_record (a_neighbourhood, 'wrong_scale')

@pytest.mark.parametrize('name_type,type_map,folder,expected_outpath',  [('accession',
                                                                          {'accession' : 'accession', 
                                                                           'organism' : 'organism'
                                                                           },
                                                                          'a_folder',
                                                                          'a_folder\\accession_(1).gbk'
                                                                          ),
                                                                        ('organism', 
                                                                         {'accession' : 'accession', 
                                                                          'organism' : 'organism'
                                                                          },
                                                                         'a_folder',
                                                                         'a_folder\\organism_(1).gbk'
                                                                         ),
                                                                        ('accession',
                                                                         {'accession' : 'accession/accession', 
                                                                          'organism' : 'organism\\organism'
                                                                          },
                                                                         'a_folder',
                                                                         'a_folder\\accession_accession.gbk'
                                                                         ),
                                                                        ('organism',
                                                                         {'accession' : 'accession/accession', 
                                                                          'organism' : 'organism\\organism'
                                                                          },
                                                                         'a_folder',
                                                                         'a_folder\\organism_organism.gbk'
                                                                         )]
                         )
def test_make_filepath(name_type : str, type_map : dict, folder : str, 
                       expected_outpath : str, monkeypatch):
    '''
    

    Parameters
    ----------
    name_type : str
        Type map key.
    type_map : dict
        Value associated with name_type that should be used to name the file 
        reached via the output of make_filepath().
    folder : str
        Folder in which the file reached via the output of make_filepath() should 
        be written.
    expected_outpath : str
        Expeted make_filepath() output.
    monkeypatch
        Pytest fixture used to replace os.listdir in thr scope of this test (no 
        type hinting as MonkeyPatch is not exposed by pytest - 
        https://github.com/pytest-dev/pytest/issues/2712).

    '''
    monkeypatch.setattr(os, "listdir", mock_listdir)
    assert make_filepath(name_type, type_map, folder) == expected_outpath
    
@pytest.mark.parametrize('filenames,scale,expected_results_folder,expected_filepath',  [('accession',
                                                                                         'genome',
                                                                                          'a_folder\\genome',
                                                                                          'a_folder\\genome\\NZ_CP042324.1.gbk'
                                                                                          ),
                                                                                        
                                                                                        ('organism',
                                                                                        'genome',
                                                                                         'a_folder\\genome',
                                                                                         'a_folder\\genome\\Streptomyces coelicolor A3(2).gbk'
                                                                                         ),
                                                                                        ('accession',
                                                                                         'neighbourhood',
                                                                                          'a_folder\\neighbourhood',
                                                                                          'a_folder\\neighbourhood\\NZ_CP042324.1.gbk'
                                                                                          ),
                                                                                        
                                                                                        ('organism',
                                                                                        'neighbourhood',
                                                                                         'a_folder\\neighbourhood',
                                                                                         'a_folder\\neighbourhood\\Streptomyces coelicolor A3(2).gbk'
                                                                                         )
                                                                                        ]
                         )
def test_write_results (filenames : str, scale : str, expected_results_folder : str,
                        expected_filepath : str, a_neighbourhood : Neighbourhood, 
                        monkeypatch):
    '''
    

    Parameters
    ----------
    filenames : str
        Record data (accession or organism) that should be used to name the written neighbourhood/genome.
    scale : str
        Scale of record to be written to genbank file ("neighbourhood" or "genome".
    expected_results_folder : str
        Folder in which the neighbourhood/genome genbank file should be written.
    expected_filepath : str
        Expected path to neighbourhood/genome genbank file.
    a_neighbourhood : Neighbourhood
        A test Neighbourhood..
    monkeypatch
        Pytest fixture used to replace os.listdir in thr scope of this test (no 
        type hinting as MonkeyPatch is not exposed by pytest - 
        https://github.com/pytest-dev/pytest/issues/2712).
    '''
    
    def mock_listdir_track_input(folder : str) -> list:
        '''
        For test simplicity/predictability, replace os.listdir() with a mocked 
        function that returns the value of mock_return when called in 
        test_write_results().  Record input folder to monkey_patch_params for 
        assertion testing as per https://stackoverflow.com/a/67074212/11357695.

        Parameters
        ----------
        folder : str
            Folder supplied to os.listdir().

        Returns
        -------
        list
            List of non-existent filenames for testing purposes.
        '''
        #
        nonlocal monkey_patch_params
        monkey_patch_params['listdir'] = folder
        return ["accession.gbk", "organism.gbk"]
    
    def mock_makedirs_track_input(folder : str, exist_ok : bool) -> None:
        '''
        To avoid dir writing, replace os.makedirs() with a mocked function that 
        returns None when called in test_write_results().  Record input folder 
        to monkey_patch_params for assertion testing as per 
        https://stackoverflow.com/a/67074212/11357695.

        Parameters
        ----------
        folder : str
            Folder supplied to os.makedirs.
        exist_ok : bool
            exist_ok parameter of os.makedirs.

        Returns
        -------
        None
        '''
        nonlocal monkey_patch_params
        monkey_patch_params['makedirs'] = folder
        return None
    
    def mock_write_genbank_file_track_input(record : SeqRecord, path : str) -> None:
        '''
        To avoid file writing, replace write_genbank_file() with a mocked function 
        that returns None when called in test_write_results(). Record input folder 
        to monkey_patch_params for assertion testing as per 
        https://stackoverflow.com/a/67074212/11357695.

        Parameters
        ----------
        record : SeqRecord
            SeqRecord to write to genbank file.
        path : str
            Path of genabnk file.

        Returns
        -------
        None
        '''
        nonlocal monkey_patch_params
        monkey_patch_params['write_genbank_file'] = {'record' : record, 
                                                     'path' : path}
        return None
    
    monkey_patch_params = {}
    monkeypatch.setattr("os.listdir", mock_listdir_track_input)
    monkeypatch.setattr("os.makedirs", mock_makedirs_track_input)
    monkeypatch.setattr("SyntenyQC.neighbourhood.write_genbank_file", 
                        mock_write_genbank_file_track_input)
    write_results(results_folder = 'a_folder', 
                  neighbourhood = a_neighbourhood, 
                  filenames = filenames, 
                  scale = scale)
    if scale == 'genome':
        expected_neighbourhood = a_neighbourhood.genome
    else:
        expected_neighbourhood = a_neighbourhood.neighbourhood
    assert expected_results_folder == monkey_patch_params['makedirs']
    assert expected_filepath == monkey_patch_params['write_genbank_file']['path']
    assert expected_neighbourhood.seq == monkey_patch_params['write_genbank_file']['record'].seq
    
    

def test_init_scrape_fail(monkeypatch):
    def mock_scrape_HTTPError(self, number_of_attempts, accession):
        raise HTTPError(None, None, None, None, None)
    def mock_scrape_IncompleteRead(self, number_of_attempts, accession):
        raise IncompleteRead(None)
    def mock_scrape_ValueError(self, number_of_attempts, accession):
        raise ValueError
    monkeypatch.setattr("SyntenyQC.neighbourhood.Neighbourhood.scrape_genome", 
                        mock_scrape_HTTPError)
    result = Neighbourhood('accession', 
                            motif_start = 1, 
                            motif_stop = 2, 
                            neighbourhood_size = 3, 
                            strict_span = False)
    assert result.scrape_error == 'HTTP error'
    assert result.overlapping_termini == False
    assert result.motif_to_long == False
    
    monkeypatch.setattr("SyntenyQC.neighbourhood.Neighbourhood.scrape_genome", 
                        mock_scrape_IncompleteRead)
    result = Neighbourhood('accession', 
                            motif_start = 1, 
                            motif_stop = 2, 
                            neighbourhood_size = 3, 
                            strict_span = False)
    assert result.scrape_error == 'IncompleteRead'
    assert result.overlapping_termini == False
    assert result.motif_to_long == False
    
    monkeypatch.setattr("SyntenyQC.neighbourhood.Neighbourhood.scrape_genome", 
                        mock_scrape_ValueError)
    result = Neighbourhood('accession', 
                            motif_start = 1, 
                            motif_stop = 2, 
                            neighbourhood_size = 3, 
                            strict_span = False)
    assert result.scrape_error == 'ValueError - No genome found'
    assert result.overlapping_termini == False
    assert result.motif_to_long == False

def test_init_define_neighbourhod_fail(monkeypatch):
    def mock_scrape_ok(self, number_of_attempts, accession):#or pass legit accession
        pass
    monkeypatch.setattr("SyntenyQC.neighbourhood.Neighbourhood.scrape_genome", 
                        mock_scrape_ok)
    result = Neighbourhood('accession', 
                            motif_start = 1, 
                            motif_stop = 4, 
                            neighbourhood_size = 2, 
                            strict_span = False)
    assert result.scrape_error == None
    assert result.overlapping_termini == False
    assert result.motif_to_long == True

def test_init_sanitise_neighbourhod(monkeypatch):
    def mock_scrape_seq(self, number_of_attempts, accession):#or pass legit accession
        return SeqRecord(seq = Seq('c'*50),
                         features = [])
    monkeypatch.setattr("SyntenyQC.neighbourhood.Neighbourhood.scrape_genome", 
                        mock_scrape_seq)
    result = Neighbourhood('accession', 
                            motif_start = 1, 
                            motif_stop = 49, 
                            neighbourhood_size = 60, 
                            strict_span = False)
    assert result.scrape_error == None
    assert result.overlapping_termini == False
    assert result.motif_to_long == False
    assert result.neighbourhood_start == 0
    assert result.neighbourhood_stop == 50
    result = Neighbourhood('accession', 
                            motif_start = 1, 
                            motif_stop = 49, 
                            neighbourhood_size = 60, 
                            strict_span = True)
    assert result.scrape_error == None
    assert result.overlapping_termini == True
    assert result.motif_to_long == False


def test_init_get_record(monkeypatch):
    def mock_scrape_seq(self, number_of_attempts, accession):#or pass legit accession
        return SeqRecord(seq = Seq(string.ascii_lowercase),
                         features = [])
    monkeypatch.setattr("SyntenyQC.neighbourhood.Neighbourhood.scrape_genome", 
                        mock_scrape_seq)
    result = Neighbourhood('accession', 
                            motif_start = 10, 
                            motif_stop = 15, 
                            neighbourhood_size = 10, 
                            strict_span = False)
    assert result.scrape_error == None
    assert result.overlapping_termini == False
    assert result.motif_to_long == False
    assert result.neighbourhood_start == 8
    assert result.neighbourhood_stop == 17
    assert result.genome.seq == string.ascii_lowercase
    assert result.neighbourhood.seq == string.ascii_lowercase[8 : 17]
    result = Neighbourhood('accession', 
                            motif_start = 1, 
                            motif_stop = 49, 
                            neighbourhood_size = 60, 
                            strict_span = True)
    assert result.scrape_error == None
    assert result.overlapping_termini == True
    assert result.motif_to_long == False
    
    
    

class TestNeighbourhood:
    '''
    Tests for the collect.Neighbourhood class
    '''
    
    Entrez.email = 'tdjkirkwood@hotmail.com'
    
    @pytest.fixture(scope = 'module')
    def read_record(self):
        '''
        Read in local genbank file to a biopython SeqRecord for downstream tests.
        Downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP042324.1/ on 10/7/2024
        '''
        test_path = 'D:/SyntenyQC/tests/NZ_CP042324.gb'
        with open(test_path, 'r') as test_file:
            test_record = SeqIO.read(test_file, 'genbank')
        return test_record
    
    
        
        
    #    do this
    #    note - cant test final exception - either it will be raised or it an error wont be raised in which case value error is raised
    
    def test_scrape_genome(self, read_record : SeqRecord):
        '''
        Test that Neighbourhood.scrape_genome() downloads a record from NCBI, 
        using a manually-downloaded local genbank file accession, that 
        has the same sequence as the local genbank. Note, SeqRecord comparison 
        is not supported by Bioython, and explicit checks of genbank SeqFeature 
        equality are not performed. Check that annotations (organism etc) are 
        identical for the scraped and local records.

        Parameters
        ----------
        read_record : SeqRecord
            Biopython SeqRecord read from local genbank file (see read_record()).
        '''
        scraped_data = Neighbourhood.scrape_genome(5, 
                                                   read_record.id)
        assert read_record.seq == scraped_data.seq
        assert read_record.annotations == scraped_data.annotations
        
    def test_scrape_genome_exceptions(self):
        '''
        Test that Neighbourhood.scrape_genome() raises a HTTPError if an accession 
        is supplied that is not identified in NCBI.  Note, Neighbourhood.scrape_genome()
        can also raise an IncompleteRead exception, but this is not tested as 
        the exact circumstances that lead to it are hard to replicate. 

        '''
        with pytest.raises(HTTPError):
            Neighbourhood.scrape_genome(5, 
                                        'wrong_accession')
            
    @pytest.mark.parametrize("motif_start,motif_stop,neighbourhood_size, expected_value", 
                             [(5, 10, 10, (3, 12)), 
                              (-10, -5, 10, (-12, -3))
                              ]
                             )
    def test_define_neighbourhood(self, motif_start : int, motif_stop : int, 
                                  neighbourhood_size : int, expected_value : tuple):
        '''
        Test that Neighbourhood.define_neighbourhood() defines neighbourhood start 
        and ends from defined inputs.  No sanity cheking is perfoemd for:
            - neighbourhood size (argparse interface checks it is >0).
            - motif_start/motif_stop (done by check_motif() to integrate with logging).  
            - returned neighbourhood values (done by Neighbourhood.sanitise_neighbourhood()). 
        ValueError should be raised if motif_stop - motif_start > neighbourhood_size. 

        Parameters
        ----------
        motif_start : int
            Start location of first cblaster hit in genome A.
        motif_stop : int
            End location of last cblaster hit in genome A.
        neighbourhood_size : int
            Size of neighbourhood defined in genome A from motif_start and motif_end.
        expected_value : tuple
            Expected value returned by Neighbourhood.define_neighbourhood() for the given inputs.
        '''
        assert expected_value == Neighbourhood.define_neighbourhood(motif_start,# = 5, 
                                                  motif_stop,# = 10, 
                                                  neighbourhood_size# = 10, 
                                                  )

    
    def test_define_neighbourhood_exceptions(self):
        '''
        Test that Neighbourhood.define_neighbourhood() raises a ValueError if 
        motif_stop - motif_start > neighbourhood_size. 
        '''
        with pytest.raises(ValueError):
            Neighbourhood.define_neighbourhood(motif_start = 5, 
                                               motif_stop = 10, 
                                               neighbourhood_size = 2 
                                               )
    
    
    @pytest.mark.parametrize("neighbourhood_start,neighbourhood_stop,genome_length,strict_span,expected_values", 
                             [(-5, 5, 10, False, (0, 5)), 
                              (0, 15, 10, False, (0, 10)),
                              (3, 4, 10, True, (3, 4)),
                              (3, 4, 10, False, (3, 4))]
                             ) 
    def test_sanitise_neighbourhood(self, neighbourhood_start : int, 
                                    neighbourhood_stop : int, genome_length : int,
                                    strict_span : bool, expected_values : tuple):
        '''
        Test that Neighbourhood.sanitise_neighbourhood() correctly sanitises 
        neighbourhood start/stop locii returned by Neighbourhood.define_neighbourhood().
        If neighbourhood_start >= 0 and neighbourhood_stop <= genome_length, they 
        should be returned with no changes. If strict_span is False and 
        neighbourhood_start < 0, it should be set to 0.  If strict_span is False 
        and neighbourhood_stop is > genome_length, it should be set to genome_length. 

        Parameters
        ----------
        neighbourhood_start : int
            Start of neighbourhood defined by Neighbourhood.define_neighbourhood().
        neighbourhood_stop : int
            End of neighbourhood defined by Neighbourhood.define_neighbourhood().
        genome_length : int
            Length of record in which the neighbourhood is being defined 
            (i.e. maximum value for neighbourhood_stop).
        strict_span : bool
            Specify if the neighbourhood should be discarded if neighbourhood_start < 0 
            or neighbourhood_stop > genome_length
        expected_values : tuple
            Expected value returned by Neighbourhood.sanitise_neighbourhood() for the given inputs..

        Returns
        -------
        None.

        '''
        
        assert expected_values == Neighbourhood.sanitise_neighbourhood(neighbourhood_start, 
                                                                    neighbourhood_stop, 
                                                                    genome_length, 
                                                                    strict_span)
         

    @pytest.mark.parametrize("neighbourhood_start,neighbourhood_stop", 
                             [(-5, 5), 
                              (0, 15),
                              ]
                             ) 
    def test_sanitise_neighbourhood_exceptions(self, neighbourhood_start : int, 
                                               neighbourhood_stop : int):
        '''
        Test that Neighbourhood.sanitise_neighbourhood() raises a ValueError if 
        neighbourhood_start < 0 or neighbourhood_stop > genome_length and strict_span
        is True.

        Parameters
        ----------
        neighbourhood_start : int
            Start of neighbourhood defined by Neighbourhood.define_neighbourhood().
        neighbourhood_stop : int
            End of neighbourhood defined by Neighbourhood.define_neighbourhood().
        '''
        with pytest.raises(ValueError):
            Neighbourhood.sanitise_neighbourhood(neighbourhood_start, 
                                                 neighbourhood_stop, 
                                                 genome_length = 10, 
                                                 strict_span = True)
            
    @pytest.mark.parametrize('neighbourhood_start,neighbourhood_stop,first_cds,last_cds',
                             [(1843790, 1860786, 'FQ762_RS08775', 'FQ762_RS08855'),
                              (1843792, 1860784, 'FQ762_RS08780', 'FQ762_RS08850')
                              ]
                             )
    def test_get_neighbourhood(self, neighbourhood_start : int, 
                               neighbourhood_stop : int, first_cds : str, 
                               last_cds : str, read_record : SeqRecord):
        '''
        Test that Neighbourhood.get_neighbourhood() successfully slices original 
        SeqRecord to get a record for the desired Neighbourhood.  No sanity checking 
        is performed for neighbourhood_start or neighbourhood_stop (this is done 
        by Neighbourhood.sanitise_neighbourhood()).  Check that if a feature is 
        not completely within the neighbourhood, it is excluded (feature locci 
        defined following manual examination of the genbank file in SnapGene).  
        Check that annotations (organism etc) are copied to the neighbourhood record.

        Parameters
        ----------
        neighbourhood_start : int
            Start of neighbourhood defined by Neighbourhood.define_neighbourhood().
        neighbourhood_stop : int
            End of neighbourhood defined by Neighbourhood.define_neighbourhood().
        first_cds : str
            first CDS locus_tag expected in the record after neighbourhood_start.
        last_cds : str
            last CDS locus_tag expected in the record after neighbourhood_stop.
        read_record : SeqRecord
            Biopython SeqRecord read from local genbank file (see read_record()).
        '''
        cut_neighbourhood = Neighbourhood.get_neighbourhood(neighbourhood_start, 
                                                            neighbourhood_stop, 
                                                            record = read_record)
        cut_cds = [f for f in cut_neighbourhood.features if f.type == 'CDS']
        assert ''.join(cut_cds[0].qualifiers['locus_tag']) == first_cds, cut_cds[0].qualifiers['locus_tag']
        assert ''.join(cut_cds[-1].qualifiers['locus_tag']) == last_cds, cut_cds[-1].qualifiers['locus_tag']
        assert cut_neighbourhood.annotations == read_record.annotations
        assert cut_neighbourhood.seq == read_record.seq[neighbourhood_start : neighbourhood_stop]



