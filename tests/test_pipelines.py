# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 17:01:26 2024

@author: u03132tk
"""
import logging
import pandas as pd
import pytest 
from SyntenyQC.pipelines import collect, sieve
from SyntenyQC.neighbourhood import Neighbourhood
from general_mocks import mock_get_gbk_files, mock_makedirs
import networkx as nx





# =============================================================================
# TestClasses
# =============================================================================
class TestCollect:
    
    @pytest.fixture
    def motif_error_params(self) -> str:
        #same params are used when checking for motif errors in input data 
        return '---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                   'strict_span: False\nneighbourhood_size: 200\nwrite_genomes: False\n'\
                       'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                           'results_dir: a_dir\n\n\n'
                           
    @pytest.fixture 
    def collect_setup(self, monkeypatch):
        def mock_good_df(path : str, sep : str) -> pd.DataFrame:
            data = [
                ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', 100, 200, 15.54, 1],#good
                ['Streptomyces anthocyanicus IPS92w', 'bad_accession', 100, 200, 15.54, 1],#scrape error
                ['Streptomyces anthocyanicus IPS92w', '', 100, 200, 15.54, 1],#empty accession
                ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', 1, 100, 15.54, 1],#overlapping termini
                ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', 1, 10**6, 15.54, 1]#motif too large - note will be extended from MIDPOINT of start/end - i.e. 5**6 
                ]
            return pd.DataFrame(data, 
                                columns = ['Organism', 'Scaffold', 'Start', 'End', 
                                           'Score', 'a_locus_tag']
                                )
                                
        def mock_write_results(results_folder : str, neighbourhood : Neighbourhood,
                               filenames : str, scale : str, logger_name : str) -> str:
            path = f'path_to_{scale}.gbk'
            logger = logging.getLogger(logger_name)
            logger.info(f'written {neighbourhood.accession} {scale.upper()} to {path}')
            return path
        
        monkeypatch.setattr('pandas.read_csv', 
                            mock_good_df)
        monkeypatch.setattr('SyntenyQC.pipelines.write_results', 
                            mock_write_results)
    
    def test_sp_wg(self, monkeypatch, caplog, log_setup, collect_setup):
        with caplog.at_level(logging.INFO):
            output = collect(binary_path = 'a_binary_path', 
                             strict_span=True, 
                             neighbourhood_size=200, 
                             write_genomes = True, 
                             email = 'tdjkirkwood@hotmail.com', 
                             filenames = 'doesnt matter', 
                             results_dir = 'a_dir')
        assert output == 'a_dir'
        message = ['---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                       'strict_span: True\nneighbourhood_size: 200\nwrite_genomes: True\n'\
                           'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                               'results_dir: a_dir\n\n\n', 
                   'Extracting 5 gbks...', 
                   'accession NZ_JAPFQR010000057.1 #0 of 4', 
                   'motif = 100 -> 200', 
                   'neighbourhood = 50 -> 250', 
                   'written NZ_JAPFQR010000057.1 NEIGHBOURHOOD to path_to_neighbourhood.gbk', 
                   'written NZ_JAPFQR010000057.1 GENOME to path_to_genome.gbk', 
                   'accession bad_accession #1 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - HTTP error - bad_accession', 
                   'accession  #2 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - ValueError - ', 
                   'accession NZ_JAPFQR010000057.1 #3 of 4', 
                   'motif = 1 -> 100', 
                   'overlapping_termini - NZ_JAPFQR010000057.1', 
                   'accession NZ_JAPFQR010000057.1 #4 of 4', 
                   'motif = 1 -> 1000000', 
                   'motif is longer than specified neighbourhood - NZ_JAPFQR010000057.1', 
                   'Written 1 records to a_dir - 2 failed scraping, '\
                       '1 had termini that exceeded the genome boundaries, '\
                           '1 had a motif that was too long']
        assert caplog.messages == message
        caplog.clear()
            
    def test_wg(self, monkeypatch, caplog, log_setup, collect_setup):
        with caplog.at_level(logging.INFO):
            output = collect(binary_path = 'a_binary_path', 
                             strict_span=False, 
                             neighbourhood_size=200, 
                             write_genomes = True, 
                             email = 'tdjkirkwood@hotmail.com', 
                             filenames = 'doesnt matter', 
                             results_dir = 'a_dir')
        assert output == 'a_dir'
        message = ['---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                       'strict_span: False\nneighbourhood_size: 200\nwrite_genomes: True\n'\
                           'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                               'results_dir: a_dir\n\n\n', 
                   'Extracting 5 gbks...', 
                   'accession NZ_JAPFQR010000057.1 #0 of 4', 
                   'motif = 100 -> 200', 
                   'neighbourhood = 50 -> 250', 
                   'written NZ_JAPFQR010000057.1 NEIGHBOURHOOD to path_to_neighbourhood.gbk', 
                   'written NZ_JAPFQR010000057.1 GENOME to path_to_genome.gbk', 
                   'accession bad_accession #1 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - HTTP error - bad_accession', 
                   'accession  #2 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - ValueError - ', 
                   'accession NZ_JAPFQR010000057.1 #3 of 4', 
                   'motif = 1 -> 100', 
                   'neighbourhood = 0 -> 150', 
                   'written NZ_JAPFQR010000057.1 NEIGHBOURHOOD to path_to_neighbourhood.gbk', 
                   'written NZ_JAPFQR010000057.1 GENOME to path_to_genome.gbk',
                   'accession NZ_JAPFQR010000057.1 #4 of 4', 
                   'motif = 1 -> 1000000', 
                   'motif is longer than specified neighbourhood - NZ_JAPFQR010000057.1', 
                   'Written 2 records to a_dir - 2 failed scraping, '\
                       '0 had termini that exceeded the genome boundaries, '\
                           '1 had a motif that was too long']
            
        assert caplog.messages == message
        caplog.clear()
            
    def test_sp(self, monkeypatch, caplog, log_setup, collect_setup):
        with caplog.at_level(logging.INFO):
            output = collect(binary_path = 'a_binary_path', 
                             strict_span=True, 
                             neighbourhood_size=200, 
                             write_genomes = False, 
                             email = 'tdjkirkwood@hotmail.com', 
                             filenames = 'doesnt matter', 
                             results_dir = 'a_dir')
        assert output == 'a_dir'
        message = ['---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                       'strict_span: True\nneighbourhood_size: 200\nwrite_genomes: False\n'\
                           'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                               'results_dir: a_dir\n\n\n', 
                   'Extracting 5 gbks...', 
                   'accession NZ_JAPFQR010000057.1 #0 of 4', 
                   'motif = 100 -> 200', 
                   'neighbourhood = 50 -> 250', 
                   'written NZ_JAPFQR010000057.1 NEIGHBOURHOOD to path_to_neighbourhood.gbk', 
                   #'written NZ_JAPFQR010000057.1 GENOME to path_to_genome.gbk', 
                   'accession bad_accession #1 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - HTTP error - bad_accession', 
                   'accession  #2 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - ValueError - ', 
                   'accession NZ_JAPFQR010000057.1 #3 of 4', 
                   'motif = 1 -> 100', 
                   'overlapping_termini - NZ_JAPFQR010000057.1', 
                   'accession NZ_JAPFQR010000057.1 #4 of 4', 
                   'motif = 1 -> 1000000', 
                   'motif is longer than specified neighbourhood - NZ_JAPFQR010000057.1', 
                   'Written 1 records to a_dir - 2 failed scraping, '\
                       '1 had termini that exceeded the genome boundaries, '\
                           '1 had a motif that was too long']
        assert caplog.messages == message
        caplog.clear()
            
    def test_minimal(self, monkeypatch, caplog, log_setup, collect_setup):
        with caplog.at_level(logging.INFO):
            output = collect(binary_path = 'a_binary_path', 
                             strict_span=False, 
                             neighbourhood_size=200, 
                             write_genomes = False, 
                             email = 'tdjkirkwood@hotmail.com', 
                             filenames = 'doesnt matter', 
                             results_dir = 'a_dir')
        assert output == 'a_dir'
        message = ['---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                       'strict_span: False\nneighbourhood_size: 200\nwrite_genomes: False\n'\
                           'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                               'results_dir: a_dir\n\n\n', 
                   'Extracting 5 gbks...', 
                   'accession NZ_JAPFQR010000057.1 #0 of 4', 
                   'motif = 100 -> 200', 
                   'neighbourhood = 50 -> 250', 
                   'written NZ_JAPFQR010000057.1 NEIGHBOURHOOD to path_to_neighbourhood.gbk', 
                   #'written NZ_JAPFQR010000057.1 GENOME to path_to_genome.gbk', 
                   'accession bad_accession #1 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - HTTP error - bad_accession', 
                   'accession  #2 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - ValueError - ', 
                   'accession NZ_JAPFQR010000057.1 #3 of 4', 
                   'motif = 1 -> 100', 
                   'neighbourhood = 0 -> 150', 
                   'written NZ_JAPFQR010000057.1 NEIGHBOURHOOD to path_to_neighbourhood.gbk', 
                   #'written NZ_JAPFQR010000057.1 GENOME to path_to_genome.gbk',
                   'accession NZ_JAPFQR010000057.1 #4 of 4', 
                   'motif = 1 -> 1000000', 
                   'motif is longer than specified neighbourhood - NZ_JAPFQR010000057.1', 
                   'Written 2 records to a_dir - 2 failed scraping, '\
                       '0 had termini that exceeded the genome boundaries, '\
                           '1 had a motif that was too long']
        assert caplog.messages == message
        caplog.clear()
            
    def test_fail_equal_motif(self, monkeypatch, caplog, log_setup,
                              motif_error_params : str):
        #setup test
        def mock_equal_start_stop(path : str, sep : str) -> pd.DataFrame:
            data = [
                ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', 1, 1, 15.54, 1],
                ]
            return pd.DataFrame(data, 
                                columns = ['Organism', 'Scaffold', 'Start', 
                                           'End', 'Score', 'a_locus_tag']
                                )
        monkeypatch.setattr('pandas.read_csv', 
                            mock_equal_start_stop)
        #run test
        with pytest.raises(ValueError):
            with caplog.at_level(logging.INFO):
                collect(binary_path = 'a_binary_path', 
                        strict_span=False, 
                        neighbourhood_size=200, 
                        write_genomes = False,
                        email = 'tdjkirkwood@hotmail.com', 
                        filenames = 'doesnt matter', 
                        results_dir = 'a_dir')
        assert caplog.messages == [motif_error_params, 
                                   'Extracting 1 gbks...', 
                                   'accession NZ_JAPFQR010000057.1 #0 of 0',
                                   'motif = 1 -> 1',
                                   'motif_start 1 is >= motif_stop 1\n']
        caplog.clear()
        
    def test_fail_bad_start_stop(self, monkeypatch, caplog, log_setup,
                                 motif_error_params : str):
        #setup test
        def mock_bad_start_stop(path : str, sep : str) -> pd.DataFrame:
            data = [
                ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', -1, 0, 15.54, 1],
                ]
            return pd.DataFrame(data, 
                                columns = ['Organism', 'Scaffold', 'Start', 'End', 
                                           'Score', 'query_locus_tag']
                                )
        monkeypatch.setattr('pandas.read_csv', 
                            mock_bad_start_stop)
        #run test
        with pytest.raises(ValueError):
            with caplog.at_level(logging.INFO):
                collect(binary_path = 'a_binary_path', 
                        strict_span=False, 
                        neighbourhood_size=200, 
                        write_genomes = False, 
                        email = 'tdjkirkwood@hotmail.com', 
                        filenames = 'doesnt matter', 
                        results_dir = 'a_dir')
        assert caplog.messages == [motif_error_params, 
                                   'Extracting 1 gbks...', 
                                   'accession NZ_JAPFQR010000057.1 #0 of 0',
                                   'motif = -1 -> 0',
                                   'motif_start -1 < 0\nmotif_stop 0 <= 0\n']
        caplog.clear()
    
    def test_fail_bad_start(self, monkeypatch, caplog, log_setup,
                            motif_error_params : str):
        #setup test
        def mock_bad_start(path : str, sep : str) -> pd.DataFrame:
            data = [
                ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', -1, 1, 15.54, 1],
                ]
            return pd.DataFrame(data, 
                                columns = ['Organism', 'Scaffold', 'Start', 
                                           'End', 'Score', 'query_locus_tag']
                                )
        monkeypatch.setattr('pandas.read_csv', 
                            mock_bad_start)
        #run test
        with pytest.raises(ValueError):
            with caplog.at_level(logging.INFO):
                collect(binary_path = 'a_binary_path', 
                        strict_span=False, 
                        neighbourhood_size=200, 
                        write_genomes = False, 
                        email = 'tdjkirkwood@hotmail.com', 
                        filenames = 'doesnt matter', 
                        results_dir = 'a_dir')
        assert caplog.messages == [motif_error_params, 
                                   'Extracting 1 gbks...', 
                                   'accession NZ_JAPFQR010000057.1 #0 of 0',
                                   'motif = -1 -> 1',
                                   'motif_start -1 < 0\n']
        caplog.clear()
        
    def test_fail_equal_bad_start_stop(self, monkeypatch, caplog, log_setup,
                                       motif_error_params : str):
        #setup test
        def mock_equal_bad_start_stop(path : str, sep : str) -> pd.DataFrame:
            data = [
                ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', -1, -1, 15.54, 1],
                ]
            return pd.DataFrame(data, 
                                columns = ['Organism', 'Scaffold', 'Start', 'End', 
                                           'Score', 'query_locus_tag']
                                )
        monkeypatch.setattr('pandas.read_csv', 
                            mock_equal_bad_start_stop)
        #run test
        with pytest.raises(ValueError):
            with caplog.at_level(logging.INFO):
                collect(binary_path = 'a_binary_path', 
                        strict_span=False, 
                        neighbourhood_size=200, 
                        write_genomes = False, 
                        email = 'tdjkirkwood@hotmail.com', 
                        filenames = 'doesnt matter', 
                        results_dir = 'a_dir')
        assert caplog.messages == [motif_error_params, 
                                   'Extracting 1 gbks...', 
                                   'accession NZ_JAPFQR010000057.1 #0 of 0',
                                   'motif = -1 -> -1',
                                   'motif_start -1 < 0\nmotif_stop -1 <= 0\n'\
                                       'motif_start -1 is >= motif_stop -1\n'  ]
        caplog.clear()
        
    def test_fail_bad_stop(self, monkeypatch, caplog, log_setup,
                           motif_error_params : str):
        #setup test
        def mock_bad_stop(path : str, sep : str) -> pd.DataFrame:
            data = [
                ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', 2, 1, 15.54, 1],
                ]
            return pd.DataFrame(data, 
                                columns = ['Organism', 'Scaffold', 'Start', 'End', 
                                           'Score', 'query_locus_tag']
                                )
        monkeypatch.setattr('pandas.read_csv', 
                            mock_bad_stop)
        #run test
        with pytest.raises(ValueError):
            with caplog.at_level(logging.INFO):
                collect(binary_path = 'a_binary_path', 
                        strict_span=False, 
                        neighbourhood_size=200, 
                        write_genomes = False, 
                        email = 'tdjkirkwood@hotmail.com', 
                        filenames = 'doesnt matter', 
                        results_dir = 'a_dir')
        assert caplog.messages == [motif_error_params, 
                                   'Extracting 1 gbks...', 
                                   'accession NZ_JAPFQR010000057.1 #0 of 0',
                                   'motif = 2 -> 1',
                                   'motif_start 2 is >= motif_stop 1\n'  ]
        caplog.clear()
        
    def test_fail_csv_format(self, monkeypatch, caplog, log_setup,
                             motif_error_params : str):
        #setup test
        def bad_format_df(path : str, sep : str) -> pd.DataFrame:
            data = [
                ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', -1, -1, 15.54],
                ]
            return pd.DataFrame(data, 
                                columns = ['Organism', 'Scaffold', 'Start', 'End', 
                                           'Score'])
        monkeypatch.setattr('pandas.read_csv', 
                            bad_format_df)
        #run test
        with pytest.raises(ValueError):
            with caplog.at_level(logging.INFO):
                collect(binary_path = 'a_binary_path', 
                        strict_span=False, 
                        neighbourhood_size=200, 
                        write_genomes = False, 
                        email = 'tdjkirkwood@hotmail.com', 
                        filenames = 'doesnt matter', 
                        results_dir = 'a_dir')
        assert caplog.messages == [motif_error_params, 
                                   "Unexpected binary file format\nexpected columns "\
                                       "- ['Organism', 'Scaffold', 'Start', 'End', "\
                                           "'Score', one or more Query Gene Names].\n"\
                                               "Actual columns: ['Organism', 'Scaffold', "\
                                                   "'Start', 'End', 'Score']"]
        caplog.clear()



    
class TestSieve:
    
    @pytest.fixture
    def sieve_setup(self, monkeypatch):
        def mock_all_vs_all_blast(input_genbank_dir : str, e_value : float, 
                                  max_target_seqs : int, output_blast_dir : str):
            pass
        
        @staticmethod
        def mock_build_neighbourhood_size_map(folder_with_genbanks : str) -> dict:
            return {'file1' : 2, 'file2' : 2, 'file3' : 3, 'file4' : 2}
        
        def mock_make_rbh_matrix(all_v_all_blast_xml : str, 
                                 min_percent_identity : int) -> dict:
            return {'file1' : {'0' : {'file1' : '0',
                                      'file2' : '0',
                                      'file3' : '0',
                                      'file4' : '0'},
                               '1' : {'file1' : '1',
                                      'file3' : '1',
                                      },
                               },
        
                    'file2' : {'0' : {'file1' : '0',
                                      'file2' : '0',
                                      'file3' : '0',
                                      'file4' : '0'},
                               '1' : {'file2' : '1'}
                               },
                    'file3' : {'0' : {'file1' : '0',
                                      'file2' : '0',
                                      'file3' : '0'},
                               '1' : {'file1' : '1',
                                      'file3' : '1'},
                               '2' : {'file3' : '2'}
                                       },
                    'file4' : {'0' : {'file1' : '0',
                                      'file2' : '0',
                                      'file4' : '0'},
                               
                               '1' : {'file4' : '1'}
                               }
                    }
        
        def mock_write_graph(graph : nx.Graph, path : str, similarity_filter : float, 
                             min_edge_view : float, logger_name : str):
            expected_edges = [('file1', 'file2', 0.5),
                              ('file1', 'file3', 1),
                              ('file1', 'file4', 0.5),
                              ('file2', 'file3', 0.5),
                              ('file2', 'file4', 0.5)
                              ]
            checked_edges = []
            
            #check all edges have expected weight
            for n1, n2, weight in expected_edges:
                if graph.has_edge(n1, n2):
                    assert graph.get_edge_data(n1, n2)['weight'] == weight
                    checked_edges += [(n1, n2)]
                else:
                    assert graph.get_edge_data(n2, n1)['weight'] == weight
                    checked_edges += [(n2, n1)]
                    
            #check there are no unexpected edges
            assert sorted(checked_edges) == sorted(graph.edges())
            logger = logging.getLogger(logger_name)
            logger.info(f'Made RBH graph of {len(graph.nodes)} unpruned '\
                            f'neighbourhoods - written to {path}')
                
        def mock_write_hist(graph : nx.Graph, path : str, logger_name : str):
            logger = logging.getLogger(logger_name)
            logger.info('Made histogram showing the distribution of RBH similarities '\
                            f'between all {len(graph.nodes)} neighbourhoods - written '\
                                f'to {path}')
        
        #general mocks
        monkeypatch.setattr('os.makedirs', 
                            mock_makedirs)
        monkeypatch.setattr('SyntenyQC.helpers.get_gbk_files', 
                            mock_get_gbk_files)
        
        #class mocks
        monkeypatch.setattr('SyntenyQC.pipelines.all_vs_all_blast', 
                            mock_all_vs_all_blast)
        monkeypatch.setattr('SyntenyQC.pipelines.make_rbh_matrix', 
                            mock_make_rbh_matrix)
        monkeypatch.setattr('SyntenyQC.pipelines.PrunedGraphWriter.build_neighbourhood_size_map', 
                            mock_build_neighbourhood_size_map)
        monkeypatch.setattr('SyntenyQC.pipelines.write_graph', 
                            mock_write_graph)
        monkeypatch.setattr('SyntenyQC.pipelines.write_hist', 
                            mock_write_hist)
    
    @pytest.fixture 
    def params(self) -> list:
        #Note - as these params have fairly minimal impact on logging behaviour 
        #(unlike the collect flags -sp and -wg) I am not checking many flag permutations 
        return ['---PARAMETERS---\nCommand: sieve\ninput_genbank_dir: a_gbk_folder\n'\
                       'e_value: 0.1\nmin_percent_identity: 50\nsimilarity_filter: 0.5\n'\
                           'results_dir: a_results_dir\noutput_blast_dir: a_results_dir\\blastp\n'\
                               'output_genbank_dir: a_results_dir\\genbank\n'\
                                   'output_vis_dir: a_results_dir\\visualisations\n'\
                                       'min_edge_view: 0.1\n\n\n']
                                       
        
    
        
    def test_write_fail(self, monkeypatch, caplog, params : list, sieve_setup,
                        log_setup):
        #setup test
        @staticmethod
        def mock_pruned_graph_write_nodes_fail(nodes : list, 
                                               input_genbank_folder : str, 
                                               output_genbank_folder : str) -> list:
            assert nodes == ['file3', 'file4']
            return ['this_node_is_new']
        monkeypatch.setattr('SyntenyQC.pipelines.PrunedGraphWriter.write_nodes', 
                            mock_pruned_graph_write_nodes_fail)
        #run test
        with pytest.raises(ValueError):
            with caplog.at_level(logging.INFO):
                sieve(input_genbank_dir = 'a_gbk_folder', 
                      e_value = 0.1, 
                      min_percent_identity = 50, 
                      max_target_seqs = 200,                      
                      similarity_filter = 0.5,
                      results_dir = 'a_results_dir',
                      min_edge_view = 0.1)
        fail_log = ['WRITTEN nodes and PRUNED node names dont match\n'\
                             "PRUNED NODES:\n['file3', 'file4']\n\n"\
                                     "WRITTEN NODES:\n['this_node_is_new']"]
        assert caplog.messages == params + fail_log
        caplog.clear()
    
    def test_sieve(self, monkeypatch, caplog, params : list, sieve_setup, 
                   log_setup):
        #setup test
        @staticmethod
        def mock_pruned_graph_write_nodes(nodes : list, input_genbank_dir : str, 
                                          output_genbank_dir : str) -> list:
            assert nodes == ['file3', 'file4']
            return nodes
        monkeypatch.setattr('SyntenyQC.pipelines.PrunedGraphWriter.write_nodes', 
                            mock_pruned_graph_write_nodes)
        #run test
        with caplog.at_level(logging.INFO):
            sieve(input_genbank_dir = 'a_gbk_folder', 
                  e_value = 0.1, 
                  min_percent_identity = 50, 
                  max_target_seqs = 200,                 
                  similarity_filter = 0.5,
                  results_dir = 'a_results_dir',
                  min_edge_view = 0.1)
            
        prune_string = 'Pruned graph - written 2 out of 4 initial neighbourhoods '\
                            'to a_results_dir\\genbank'
        graph_string = 'Made RBH graph of 4 unpruned neighbourhoods - written to '\
                           'a_results_dir\\visualisations\\RBH_graph.html'
        hist_string = 'Made histogram showing the distribution of RBH similarities '\
                          'between all 4 neighbourhoods - written to '\
                              'a_results_dir\\visualisations\\RBH_histogram.html'
        assert caplog.messages == params + [prune_string, graph_string, hist_string]
        caplog.clear()
