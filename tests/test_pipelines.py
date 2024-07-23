# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 17:01:26 2024

@author: u03132tk
"""
import logging
import pandas as pd
import pytest 
from SyntenyQC.pipelines import collect, sieve
from general_mocks import mock_get_gbk_files
                        


                        
def mock_get_log_handlers(path):
        #nonlocal test_stream
        return [logging.StreamHandler()]
def test_collect(monkeypatch, caplog):
    def good_df(path, sep):
        data = [
            ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', 100, 200, 15.54, 1],#good
            ['Streptomyces anthocyanicus IPS92w', 'bad_accession', 100, 200, 15.54, 1],#scrape error
            ['Streptomyces anthocyanicus IPS92w', '', 100, 200, 15.54, 1],#empty accession
            ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', 1, 100, 15.54, 1],#overlapping termini
            ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', 1, 10**6, 15.54, 1]#motif too large - note will be extended from MIDPOINT of start/end - i.e. 5**6 
            ]
        return pd.DataFrame(data, 
                            columns = ['Organism', 'Scaffold', 'Start', 'End', 'Score', 'a_locus_tag']
                            )
                            
    def mock_write_results(results_folder, neighbourhood, filenames, scale):
        return f'path_to_{scale}.gbk'
    

    #test_stream = StringIO()
    monkeypatch.setattr('SyntenyQC.pipelines.get_log_handlers', mock_get_log_handlers)
    monkeypatch.setattr('pandas.read_csv', good_df)
    monkeypatch.setattr('SyntenyQC.pipelines.write_results', mock_write_results)
    with caplog.at_level(logging.INFO):
        output = collect(binary_path = 'a_binary_path', strict_span=True, neighbourhood_size=200, 
                         write_genomes = True, email = 'tdjkirkwood@hotmail.com', 
                         filenames = 'doesnt matter', results_dir = 'a_dir')
        assert output == 'a_dir'
        #print (caplog.messages)
        message = ['---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                       'strict_span: True\nneighbourhood_size: 200\nwrite_genomes: True\n'\
                           'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                               'results_dir: a_dir\n\n\n', 
                   'Extracting 5 gbks...', 
                   'accession NZ_JAPFQR010000057.1 #0 of 4', 
                   'motif = 100 -> 200', 'neighbourhood = 50 -> 250', 
                   'written NZ_JAPFQR010000057.1 NEIGHBOURHOOD to path_to_neighbourhood.gbk', 
                   'written NZ_JAPFQR010000057.1 GENOME to path_to_genome.gbk', 
                   'accession bad_accession #1 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - HTTP error - bad_accession', 
                   'accession  #2 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - ValueError - No genome found - ', 
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
    with caplog.at_level(logging.INFO):
        output = collect(binary_path = 'a_binary_path', strict_span=False, neighbourhood_size=200, 
                         write_genomes = True, email = 'tdjkirkwood@hotmail.com', 
                         filenames = 'doesnt matter', results_dir = 'a_dir')
        assert output == 'a_dir'
        #print (caplog.messages)
        message = ['---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                       'strict_span: False\nneighbourhood_size: 200\nwrite_genomes: True\n'\
                           'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                               'results_dir: a_dir\n\n\n', 
                   'Extracting 5 gbks...', 
                   'accession NZ_JAPFQR010000057.1 #0 of 4', 
                   'motif = 100 -> 200', 'neighbourhood = 50 -> 250', 
                   'written NZ_JAPFQR010000057.1 NEIGHBOURHOOD to path_to_neighbourhood.gbk', 
                   'written NZ_JAPFQR010000057.1 GENOME to path_to_genome.gbk', 
                   'accession bad_accession #1 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - HTTP error - bad_accession', 
                   'accession  #2 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - ValueError - No genome found - ', 
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
    with caplog.at_level(logging.INFO):
        output = collect(binary_path = 'a_binary_path', strict_span=True, neighbourhood_size=200, 
                         write_genomes = False, email = 'tdjkirkwood@hotmail.com', 
                         filenames = 'doesnt matter', results_dir = 'a_dir')
        assert output == 'a_dir'
        #print (caplog.messages)
        message = ['---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                       'strict_span: True\nneighbourhood_size: 200\nwrite_genomes: False\n'\
                           'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                               'results_dir: a_dir\n\n\n', 
                   'Extracting 5 gbks...', 
                   'accession NZ_JAPFQR010000057.1 #0 of 4', 
                   'motif = 100 -> 200', 'neighbourhood = 50 -> 250', 
                   'written NZ_JAPFQR010000057.1 NEIGHBOURHOOD to path_to_neighbourhood.gbk', 
                   #'written NZ_JAPFQR010000057.1 GENOME to path_to_genome.gbk', 
                   'accession bad_accession #1 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - HTTP error - bad_accession', 
                   'accession  #2 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - ValueError - No genome found - ', 
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
    with caplog.at_level(logging.INFO):
        output = collect(binary_path = 'a_binary_path', strict_span=False, neighbourhood_size=200, 
                         write_genomes = False, email = 'tdjkirkwood@hotmail.com', 
                         filenames = 'doesnt matter', results_dir = 'a_dir')
        assert output == 'a_dir'
        #print (caplog.messages)
        message = ['---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                       'strict_span: False\nneighbourhood_size: 200\nwrite_genomes: False\n'\
                           'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                               'results_dir: a_dir\n\n\n', 
                   'Extracting 5 gbks...', 
                   'accession NZ_JAPFQR010000057.1 #0 of 4', 
                   'motif = 100 -> 200', 'neighbourhood = 50 -> 250', 
                   'written NZ_JAPFQR010000057.1 NEIGHBOURHOOD to path_to_neighbourhood.gbk', 
                   #'written NZ_JAPFQR010000057.1 GENOME to path_to_genome.gbk', 
                   'accession bad_accession #1 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - HTTP error - bad_accession', 
                   'accession  #2 of 4', 
                   'motif = 100 -> 200', 
                   'scrape_fail - ValueError - No genome found - ', 
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
    
def test_collect_motif_exceptions(monkeypatch, caplog):
    def error_df_1(path, sep):
        
        data = [
            ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', 1, 1, 15.54, 1],
            ]
        return pd.DataFrame(data, 
                            columns = ['Organism', 'Scaffold', 'Start', 'End', 'Score', 'a_locus_tag']
                            )
    def error_df_2(path, sep):
        data = [
            ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', -1, 0, 15.54, 1],
            ]
        return pd.DataFrame(data, 
                            columns = ['Organism', 'Scaffold', 'Start', 'End', 'Score', 'query_locus_tag']
                            )
    def error_df_3(path, sep):
        data = [
            ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', -1, 1, 15.54, 1],
            ]
        return pd.DataFrame(data, 
                            columns = ['Organism', 'Scaffold', 'Start', 'End', 'Score', 'query_locus_tag']
                            )
    def error_df_4(path, sep):
        data = [
            ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', -1, -1, 15.54, 1],
            ]
        return pd.DataFrame(data, 
                            columns = ['Organism', 'Scaffold', 'Start', 'End', 'Score', 'query_locus_tag']
                            )
    
    monkeypatch.setattr('SyntenyQC.pipelines.get_log_handlers', mock_get_log_handlers)
    param_string = '---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                       'strict_span: False\nneighbourhood_size: 200\nwrite_genomes: False\n'\
                           'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                               'results_dir: a_dir\n\n\n'
    
    funcs = [error_df_1, error_df_2, error_df_3, error_df_4]
    start_stop = [(1, 1), (-1, 0), (-1, 1), (-1, -1)]
    messages = ['motif_start 1 is >= motif_stop 1\n',
                'motif_start -1 < 0\nmotif_stop 0 <= 0\n',
                'motif_start -1 < 0\n',
                'motif_start -1 < 0\nmotif_stop -1 <= 0\nmotif_start -1 is >= motif_stop -1\n'                
                ]
    
    
    for func, message, (motif_start, motif_stop) in zip(funcs, messages, start_stop):
        monkeypatch.setattr('pandas.read_csv', func)
        with pytest.raises(ValueError):
            with caplog.at_level(logging.INFO):
                collect(binary_path = 'a_binary_path', strict_span=False, neighbourhood_size=200, 
                        write_genomes = False, email = 'tdjkirkwood@hotmail.com', 
                        filenames = 'doesnt matter', results_dir = 'a_dir')
        assert caplog.messages == [param_string, 
                                   'Extracting 1 gbks...', 
                                   'accession NZ_JAPFQR010000057.1 #0 of 0',
                                   f'motif = {motif_start} -> {motif_stop}',
                                   message]
        caplog.clear()

def test_collect_binary_exceptions(monkeypatch, caplog):
    def bad_format_df(path, sep):
        data = [
            ['Streptomyces anthocyanicus IPS92w', 'NZ_JAPFQR010000057.1', -1, -1, 15.54],
            ]
        return pd.DataFrame(data, 
                            columns = ['Organism', 'Scaffold', 'Start', 'End', 'Score'])
    monkeypatch.setattr('SyntenyQC.pipelines.get_log_handlers', mock_get_log_handlers)
    monkeypatch.setattr('pandas.read_csv', bad_format_df)
    param_string = '---PARAMETERS---\nCommand: collect\nbinary_path: a_binary_path\n'\
                       'strict_span: False\nneighbourhood_size: 200\nwrite_genomes: False\n'\
                           'email: tdjkirkwood@hotmail.com\nfilenames: doesnt matter\n'\
                               'results_dir: a_dir\n\n\n'
    with pytest.raises(ValueError):
        with caplog.at_level(logging.INFO):
            collect(binary_path = 'a_binary_path', strict_span=False, neighbourhood_size=200, 
                    write_genomes = False, email = 'tdjkirkwood@hotmail.com', 
                    filenames = 'doesnt matter', results_dir = 'a_dir')
    assert caplog.messages == [param_string, 
                               "Unexpected binary file format\nexpected columns - ['Organism', 'Scaffold', 'Start', 'End', 'Score', one or more Query Gene Names].\nActual columns: ['Organism', 'Scaffold', 'Start', 'End', 'Score']"]
    caplog.clear()


        

def test_sieve_empty_gbk(monkeypatch, caplog):
    
    def mock_get_empty_gbk(path):
        return []
    #empty genbank foldr
    #does it fail 
    #does it log 
    monkeypatch.setattr('SyntenyQC.pipelines.get_log_handlers', mock_get_log_handlers)
    monkeypatch.setattr('SyntenyQC.pipelines.get_gbk_files', mock_get_empty_gbk)
    param_string = '---PARAMETERS---\nCommand: sieve\ngenbank_folder: a_gbk_folder\n'\
                       'e_value: 0.1\nmin_percent_identity: 50\nsimilarity_filter: 0.7\n'\
                           'results_dir: a_results_dir\nmin_edge_view: 0.6\n\n\n'
                               
    with pytest.raises(ValueError):
        with caplog.at_level(logging.INFO):
            sieve(genbank_folder = 'a_gbk_folder', 
                  e_value = 0.1, 
                  min_percent_identity = 50, 
                  max_target_seqs = 200,
                  
                  similarity_filter = 0.7,
                  results_dir = 'a_results_dir',
                  min_edge_view = 0.6)
    #print (caplog.messages)
    assert caplog.messages == [param_string, 
                               'No genbank (.gbk, .gb) files in a_gbk_folder']
    caplog.clear()

def mock_all_vs_all_blast(genbank_folder, 
                          e_value,
                          max_target_seqs):
    return None

def mock_build_neighbourhood_size_map(self, folder_with_genbanks):
    return {'file1' : 2, 'file2' : 2, 'file3' : 3, 'file4' : 2}
def mock_make_rbh_matrix(all_v_all_blast_xml, 
                         min_percent_identity):
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

def test_sieve_write_fail(monkeypatch, caplog):
    
    
    
    def mock_pruned_graph_write_nodes_fail(self, nodes, 
                                     folder_with_genbanks, 
                                     results_folder):
        return ['this_node_is_new']
    
    #empty genbank foldr
    #does it fail 
    #does it log 
    monkeypatch.setattr('SyntenyQC.pipelines.get_log_handlers', mock_get_log_handlers)
    monkeypatch.setattr('SyntenyQC.pipelines.get_gbk_files', mock_get_gbk_files)
    monkeypatch.setattr('SyntenyQC.pipelines.all_vs_all_blast', mock_all_vs_all_blast)
    monkeypatch.setattr('SyntenyQC.pipelines.make_rbh_matrix', mock_make_rbh_matrix)
    monkeypatch.setattr('SyntenyQC.pipelines.PrunedGraphWriter.write_nodes', 
                        mock_pruned_graph_write_nodes_fail)
    monkeypatch.setattr('SyntenyQC.pipelines.PrunedGraphWriter.build_neighbourhood_size_map', 
                        mock_build_neighbourhood_size_map)
    param_string = '---PARAMETERS---\nCommand: sieve\ngenbank_folder: a_gbk_folder\n'\
                       'e_value: 0.1\nmin_percent_identity: 50\nsimilarity_filter: 0.5\n'\
                           'results_dir: a_results_dir\nmin_edge_view: 0\n\n\n'
    error_message = 'WRITTEN nodes and PRUNED node names dont match\n'\
                         'PRUNED NODES:\n'\
                             "['file3', 'file4']\n\n"\
                                 'WRITTEN NODES:\n'\
                                     "['this_node_is_new']"
    with pytest.raises(ValueError):
        with caplog.at_level(logging.INFO):
            sieve(genbank_folder = 'a_gbk_folder', 
                  e_value = 0.1, 
                  min_percent_identity = 50, 
                  max_target_seqs = 200,
                  
                  similarity_filter = 0.5,
                  results_dir = 'a_results_dir',
                  min_edge_view = 0)
    assert caplog.messages == [param_string, 
                               error_message]
    caplog.clear()







def test_sieve(monkeypatch, caplog):
    def mock_pruned_graph_write_nodes(self, nodes, 
                                     folder_with_genbanks, 
                                     results_folder):
        assert nodes == ['file3', 'file4']
        return nodes
    def mock_write_graph(graph, 
                         graph_path, 
                         similarity_filter,
                         min_edge_view):
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
        
    def mock_write_hist(graph, 
                        path):
        return None
    
    #empty genbank foldr
    #does it fail 
    #does it log 
    monkeypatch.setattr('SyntenyQC.pipelines.get_log_handlers', 
                        mock_get_log_handlers)
    monkeypatch.setattr('SyntenyQC.pipelines.get_gbk_files', 
                        mock_get_gbk_files)
    monkeypatch.setattr('SyntenyQC.pipelines.all_vs_all_blast', 
                        mock_all_vs_all_blast)
    monkeypatch.setattr('SyntenyQC.pipelines.make_rbh_matrix', 
                        mock_make_rbh_matrix)
    monkeypatch.setattr('SyntenyQC.pipelines.PrunedGraphWriter.build_neighbourhood_size_map', 
                        mock_build_neighbourhood_size_map)
    monkeypatch.setattr('SyntenyQC.pipelines.PrunedGraphWriter.write_nodes', 
                        mock_pruned_graph_write_nodes)
    monkeypatch.setattr('SyntenyQC.pipelines.write_graph', 
                        mock_write_graph)
    monkeypatch.setattr('SyntenyQC.pipelines.write_hist', 
                        mock_write_hist)
    with caplog.at_level(logging.INFO):
        sieve(genbank_folder = 'a_gbk_folder', 
                          e_value = 0.1, 
                          min_percent_identity = 50, 
                          max_target_seqs = 200,
                          
                          similarity_filter = 0.5,
                          results_dir = 'a_results_dir',
                          min_edge_view = 0.1)
    #assert output.nodes == ['file3', 'file4'], output
    param_string = '---PARAMETERS---\nCommand: sieve\ngenbank_folder: a_gbk_folder\n'\
                       'e_value: 0.1\nmin_percent_identity: 50\nsimilarity_filter: 0.5\n'\
                           'results_dir: a_results_dir\nmin_edge_view: 0.1\n\n\n'
    prune_string = 'Pruned graph - written 2 out of 4 initial neighbourhoods to a_results_dir'
    graph_string = 'Made RBH graph of 4 unpruned neighbourhoods - written to '\
                       'a_results_dir\\RBH_graph.html'
    hist_string = 'Made histogram showing the distribution of RBH similarities '\
                      'between all 4 neighbourhoods - written to a_results_dir\\RBH_histogram.html'
    assert caplog.messages == [param_string, 
                               prune_string,
                               graph_string,
                               hist_string
                               ]
    caplog.clear()
