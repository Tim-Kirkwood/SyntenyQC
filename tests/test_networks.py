# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:31:24 2024

@author: u03132tk
"""
from SyntenyQC.networks import PrunedGraphWriter
from general_mocks import mock_read_good_gbk, mock_get_gbk_files
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import pytest
import networkx as nx
import logging


@pytest.fixture
def reciprocal_best_hit_matrix() -> dict:
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

@pytest.fixture
def setup_get_gbk_files(monkeypatch):
    monkeypatch.setattr('SyntenyQC.networks.get_gbk_files', 
                        mock_get_gbk_files)
    
class TestBuildNeibourhoodSizeMap:
    
    @pytest.fixture
    def setup_read_good_gbk(self, monkeypatch, setup_get_gbk_files):   
        monkeypatch.setattr('SyntenyQC.networks.read_gbk', 
                            mock_read_good_gbk)
        
    @pytest.fixture
    def setup_read_empty_gbk(self, monkeypatch, setup_get_gbk_files):  
        def mock_read_empty_gbk(path : str) -> SeqRecord:
            feature = SeqFeature(type="not_CDS",
                                 qualifiers = {'translation' : ['not ',
                                                                'a ', 
                                                                'protein ', 
                                                                'sequence']
                                               }
                                 )
            return SeqRecord(seq = None,
                             features = [feature])
        monkeypatch.setattr('SyntenyQC.networks.read_gbk', 
                            mock_read_empty_gbk)
        
    def test_normal(self, setup_read_good_gbk):
        size_map = PrunedGraphWriter.build_neighbourhood_size_map('a_folder')
        assert size_map == {'file1' : 2,
                            'file2' : 2}
    
    def test_exceptions(self, setup_read_empty_gbk):
        with pytest.raises(ValueError):
            PrunedGraphWriter.build_neighbourhood_size_map('a_folder')

class TestMakeGraph:
    
    def has_equal_edges(self, graph : nx.Graph, expected_edges : list) -> bool:
        '''
        Cannot directly comapre nx graphs - no nx.Graph.__eq__() method.  So this 
        function compares graph edges and says if all edges (n1, n2, weight) 
        are in graph and graph has no other edges.  Note, this will not check 
        singlet nodes (which have no edges) so this function should be accompanied 
        by an explicit check for expected nodes if checking graph equality.

        Parameters
        ----------
        graph : nx.Graph
            Graph being assessed.
        expected_edges : list
            All edges expected in graph - format: [(n1, n2, weight), ...].

        Returns
        -------
        bool
            All edges in expected_edges are in graph.

        '''
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
        return True
    
    @pytest.mark.parametrize('neighbourhood_size_map,min_edge_view,expected_edges', 
                             [({'file1' : 2, 'file2' : 2, 'file3' : 3, 'file4' : 2},
                               0.5,
                               [('file1', 'file2', 0.5),
                                ('file1', 'file3', 1),
                                ('file1', 'file4', 0.5),
                                ('file2', 'file3', 0.5),
                                ('file2', 'file4', 0.5)
                                ]
                               ),
                              ({'file1' : 2, 'file2' : 2, 'file3' : 3, 'file4' : 2},
                               1,
                               [('file1', 'file3', 1)
                                ]
                               ),
                              ({'file1' : 5, 'file2' : 5, 'file3' : 5, 'file4' : 5},
                               0.2,                                                    
                               [('file1', 'file2', 0.2),
                                ('file1', 'file3', 0.4),
                                ('file1', 'file4', 0.2),
                                ('file2', 'file3', 0.2),
                                ('file2', 'file4', 0.2)
                                ]
                               ),
                              ({'file1' : 5, 'file2' : 5, 'file3' : 5, 'file4' : 5},
                               0.4,
                               [('file1', 'file3', 0.4)
                                ]
                               ),
                              ({'file1' : 5, 'file2' : 5, 'file3' : 5, 'file4' : 5},
                               0.6,
                               []
                               ),
                              ({'file1' : 4, 'file2' : 3, 'file3' : 2, 'file4' : 1},
                               0.3,
                               [('file1', 'file2', 1/3),
                                ('file1', 'file3', 1),
                                ('file1', 'file4', 1),
                                ('file2', 'file3', 0.5),
                                ('file2', 'file4', 1)
                                ]
                               )
                              ]
                         )
    def test_run(self, neighbourhood_size_map : dict, min_edge_view : float, 
                 expected_edges : list, reciprocal_best_hit_matrix : dict):
        G = PrunedGraphWriter.make_graph(reciprocal_best_hit_matrix, 
                                         neighbourhood_size_map, 
                                         min_edge_view)
        assert sorted(G.nodes()) == sorted(neighbourhood_size_map.keys())
        #cannot directly comapre graphs - no nx.Graph.__eq__()
        assert self.has_equal_edges(G, expected_edges)

class TestPruneGraph:
    @pytest.fixture
    def similarity_graph(self) -> nx.Graph:
        G = nx.Graph()
        G.add_nodes_from([f'node{i}' for i in range(1,6)])
        G.add_weighted_edges_from([('node1', 'node5', 1),
                                   ('node2', 'node3', 1),
                                   ('node2', 'node4', 0.5),
                                   ('node2', 'node5', 1),
                                   ('node3', 'node4', 1),
                                   ('node3', 'node5', 1),
                                   ('node4', 'node5', 1)
                                   ]
                                  )
        return G
    
    def test_prune_graph(self, similarity_graph : nx.Graph):
        pruned_graph = PrunedGraphWriter.prune_graph(similarity_graph, 
                                                     1)
        assert pruned_graph == ['node1','node2','node4']

        #Tie break - in this graph structure: 
        #    - Node 5 is the most highly connected and so is always deleted
        #    - Nodes 2-4 are considered equivalent as all degree 3. 
        #      Deleted node is purely dependent on order of nodes in graph.nodes
        pruned_graph = PrunedGraphWriter.prune_graph(similarity_graph, 
                                                     0.5)
        assert pruned_graph in [['node1','node2'], 
                                ['node1','node3'], 
                                ['node1','node4']
                                ]

class TestWriteNodes:
    
    @pytest.fixture 
    def setup_track_copy(self, monkeypatch, setup_get_gbk_files) -> dict:
        '''setup a mocked shutil.copy that records supplied src and dst in a 
        non-local dict, returned by this function'''
        def mock_copy(src : str, dst : str):
            nonlocal monkeypatch_params
            monkeypatch_params['src'] += [src]
            monkeypatch_params['dst'] += [dst]
        monkeypatch_params = {'src' : [],
                              'dst' : []}
        monkeypatch.setattr('shutil.copy', 
                            mock_copy)
        return monkeypatch_params
    
    def test_run(self, setup_track_copy):        
        nodes = ['file1', 'file2',
                 #should ignore any files that are not a node
                 'ignore1', 'ignore2', 'file1.ignore', 'file2.ignore']
        written_nodes = PrunedGraphWriter.write_nodes(nodes, 
                                             'a_genbank_folder', 
                                             'a_results_folder')
        assert written_nodes == ['file1', 'file2']
        assert setup_track_copy['src'] == ['a_genbank_folder\\file1.gbk', 
                                           'a_genbank_folder\\file2.gb']
        assert setup_track_copy['dst'] == ['a_results_folder\\file1.gbk', 
                                           'a_results_folder\\file2.gb']

class TestLogResults:
    
    def test_run(log_setup, caplog):
        with caplog.at_level(logging.INFO):
            PrunedGraphWriter.log_results(nodes_to_write = ['node1', 'node2'], 
                                          written_nodes = ['node2', 'node1'], 
                                          number_of_nodes = 3, 
                                          output_genbank_dir = 'output_path', 
                                          logger_name = 'collect')
        assert caplog.messages == ['Pruned graph - written 2 out of 3 initial '\
                                     'neighbourhoods to output_path']
        
    def test_no_nodes_to_write(log_setup, caplog):
        with caplog.at_level(logging.INFO):
            with pytest.raises(ValueError):
                PrunedGraphWriter.log_results(nodes_to_write = [], 
                                              written_nodes = ['node1', 'node2'], 
                                              number_of_nodes = 3, 
                                              output_genbank_dir = 'output_path', 
                                              logger_name = 'collect')
        assert caplog.messages == ['No nodes were available to write']
        
    def test_no_nodes_written(log_setup, caplog):
        with caplog.at_level(logging.INFO):
            with pytest.raises(ValueError):
                PrunedGraphWriter.log_results(nodes_to_write = ['node1', 'node2'], 
                                              written_nodes = [], 
                                              number_of_nodes = 3, 
                                              output_genbank_dir = 'output_path', 
                                              logger_name = 'collect')
        assert caplog.messages == ['No nodes were written']
        
    def test_different_nodes_written(log_setup, caplog):
        with caplog.at_level(logging.INFO):
            with pytest.raises(ValueError):
                PrunedGraphWriter.log_results(nodes_to_write = ['node1', 'node2'], 
                                              written_nodes = ['node1', 'node2', 'node3'], 
                                              number_of_nodes = 3, 
                                              output_genbank_dir = 'output_path', 
                                              logger_name = 'collect')
                
        message = 'WRITTEN nodes and PRUNED node names dont match\n'\
                     "PRUNED NODES:\n['node1', 'node2']\n\n"\
                         "WRITTEN NODES:\n['node1', 'node2', 'node3']"
        assert caplog.messages == [message]
            
    def test_different_nodes_write(log_setup, caplog):
        with caplog.at_level(logging.INFO):
            with pytest.raises(ValueError):
                PrunedGraphWriter.log_results(nodes_to_write = ['node1', 'node2', 'node3'], 
                                              written_nodes = ['node1', 'node2'], 
                                              number_of_nodes = 3, 
                                              output_genbank_dir = 'output_path', 
                                              logger_name = 'collect')
                
        message = 'WRITTEN nodes and PRUNED node names dont match\n'\
                     "PRUNED NODES:\n['node1', 'node2', 'node3']\n\n"\
                         "WRITTEN NODES:\n['node1', 'node2']"
        assert caplog.messages == [message]


class TestPrunedGraphWriter:

    @pytest.fixture
    def function_setup(self, monkeypatch, log_setup):
        @staticmethod
        def mock_build_neighbourhood_size_map(genbank_folder : str) -> dict:
            return {'file1' : 2, 'file2' : 2, 'file3' : 3, 'file4' : 2}
        def mock_copy(src : str, dst : str):
            pass 
        monkeypatch.setattr('shutil.copy', 
                            mock_copy)
        monkeypatch.setattr('SyntenyQC.networks.PrunedGraphWriter.build_neighbourhood_size_map',
                            mock_build_neighbourhood_size_map)
        
    @pytest.fixture
    def setup_nodes_written_ok(self, function_setup, monkeypatch):
        @staticmethod
        def mock_write_nodes(nodes : list, input_genbank_dir : str, 
                             output_genbank_dir : str) -> list:
            return nodes
        monkeypatch.setattr('SyntenyQC.networks.PrunedGraphWriter.write_nodes',
                            mock_write_nodes)
    
    @pytest.fixture
    def setup_no_nodes_to_write(self, function_setup, monkeypatch):
        @staticmethod
        def mock_no_pruned_nodes(graph : nx.Graph, similarity_filter : str) -> list:
            return []
        @staticmethod
        def mock_write_nodes(nodes : list, input_genbank_dir : str, 
                             output_genbank_dir : str) -> list:
            return ['file3', 'file4']
        monkeypatch.setattr('SyntenyQC.networks.PrunedGraphWriter.prune_graph',
                            mock_no_pruned_nodes)
        monkeypatch.setattr('SyntenyQC.networks.PrunedGraphWriter.write_nodes',
                            mock_write_nodes)
        
    @pytest.fixture
    def setup_no_nodes_written(self, function_setup, monkeypatch):
        @staticmethod
        def mock_no_written_nodes(nodes : list, input_genbank_dir : str, 
                                  output_genbank_dir : str) -> list:
            return []
        monkeypatch.setattr('SyntenyQC.networks.PrunedGraphWriter.write_nodes',
                            mock_no_written_nodes)
        
    @pytest.fixture
    def setup_different_nodes_written(self, function_setup, monkeypatch):
        @staticmethod
        def mock_different_written_nodes(nodes : list, input_genbank_dir : str, 
                                         output_genbank_dir : str) -> list:
            return ['file1', 'file2', 'file3', 'file4']
        #note the function_setup mock_copy() becomes irrelevant now
        monkeypatch.setattr('SyntenyQC.networks.PrunedGraphWriter.write_nodes',
                            mock_different_written_nodes)
    
    def test_run(self, reciprocal_best_hit_matrix : dict, 
                 setup_nodes_written_ok, caplog):
        with caplog.at_level(logging.INFO): 
            obj = PrunedGraphWriter(input_genbank_dir = 'doesnt matter', 
                                    reciprocal_best_hit_matrix = reciprocal_best_hit_matrix, 
                                    similarity_filter = 0.5,
                                    min_edge_view = 0.5,
                                    output_genbank_dir = 'output_dir',
                                    logger_name = 'collect')
        assert obj.nodes == ['file3', 'file4']
        message = 'Pruned graph - written 2 out of 4 initial neighbourhoods to output_dir'
        assert caplog.messages == [message]
    
    def test_no_nodes_to_write(self, reciprocal_best_hit_matrix : dict, 
                               setup_no_nodes_to_write, caplog):
        with pytest.raises(ValueError): 
            with caplog.at_level(logging.INFO): 
                PrunedGraphWriter(input_genbank_dir = 'doesnt matter', 
                                        reciprocal_best_hit_matrix = reciprocal_best_hit_matrix, 
                                        similarity_filter = 0.5,
                                        min_edge_view = 0.5,
                                        output_genbank_dir = 'output_dir',
                                        logger_name = 'collect')
        assert caplog.messages == ['No nodes were available to write']
    
    def test_no_nodes_written(self, reciprocal_best_hit_matrix : dict, 
                              setup_no_nodes_written, caplog):
        with pytest.raises(ValueError): 
            with caplog.at_level(logging.INFO): 
                PrunedGraphWriter(input_genbank_dir = 'doesnt matter', 
                                        reciprocal_best_hit_matrix = reciprocal_best_hit_matrix, 
                                        similarity_filter = 0.5,
                                        min_edge_view = 0.5,
                                        output_genbank_dir = 'output_dir',
                                        logger_name = 'collect')
        assert caplog.messages == ['No nodes were written']
        
    def test_different_nodes_written(self, reciprocal_best_hit_matrix : dict, 
                                     setup_different_nodes_written, caplog):
        with pytest.raises(ValueError): 
            with caplog.at_level(logging.INFO): 
                PrunedGraphWriter(input_genbank_dir = 'doesnt matter', 
                                        reciprocal_best_hit_matrix = reciprocal_best_hit_matrix, 
                                        similarity_filter = 0.5,
                                        min_edge_view = 0.5,
                                        output_genbank_dir = 'output_dir',
                                        logger_name = 'collect')
        message = 'WRITTEN nodes and PRUNED node names dont match\n'\
                     "PRUNED NODES:\n['file3', 'file4']\n\n"\
                         "WRITTEN NODES:\n['file1', 'file2', 'file3', 'file4']"
        assert caplog.messages == [message]