# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:31:24 2024

@author: u03132tk
"""
# import sys
# sys.path.insert(0, 'D:\\')
from SyntenyQC.networks import PrunedGraphWriter
#from networks import PrunedGraph
from general_mocks import mock_read_good_gbk, mock_get_gbk_files
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import pytest
import networkx as nx
def test_build_neighbourhood_size_map(monkeypatch):
    monkeypatch.setattr('SyntenyQC.networks.get_gbk_files', mock_get_gbk_files)
    monkeypatch.setattr('SyntenyQC.networks.read_gbk', mock_read_good_gbk)
    assert PrunedGraphWriter.build_neighbourhood_size_map('a_folder') == {'file1' : 2,
                                                                    'file2' : 2}

def test_build_neighbourhood_size_map_exceptions(monkeypatch):
    def mock_read_empty_gbk(path):
        return SeqRecord(seq = None,
                         features = [SeqFeature(type="not_CDS",
                                                qualifiers = {'translation' : ['not ',
                                                                               'a ', 
                                                                               'protein2 ', 
                                                                               'sequence']}
                                                )
                                     ]
                         )
    monkeypatch.setattr('SyntenyQC.networks.get_gbk_files', mock_get_gbk_files)
    monkeypatch.setattr('SyntenyQC.networks.read_gbk', mock_read_empty_gbk)
    with pytest.raises(ValueError):
        PrunedGraphWriter.build_neighbourhood_size_map('a_folder')
        
@pytest.fixture
def reciprocal_best_hit_matrix():
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
def test_make_graph(neighbourhood_size_map,min_edge_view,expected_edges, reciprocal_best_hit_matrix):
    #note, cannot directly comapre graphs - no nx.Graph.__eq__()
    G = PrunedGraphWriter.make_graph(reciprocal_best_hit_matrix, neighbourhood_size_map, min_edge_view)
    assert sorted(G.nodes()) == sorted(neighbourhood_size_map.keys())
    checked_edges = []
    #check all edges have expected weight
    for n1, n2, weight in expected_edges:
        if G.has_edge(n1, n2):
            assert G.get_edge_data(n1, n2)['weight'] == weight
            checked_edges += [(n1, n2)]
        else:
            assert G.get_edge_data(n2, n1)['weight'] == weight
            checked_edges += [(n2, n1)]
    #check there are no unexpected edges
    assert sorted(checked_edges) == sorted(G.edges())
@pytest.fixture
def similarity_graph():
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

def test_prune_graph(similarity_graph):
    #TODO test tie breaks (nodes with same degree > 0).  Should dleete one and leave the other 
    assert PrunedGraphWriter.prune_graph(similarity_graph, 1) == ['node1','node2','node4']
    #depends on order of node deletion - in this graph structure nodes 2-4 are considered equivalent i.e. interchangeable
    assert PrunedGraphWriter.prune_graph(similarity_graph, 0.5) in [['node1','node2'], 
                                                              ['node1','node3'], 
                                                              ['node1','node4']
                                                              ]

def test_write_nodes(monkeypatch):
    def mock_copy(src, dst):
        nonlocal monkeypatch_params
        monkeypatch_params['src'] += [src]
        monkeypatch_params['dst'] += [dst]
        return None
    monkeypatch_params = {'src' : [],
                          'dst' : []}
    monkeypatch.setattr('SyntenyQC.networks.get_gbk_files', mock_get_gbk_files)
    monkeypatch.setattr('shutil.copy', mock_copy)
    nodes = ['ignore1', 'ignore2', 'file1', 'file2', 'file1.ignore', 'file2.ignore']
    assert PrunedGraphWriter.write_nodes(nodes, 'a_genbank_folder', 'a_results_folder') == ['file1', 'file2']
    assert monkeypatch_params['src'] == ['a_genbank_folder\\file1.gbk', 'a_genbank_folder\\file2.gb']
    assert monkeypatch_params['dst'] == ['a_results_folder\\file1.gbk', 'a_results_folder\\file2.gb']
    #output nodes in and out of genbank files.  check shutil.copy has nodes in gbk files being copies.  check written nodes  
