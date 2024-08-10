# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 13:31:39 2024

@author: u03132tk
"""
import networkx as nx
import shutil
from SyntenyQC.helpers import get_gbk_files, read_gbk, get_protein_seq
import logging 
 
'''
this module contains code for RBH similarity network construction and 
manipulation - used by pipelines.sieve().

Note - should probably just make into a module...
'''
class PrunedGraphWriter:
    'Class that builds, prunes and writes files based on a similarity graph'
    def __init__(self, 
                 input_genbank_dir : str, 
                 reciprocal_best_hit_matrix : dict, 
                 similarity_filter : float,
                 min_edge_view : float,
                 output_genbank_dir : str,
                 logger_name : str) -> None:
        '''
        Build a network where gbk files in folder_with_genbanks are nodes and 
        edges are the number of shared RBHs between two files as a proportion 
        of the number of proteins in the file with the fewest proteins.  Prune 
        the graph, copy files corresponding to the remaining nodes from 
        folder_with_genbanks to results_folder.

        Parameters
        ----------
        input_genbank_folder : str
            Folder containing gbk files for sieveing (and maybe other file 
                                                      types/folders).
        reciprocal_best_hit_matrix : dict
            Reciprocal best hit matix from blast_functions.make_rbh_matrix().
            Format (numbers are CDS):
                    
                    {'file1' : {'0' : {'file1' : '0',
                                       'file2' : '0'},
                                '1' : {'file1' : '1'},
                                },

                     'file2' : {'0' : {'file1' : '0',
                                       'file2' : '0'},
                                '1' : {'file2' : '1'}
                                }
                     }
                    
        similarity_filter : float
            Minimum graph edge weight for two nodes to be considered equivalent.
        min_edge_view : float
            Minimmum edge weight for an edge to be included in the graph. This 
            setting is purely used for viusalisation and has no impact on graph 
            pruning.
        output_genbank_folder : str
            Folder to copy post-pruning genbank files into.
        logger_name : 
            name of logger to record results
        '''
        
        #get number of proteins in each gbk file
        neighbourhood_size_map = self.build_neighbourhood_size_map(input_genbank_dir)
        
        #make a similarity graph
        self.raw_graph = self.make_graph(reciprocal_best_hit_matrix, 
                                         neighbourhood_size_map, 
                                         min_edge_view
                                         )
        
        #prune the similarity graph
        self.nodes = self.prune_graph(self.raw_graph,
                                      similarity_filter
                                      )
        
        #copy the files that survive filtering
        self.written_nodes = self.write_nodes(self.nodes, 
                                              input_genbank_dir, 
                                              output_genbank_dir
                                              )
        
        #log summary of pruning outcome
        self.log_results(nodes_to_write = self.nodes,
                         written_nodes = self.written_nodes,
                         number_of_nodes = len(self.raw_graph.nodes),
                         output_genbank_dir = output_genbank_dir,
                         logger_name = logger_name)
        
        
    @staticmethod 
    def build_neighbourhood_size_map(genbank_folder : str) -> dict:
        '''
        Return a dict mapping filenames (keys) to the number of protein features 
        (values) in each file. 

        Parameters
        ----------
        genbank_folder : str
            Folder with genbank files(and potentialy other filetypes/folders).

        Raises
        ------
        ValueError
            >=1 file has no annotated proteins.

        Returns
        -------
        size_map : dict
            Dictionary mapping filenames to number of proteins.  E.g. if
            genbank_folder has one genbanl file called filename, and filename has  
            has 3 proteins, size_map will be {filename : 3}.

        '''
        genbank_files = get_gbk_files(genbank_folder)
        size_map = {}
        for file in genbank_files:
            
            #get a filename key
            record = read_gbk(f'{genbank_folder}\\{file}')
            record_id = file[0:file.rindex('.')]
            
            #initialise count
            size_map[record_id] = 0
            
            #increment count for features with a translation 
            for feature in record.features:
                if feature.type == 'CDS':
                    translation = get_protein_seq(feature)
                    if translation != '':
                        size_map[record_id] += 1
            
            #all files should have at least one protein 
            if size_map[record_id] == 0:
                raise ValueError(f'{file} has no annotated proteins')
                
        return size_map
    
    @staticmethod 
    def make_graph(reciprocal_best_hit_matrix : dict, 
                   neighbourhood_size_map : dict, 
                   min_edge_view : float) -> nx.Graph:
        '''
        Make a graph where nodes are file names and edges are the number of shared 
        RBH sequences as a proportion of the number of proteins in the file with 
        the fewest proteins.

        Parameters
        ----------
        reciprocal_best_hit_matrix : dict
            Reciprocal best hit matix from blast_functions.make_rbh_matrix().
            Format (numbers are CDS):
                    
                    {'file1' : {'0' : {'file1' : '0',
                                       'file2' : '0'},
                                '1' : {'file1' : '1'},
                                },

                     'file2' : {'0' : {'file1' : '0',
                                       'file2' : '0'},
                                '1' : {'file2' : '1'}
                                }
                     }
                    
        neighbourhood_size_map : dict
            Dictionary mapping filenames to number of proteins.  E.g. if
            genbank_folder has one genbanl file called filename, and filename has  
            has 3 proteins, size_map will be {filename : 3}.
        min_edge_view : float
            Minimmum edge weight for an edge to be included in the graph. This 
            setting is purely used for viusalisation and has no impact on graph 
            pruning.

        Returns
        -------
        G : nx.Graph
            nx.Graph object.

        '''
        
        #Empty edges to start with
        edge_map = {}
        for accession in reciprocal_best_hit_matrix.keys():
            edge_map[accession] = {}
            for other_accession in reciprocal_best_hit_matrix.keys():
                
                #assess edge between accession and other_accession if it is not 
                #a self-edge
                if accession != other_accession:
                    
                    #define the smallest number of preoteins 
                    record_size = neighbourhood_size_map[accession]
                    other_record_size = neighbourhood_size_map[other_accession]
                    smallest_record_protein_count = min(record_size, 
                                                        other_record_size)
                    
                    #count number of RBH between accession and other_accession
                    number_of_rbh = 0
                    for query_protein, rbh_scaffolds in reciprocal_best_hit_matrix[other_accession].items():
                        if accession in rbh_scaffolds.keys():
                            number_of_rbh += 1
                    
                    #define edge and add to edge map if it meeds min_edge_view threshold
                    similarity_score = number_of_rbh/smallest_record_protein_count
                    if similarity_score >= min_edge_view:
                        edge_map[accession][other_accession] = {'weight' : similarity_score} 
        
        #make a graph from the populated edge map
        G = nx.Graph(edge_map).to_undirected()
        
        #add any singlet nodes to the graph and return.  Note, every query will 
        #be in the RBH matrix, as it will be a best hit to itself at least.
        all_possible_nodes = set(reciprocal_best_hit_matrix.keys())
        G.add_nodes_from(all_possible_nodes)
        return G
    
    
    @staticmethod 
    def prune_graph(graph : nx.Graph, similarity_filter : float) -> list:
        '''
        Prune graph by removing most connected node and recomputing graph from
        remaining nodes, until no node degree is > 0. 

        Parameters
        ----------
        graph : nx.Graph
            Similarity graph.
        similarity_filter : float
            Minimum edge weight for two nodes to be considered equivalent and 
            thus prune-able.

        Returns
        -------
        nodes : list
            List of nodes in the final pruned graph.

        '''
        #initialise copy of graph - ignore edges less than similarity filter.  
        #Do not change input graph object.
        copy_graph = nx.Graph()
        copy_graph.add_nodes_from(graph.nodes)
        edges = [(n1, n2, e) for n1, n2, e in graph.edges.data('weight') if e>=similarity_filter]
        copy_graph.add_weighted_edges_from(edges)
        
        #initialise nodes
        nodes = list(copy_graph.nodes) 
        
        #prune
        while True:
            
            #make temporary subgraph of input graph using nodes 
            temp_G = copy_graph.subgraph(nodes)
            
            #get degrees of each nodein subgraph
            temp_nodes = []
            degrees = []
            for node, degree in temp_G.degree:
                temp_nodes += [node]
                degrees += [degree]
            
            #get maximum degree of nodes and exit if 0
            max_degree = max(degrees)
            if max_degree == 0:
                break
            
            #delete first node that has a degree equal to the maximum
            delete_index = degrees.index(max_degree)
            
            #reset nodes to nodes in the subgraph
            nodes = [n for i, n in enumerate(temp_nodes) if i != delete_index]
        
        return nodes
    
    
    @staticmethod 
    def write_nodes(nodes : list, input_genbank_dir : str, 
                    output_genbank_dir : str) -> list:
        '''
        Copy files with a filename in nodes from genbank_folder to results_folder 

        Parameters
        ----------
        nodes : list
            Target filenames to copy.
        input_genbank_dir : str
            Source folder with files to copy.
        output_genbank_folder : str
            Destination folder of files to copy.

        Returns
        -------
        written_nodes : list
            List of files (with no suffix) that have been written.

        '''
        gbk_files = get_gbk_files(input_genbank_dir)
        written_nodes = []
        for file in gbk_files:
            
            #get a filename in nodes
            record_id = file[0:file.rindex('.')]
            if record_id in nodes:
                
                #copy file from genbank_foilder to results_folder
                #TODO add a check for the copied file being in the results_folder 
                #already (should not be possible though, results folder will be newly
                #made in app.py).
                this_path = f'{input_genbank_dir}\\{file}'
                new_path = f'{output_genbank_dir}\\{file}'
                shutil.copy(src = this_path, 
                            dst = new_path)
                
                #add written nodes to list for return
                written_nodes += [record_id]

        return written_nodes
        
    @staticmethod
    def log_results(nodes_to_write : list, written_nodes : list, 
                    number_of_nodes : int, output_genbank_dir : str, 
                    logger_name : str):
        '''
        Logs pruning summary and says where files are written.  If there is an 
        error, describes error before raising ValueError exception. 

        Parameters
        ----------
        nodes_to_write : list
            Nodes that survive graph pruning.
        written_nodes : list
            Nodes (files) that have been copied to output_genbank_dir.
        number_of_nodes : int
            Number of nodes in graph before pruning.
        output_genbank_dir : str
            Folder where genbank files that survive pruning are written.
        logger_name : str
            name of logger to record results.

        Raises
        ------
        ValueError
            - nodes_to_write != written_nodes
            - no nodes_to_write or written_nodes (should always be >= 1).



        '''
        logger = logging.getLogger(logger_name)
        if sorted(nodes_to_write) != sorted(written_nodes):
            if nodes_to_write == []:
                error = 'No nodes were available to write'
            elif written_nodes == []:
                error = 'No nodes were written'
            else:
                error = 'WRITTEN nodes and PRUNED node names dont match\n'\
                             f'PRUNED NODES:\n{nodes_to_write}\n\n'\
                                 f'WRITTEN NODES:\n{written_nodes}'
            logger.error(error)
            raise ValueError(error)   
        else:
            logger.info(f'Pruned graph - written {len(written_nodes)} out of '\
                            f'{number_of_nodes} initial neighbourhoods '\
                                f'to {output_genbank_dir}')