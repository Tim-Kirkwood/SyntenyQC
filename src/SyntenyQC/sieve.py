# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:08:28 2024

@author: u03132tk
"""
from SyntenyQC.blast_functions import all_vs_all_blast, make_rbh_matrix 
from SyntenyQC.networks import PrunedGraph 
from SyntenyQC.visualisations import write_graph, write_hist 
import logging
 
def sieve(genbank_folder, 
          e_value, 
          min_percent_identity, 
          max_target_seqs,
           
          similarity_filter,
          results_dir,
          min_edge_view):
    logging.basicConfig(
                        level=logging.INFO,
                        format='%(asctime)s %(levelname)s:  %(message)s',
                        handlers=[
                            logging.FileHandler(f'{results_dir}\\log.txt', 
                                                'w'),
                            logging.StreamHandler()
                        ]
                    )
    
    logger = logging.getLogger()

    logger.info(f'\n---PARAMETERS---\n{locals()}\n\n')
    all_v_all_blast_xml = all_vs_all_blast(genbank_folder, 
                                           e_value,
                                           max_target_seqs)
    rbh_matrix = make_rbh_matrix(all_v_all_blast_xml, 
                            min_percent_identity)
    pruned_graph = PrunedGraph(genbank_folder, 
                                   rbh_matrix, 
                                    
                                   similarity_filter,
                                   min_edge_view,
                                   results_dir)
    if sorted(pruned_graph.written_nodes) != sorted(pruned_graph.nodes):
        logger.warn('WRITTEN nodes and PRUNED node names dont match')
        logger.warn(f'PRUNED NODES:  \n{pruned_graph.nodes}\n\n')
        logger.warn(f'WRITTEN NODES:  \n{pruned_graph.nodes}')
        raise ValueError('WRITTEN nodes and PRUNED node names dont match')
    logger.info(f'Pruned graph - written {len(pruned_graph.nodes)} out of {len(pruned_graph.raw_graph.nodes)} initial neighbourhoods to {results_dir}')
    graph_path = f'{results_dir}\\RBH_graph.html'
    write_graph(pruned_graph.raw_graph, 
                          graph_path, 
                    similarity_filter,
                    min_edge_view)
    logger.info(f'Made RBH graph of {len(pruned_graph.raw_graph.nodes)} unpruned neighbourhoods - written to {graph_path}')
    hist_path = f'{results_dir}\\RBH_histogram.html'
    write_hist(pruned_graph.raw_graph,
                  hist_path)
    logger.info(f'Made histogram showing the distribution of RBH similarities between all {len(pruned_graph.raw_graph)} neighbourhoods - written to {hist_path}')
    
    
    #TODO - other blastp params
    #TODO - queries/other folders/filesneed to be written
    
    
    