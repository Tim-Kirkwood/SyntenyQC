# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 17:01:26 2024

@author: u03132tk
"""
from SyntenyQC.blast_functions import all_vs_all_blast, make_rbh_matrix 
from SyntenyQC.networks import PrunedGraphWriter 
from SyntenyQC.visualisations import write_graph, write_hist 
import logging
import pandas as pd
from SyntenyQC.neighbourhood import Neighbourhood, write_results
from Bio import Entrez
from SyntenyQC.helpers import check_motif, get_gbk_files

'''
This module outlines the collect() and sieve() pipelines
'''

def get_params_string(local_vars : dict, command : str) -> str:
    '''
    Return a string showing the command and associated parameters for a given 
    collect() call.

    Parameters
    ----------
    command : str
        Which command was called (collect or sieve).

    local_vars : dict
        dict of local variables in env that get_params_string is called 
        (must include log keys of given command).

    Raises
    ------
    KeyError
        Command not recognised or variable not defined in local scope (likely 
        due to not being in parser).

    Returns
    -------
    str
        A string outlineing parameters used.

    '''
    #local_vars = locals()
    if command == 'collect':
        log_keys = ['binary_path', 'strict_span', 'neighbourhood_size', 
                    'write_genomes', 'email', 'filenames', 'results_dir']
    elif command == 'sieve':
        log_keys = ['genbank_folder', 'e_value', 'min_percent_identity', 
                    'similarity_filter', 'results_dir', 'min_edge_view']
    else:
        raise KeyError(f'Command {command} must be collect or sieve')
    log_string = f'Command: {command}\n'
    for log_key in log_keys:
        try:
            log_string += f'{log_key}: {local_vars[log_key]}\n'
        except KeyError:
            raise KeyError(f'{log_key} not in locals - {locals()}')
    return log_string

def get_log_handlers(filepath : str) -> list:
    '''
    Return list of handlers for logging.

    Parameters
    ----------
    filepath : str
        Log filepath.

    Returns
    -------
    list
        A file handler (so there is a local log file for posterity) and a stream 
        handler (so logs are reported at the command line for a given run).
    '''
    return [logging.FileHandler(filepath, 
                                'w'),
            logging.StreamHandler()
            ]

def sieve(genbank_folder : str, 
          e_value : float, 
          min_percent_identity : int, 
          max_target_seqs : int,
           
          similarity_filter : float,
          results_dir : str,
          min_edge_view : float) -> str:
    '''
    Run sieve pipeline to filter redundant neighbourhoods from genbank_folder.
    Redundant neighbourhoods share >= similarity_filter proportion of BLASTP 
    reciprocal best hits.  Hit alignmentss must have >= min_percent_identity 
    and <= e_value.  Only the top max_target_seqs hits are considered for a given 
    query protein. Unfiltered neighbourhoods are written in gbk format to results_dir.
    A network shouing all neighbourhoods and their pairwise similarities is 
    written to results_dir/RBH_graph.html (for similariites >= min_edge_view).
    A histogram showing the distribution of all edge weights > min_edge_view is 
    written to results_dir/RBH_histogram.html.  

    Parameters
    ----------
    genbank_folder : str
        Folder with genbanks files represnting neighburhoods to filter.
    e_value : float
        BLASTP alignment evalue threshold.
    min_percent_identity : int
        BLASTP alignment percent identity threshold.
    max_target_seqs : int
        BLASTP max_target_Seqs parameter.
    similarity_filter : float
        Filtering similarity threshold.
    results_dir : str
        Folder for results.
    min_edge_view : float
        Minimum edge weight to show in graph/histogram visualisations. If 
        min_edge_view < similarity_filter, edges will be red 
        (>= similarity_filter) or black (< similarity_filter).  Otherwise, edges 
        will be black.  This setting is purely used for viusalisation and has no 
        impact on graph pruning.

    Raises
    ------
    ValueError
        If there are no genbank files in genbank_folder or there is an issue 
        when writing the unfiltered neighbourhoods to results_dir.

    Returns
    -------
    str
        Path to results dir.

    '''
    
    #Set up logging
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s:  %(message)s',
                        handlers=get_log_handlers(f'{results_dir}\\log.txt')
                        )
    logger = logging.getLogger()
    
    #Write params to log file (can be used to debug if the run fails)
    params = get_params_string(local_vars = locals(),
                               command='sieve')
    logger.info(f'---PARAMETERS---\n{params}\n\n')
    
    #check genbank folder has >= 1 gbk file
    if get_gbk_files(genbank_folder) == []:
        error= f'No genbank (.gbk, .gb) files in {genbank_folder}'
        logger.error(error)
        raise ValueError(error)
        
    #perform RBH
    all_v_all_blast_xml = all_vs_all_blast(genbank_folder, 
                                           e_value,
                                           max_target_seqs)
    rbh_matrix = make_rbh_matrix(all_v_all_blast_xml, 
                                 min_percent_identity)
    
    #build and prune a similarity network from the RBH results, and copy acceptable 
    #neighbourhoods from genbank_folder to results_dir, and log associated information
    pruned_graph = PrunedGraphWriter(genbank_folder, 
                                     rbh_matrix,     
                                     similarity_filter,
                                     min_edge_view,
                                     results_dir)
    if sorted(pruned_graph.written_nodes) != sorted(pruned_graph.nodes):
        error = 'WRITTEN nodes and PRUNED node names dont match\n'\
                     'PRUNED NODES:\n'\
                         f'{pruned_graph.nodes}\n\n'\
                             'WRITTEN NODES:\n'\
                                 f'{pruned_graph.written_nodes}'
        logger.error(error)
        raise ValueError(error)    
    logger.info(f'Pruned graph - written {len(pruned_graph.nodes)} out of {len(pruned_graph.raw_graph.nodes)} initial neighbourhoods to {results_dir}')
    
    #write graph/histogram visualisations, save to graph/hist__path, and log associated information
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
    
    return results_dir
    

def collect(binary_path : str, 
            strict_span : bool,
            neighbourhood_size : int, 
            write_genomes : bool, 
            email : str,
            filenames : str,
            results_dir : str) -> str:
    '''
    Binary_file is a binary file from cblaster containing cblaster hit details 
    (record and locii of up/downstream terminii of first/last cblaster hits in 
    record).  This is parsed to define neighbourhood_size bp neighbourhoods, 
    which are written to local files if they can be scraped, written and extracted 
    properly.  Neighbourhood files are named according to either accession or 
    organism of parent record (which can also be written to local files if 
    write_genomes is True).  

    Parameters
    ----------
    binary_path : str
        Path to cblaster binary_path.
    strict_span : bool
        If True, reject neighbourhoods that have a terminus <0 or >length of 
        record.  Otherwise, correct to 0/length of record.
    neighbourhood_size : int
        Size of neighbourhood to build.
    write_genomes : bool
        If true, write parent record from which neighbourhood was defined to a 
        local file.
    email : str
        Required by NCBI for downlaoding files via entrez API.
    filenames : str
        Data to use for naming genome/neighbourhood files (organism or accession).
    results_dir : str
        Folder where results will be written.

    Raises
    ------
    ValueError
        If the cblaster binary file has a motif issue (which may in turn suggest 
        a data integrity issue and should not just be corrected as described for 
        strict_span).

    Returns
    -------
    str
        Folder where results are written.

    '''
    
    #set up logging 
    log_handlers= get_log_handlers(f'{results_dir}\\log.txt')
    logging.basicConfig(
                        level=logging.INFO,
                        format='%(asctime)s %(levelname)s:  %(message)s',
                        handlers=log_handlers
                    )
    logger = logging.getLogger()
    
    #NCBI requires a user email - set up before logging params so you can include 
    #in the log
    Entrez.email = email
    
    #Write params to log file (can be used to debug if the run fails)
    params = get_params_string(local_vars = locals(),
                               command='collect')
    logger.info(f'---PARAMETERS---\n{params}\n\n')
    
    #read in cblaster binary file, check it has expected format (5 columns plus
    #>=1 hit column) and parse aceession, motif start, motif end.  Define number 
    #of neighbourhoods for various print/log calls 
    location_data = pd.read_csv(binary_path, sep = ',')
    if len(location_data.columns)<=5:
        error = f"Unexpected binary file format\nexpected columns - ['Organism', 'Scaffold', 'Start', 'End', 'Score', one or more Query Gene Names].\nActual columns: {list(location_data.columns)}"
        logger.error(error)
        raise ValueError(error)   
    number_of_neighbourhoods = len(location_data['Scaffold'])
    logger.info(f"Extracting {number_of_neighbourhoods} gbks...")
    neighbourhood_data = zip (location_data['Scaffold'],
                             location_data['Start'],
                             location_data['End'])
    
    #Process each neithbourhood, keeping track of neighbourhood rejection causes 
    #for a final summary string. 
    writtten = 0
    skip_termini = 0
    skip_scrape = 0
    skip_motif = 0
    for index, (accession, motif_start, motif_stop) in enumerate(neighbourhood_data):
        
        #log aceesion, progress and motif details
        logger.info(f'accession {accession} #{index} of {number_of_neighbourhoods - 1}')  
        logger.info(f'motif = {motif_start} -> {motif_stop}')
        
        #sanity check motif read from cblaster file
        error_msg = check_motif(motif_start, motif_stop)
        if error_msg != '':
            logger.error(error_msg) 
            raise ValueError (error_msg)
            
        #define a Neighbourhood
        neighbourhood = Neighbourhood(accession, 
                                      motif_start, 
                                      motif_stop, 
                                      neighbourhood_size, 
                                      strict_span)
        
        #if the Neighbourhood could not be built, record why and continue 
        if neighbourhood.scrape_error != None:
            logger.warning(f'scrape_fail - {neighbourhood.scrape_error} - {accession}')
            skip_scrape += 1
            continue
        elif neighbourhood.overlapping_termini:
            logger.warning(f'overlapping_termini - {accession}')
            skip_termini += 1
            continue
        elif neighbourhood.motif_to_long:
            logger.warning(f'motif is longer than specified neighbourhood - {accession}')
            skip_motif += 1
            continue
        
        #log neighbourhood span and write neighbourhhood to local file with 
        #user-specified file name
        logger.info(f'neighbourhood = {neighbourhood.neighbourhood_start} -> '\
                        f'{neighbourhood.neighbourhood_stop}')
        path = write_results(results_folder=results_dir, 
                      neighbourhood=neighbourhood, 
                      filenames=filenames, 
                      scale = 'neighbourhood')
        logger.info(f'written {neighbourhood.accession} NEIGHBOURHOOD to {path}')
        
        #for a typical usecase (i.e. synteny plot), you only want neighbourhoods.
        #However, if you want the genome as well (e.g. to check for distal associateions 
        #shared between neighbourhoods) you can also get the entire record from 
        #which the neighbourhood was defined. THIS RECORD IS NOT CHECKED - it 
        #will be a contig buyt may not be an entire genome.
        if write_genomes:
            path = write_results(results_folder=results_dir, 
                          neighbourhood=neighbourhood, 
                          filenames=filenames, 
                          scale = 'genome')
            logger.info(f'written {neighbourhood.accession} GENOME to {path}')
        
        #track number of successful writes for final log entry
        writtten += 1
    
    #summarise run
    logger.info(f'Written {writtten} records to {results_dir} - {skip_scrape} '\
                    f'failed scraping, {skip_termini} had termini that exceeded the '\
                        f'genome boundaries, {skip_motif} had a motif that was too long'
                        )
    
    return results_dir

