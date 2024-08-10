# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 20:36:03 2024

@author: u03132tk
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature 
import os
import logging
import pandas as pd

def read_gbk(path : str) -> SeqRecord:
    '''
    Read a genbank file to a biopython SeqRecord.

    Parameters
    ----------
    path : str
        Path to genbank file.

    Returns
    -------
    SeqRecord
        Biopython SeqRecord.

    '''
    with open(path, 'r') as handle:
        return SeqIO.read(handle, 
                          "genbank")

def get_gbk_files(folder : str) -> list:
    '''
    Return list of genbank filenames from folder.

    Parameters
    ----------
    folder : str
        Folder with genbank files (and perhaps other fies/folders).

    Returns
    -------
    gbk_files : list
        List of genbank file names.

    '''
    files = os.listdir(folder)
    suffixes = ['.gbk','.gb'] 
    gbk_files = []
    for file in files:
        #os.listdir also returns dirs, which wont have a suffix
        if os.path.isfile(f'{folder}\\{file}'):        
            if file[file.rindex('.') : ] in suffixes:
                gbk_files += [file]
    return gbk_files

def get_cds_count(genbank_folder : str) -> int:
    '''
    Count cds in all genbank files in genbank_folder

    Parameters
    ----------
    genbank_folder : str
        Folder with genbank files (and perhaps other fies/folders)..

    Returns
    -------
    int
        Number of cds across all genbank files in genbank_folder.

    '''
    genbank_files = get_gbk_files(genbank_folder)
    cds_count = 0
    for file in genbank_files:
        record = read_gbk(f'{genbank_folder}\\{file}')
        for f in record.features:
            if f.type == 'CDS':
                #not psuedo etc
                if get_protein_seq(f) != '':
                    cds_count += 1
    return cds_count


def check_motif(motif_start : int, motif_stop : int) -> str:
    '''
    Check that a cblaster motif meets expected formatting.

    Parameters
    ----------
    motif_start : int
        Genomic start of first cblaster hit.
    motif_stop : int
        Genomic end of last cblaster hit.

    Returns
    -------
    msg : str
        Error message ('' if no errors).

    '''
    msg = ''
    if motif_start < 0:
        msg += f'motif_start {motif_start} < 0\n'
    if motif_stop <= 0:
        msg += f'motif_stop {motif_stop} <= 0\n'
    if motif_start >= motif_stop:
        msg += f'motif_start {motif_start} is >= motif_stop {motif_stop}\n'
    return msg

def get_protein_seq(cds_feature : SeqFeature) -> str:
    '''
    Extract protein sequence as a string from a give SeqFeature.  If there is 
    no sequence and the protein is annotated as pseudo, return a empty string.  
    If there is no sequence and no pseudo annotation, raise KeyError.

    Parameters
    ----------
    cds_feature : SeqFeature
        Protein feature from SeqRecord.

    Raises
    ------
    ValueError
        cds_feature type attribute != 'CDS'.
    KeyError
        There is no protein sequence and no pseudo annotation.

    Returns
    -------
    str
        A protein sequence.

    '''
    if cds_feature.type != 'CDS':
        raise ValueError('Not a CDS feature')
    
    #Get the translation associated with the protein
    try:
        protein_seq = ''.join(cds_feature.qualifiers['translation'])
    
    #If there is no translation, check if pseudo (return '') or not (raise KeyError)
    except KeyError:
        if 'pseudo' in cds_feature.qualifiers.keys():
            return ''
        else:
            raise KeyError
    
    #If there is a translation but it is empty, check if pseudo (return '') or 
    #not (raise KeyError)
    if protein_seq == '':
        if 'pseudo' in cds_feature.qualifiers.keys():
            return ''
        else:
            raise KeyError 
    
    #Return non-empty translation
    return protein_seq



def log_params(local_vars : dict, command : str, logger_name : str) -> str:
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
    
    logger_name : str
        name of logger with which to record params.

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
    logger = logging.getLogger(logger_name)
    if command == 'collect':
        log_keys = ['binary_path', 'strict_span', 'neighbourhood_size', 
                    'write_genomes', 'email', 'filenames', 'results_dir']
    elif command == 'sieve':
        log_keys = ['input_genbank_dir', 'e_value', 'min_percent_identity', 
                    'similarity_filter', 'results_dir', 'output_blast_dir',
                    'output_genbank_dir', 'output_vis_dir', 'min_edge_view']

         
    else:
        raise KeyError(f'Command {command} must be collect or sieve')
    log_string = f'Command: {command}\n'
    for log_key in log_keys:
        try:
            log_string += f'{log_key}: {local_vars[log_key]}\n'
        except KeyError:
            raise KeyError(f'{log_key} not in locals - {locals()}')
            
    logger.info(f'---PARAMETERS---\n{log_string}\n\n')

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

def initalise_log(log_file : str, logger_name : str) -> logging.Logger:    
    '''
    Set up logger

    Parameters
    ----------
    log_file : str
        File in which to write log.
    logger_name : str
        Name of logger.

    Returns
    -------
    logger : logging.Logger
        Logger object.

    '''
    log_handlers = get_log_handlers(log_file)
    logging.basicConfig(
                        level=logging.INFO,
                        format='%(asctime)s %(levelname)s:  %(message)s',
                        handlers=log_handlers
                    )
    logger = logging.getLogger(logger_name)
    return logger

def ingest_binary_csv(binary_path : str, logger_name : str) -> pd.DataFrame:
    '''
    Process binary file to pandas dataframe and check data integrity.

    Parameters
    ----------
    binary_path : str
        Path to cblaster binary file.
    logger_name : str
        Name of logger to save results.

    Raises
    ------
    ValueError
        Unexpected CSV format.

    Returns
    -------
    location_data : pd.Dataframe
        Dataframe containing cblaster binary file data.

    '''
    location_data = pd.read_csv(binary_path, sep = ',')
    if len(location_data.columns)<=5:
        logger = logging.getLogger(logger_name)
        error = "Unexpected binary file format\nexpected columns - "\
                    "['Organism', 'Scaffold', 'Start', 'End', 'Score', one or "\
                         "more Query Gene Names].\n"\
                             f"Actual columns: {list(location_data.columns)}"
        logger.error(error)
        raise ValueError(error) 
    return location_data