# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 20:36:03 2024

@author: u03132tk
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature 
import os

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
        cds_count += len([f for f in record.features if f.type == 'CDS'])
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

