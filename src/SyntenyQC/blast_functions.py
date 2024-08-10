# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 13:12:57 2024

@author: u03132tk
"""
from Bio.Blast import NCBIXML
from Bio.Blast.Record import Alignment, HSP
import shutil
import subprocess
from SyntenyQC.helpers import get_cds_count, get_gbk_files, read_gbk, get_protein_seq
import os
from typing import TextIO

'''
This module outlines code for performing reciprocal best hits (used by pipelines.sieve())
'''

def run_blast_process(cmd : list):
    r'''
    Run command via the subprocess module in a blastp-compatible environment.
    For discussion re env, see https://www.biostars.org/p/413294/#415002 and
    https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.makeblastdb_application_opt/

    Parameters
    ----------
    cmd : list
        Command line options (e.g. ['path/to/file.exe', '-flag1', 'paramter1']).

    '''
    try:
        return subprocess.run(cmd, 
                              check=True, 
                              capture_output=True, 
                              env={'BLASTDB_LMDB_MAP_SIZE':'1000000000'})
    except FileNotFoundError as e1:
        print ('Command:  ', cmd)
        print ('FileNotFoundError:\n', e1)
        return e1
    except subprocess.CalledProcessError as e2:
        print ('Command:  ', cmd)
        print ('CalledProcessError:')
        print ('\nError object:\n',e2)
        print ('\nStd_err:\n',e2.stderr)
        print ('\nStd_out:\n',e2.stdout)
        return e2

def makeblastdb_subprocess(makeblastdb_exe_path : str, input_fasta_path : str, 
                           db_out_path : str):
    r'''
    This runs the BLAST+ makeblastdb exe using Python via the subprocess module.  
    Input fasta (path = input_fasta_path) containing database sequences is used 
    to make several BLAST db files with the same name and different filetypes
    (.pdb, .phr, .pin, .pot, .psq, .ptf, .pto). db_out_path is where you want 
    the database files to be placed. NB - db_out_path should have the full path 
    up until where you want the files to be made, including the name shared by 
    all database files - do NOT include a file suffix.
    
    e.g. db_out_path = path\to\dir\name will generate database files 
    path\to\dir\name.pdb, path\to\dir\name.phr, etc...

    Parameters
    ----------
    makeblastdb_exe_path : str
        Path to local NCBI makeblastdb executable.
    input_fasta_path : str
        Path to fasta file with sequences used to make BLASTP database.
    db_out_path : str
        File path used to write database files (do not include a file suffix).
    '''
    return run_blast_process([makeblastdb_exe_path, 
                              '-in', input_fasta_path, 
                              '-out', db_out_path, 
                              '-dbtype', 'prot'])
                          


def blastP_subprocess(evalue_user : float, query_user : str, blastp_exe_path : str, 
                      results_out_path : str, db_user : str, thread_num : int, 
                      max_target_seqs : int):
    r'''
    This runs the BLAST+ blastp exe using Python via the subprocess module.  
    
    Parameters
    ----------
    evalue_user : float
        BLASTP maximum evalue.
    query_user : str
        Path to query fasta file.
    blastp_exe_path : str
        Path to local BLASTP exe file..
    results_out_path : str
        Path used to write BLASTP results XML file.
    db_user : str
        Path to BLASTP database files. Do not include a filetype suffix - if 
        your database files are path\to\dir\name.pdb, path\to\dir\name.phr, etc, 
        set db_user to path\to\dir\name.
    thread_num : int
        BLASTP thread number.
    max_target_seqs : int
        BLASTP max_target_seqs.
    '''
    
    return run_blast_process([blastp_exe_path, 
                              '-out', fr'{results_out_path}', 
                              '-query', fr'{query_user}', 
                              '-db', fr'{db_user}', 
                              '-evalue', fr'{evalue_user}', 
                              '-outfmt', '5', 
                              '-num_threads', fr'{thread_num}',
                              '-max_target_seqs', str(max_target_seqs)
                              ])



class FastaWriter:
    '''
    Class to write all proteins in a collection of genbank files to a single 
    fasta format file.
    '''
    
    def __init__(self, genbank_folder : str, output_filepath : str):
        '''
        Open output_filepath in write mode and write data as described in 
        write_fasta().

        Parameters
        ----------
        genbank_folder : str
            Folder with genbank files.
        output_filepath : str
            Filepath to fasta.

        '''
        with open (output_filepath, 'w') as outfile:
            self.write_fasta(genbank_folder, outfile)
            
    @staticmethod
    def write_fasta(genbank_folder : str, outfile_handle : TextIO) -> TextIO:
        '''
        Extract protein sequences from gbk files in genbank_folder and write to 
        outfile_handle.  Ignore features annotated as pseudo unless they have a 
        sequence.

        Parameters
        ----------
        genbank_folder : str
            Folder with gbk files.
        outfile_handle : TextIO
            Handle to write protine sequences in fasta format.

        Raises
        ------
        KeyError
            There is a non-pseudo protein feature in >=1 gbk file that has no 
            protein sequence.
        ValueError
            There are no CDS features with a translation in >= 1 gbk file(s) in 
            genbank_folder
            
        Returns
        -------
        outfile_handle : TextIO
            The input handle, populated with fasta data.

        '''
        fasta = []
        for file in get_gbk_files(genbank_folder):
            record = read_gbk(f'{genbank_folder}\\{file}')
            
            #Record Id will be used to split BLASTP results by input record
            #Take from filename rather than genbank annotations (e.g. organism 
            #or accession) as these may not be unique
            record_id = file[0:file.rindex('.')]
            #use cds count rather thn index so you only count legitimate proteins
            cds_count = 0
            for feature in record.features:
                if feature.type == 'CDS':
                    protein_seq = get_protein_seq(feature)
                    if protein_seq == '':
                        continue
                    else:
                        fasta += [(f'>{record_id}__{cds_count}', 
                                   protein_seq)
                                  ]
                        cds_count += 1
            if cds_count == 0:
                raise ValueError(f'{file} has no CDS features with a protein sequence')
        #write fasta data to handle
        for index, (defline, seq) in enumerate(fasta):
            assert defline[0] == '>'
            outfile_handle.write(defline + '\n')
            outfile_handle.write(seq)
            if index < len(fasta) - 1:
                outfile_handle.write('\n')
                
        return outfile_handle
    

    

def all_vs_all_blast(folder_with_genbanks : str, e_value : float, 
                     max_target_seqs : int, blast_dir : str) -> str:
    '''
    Run all v all blast, return XML results path.

    Parameters
    ----------
    folder_with_genbanks : str
        Folder with genbanks that will be used to write the BLASTP query fasta.
    e_value : float
        BLASTP evalue threshold.
    max_target_seqs : int
        BLASTP max_target_seqs.
    blast_dir : str
        folder to write blast database and result files
    Returns
    -------
    str
        Path to BLASTP results xml.
    '''
    
    #check your results directory exists
    if not os.path.isdir(blast_dir):
        raise ValueError (f'blast_dir does not exist - {blast_dir}')
        
    #define filepaths
    all_proteins = f'{blast_dir}\\all_proteins.txt' 
    all_proteins_db = f'{blast_dir}\\all_proteins_db'
    results_out_path =  f'{blast_dir}\\results.xml' 
    
    #find blast+ blastp and makeblast db executables
    makeblastdb_exe_path = shutil.which('makeblastdb')
    blastp_exe_path = shutil.which('blastp')#'blast-2.10.1+')

    #write fasta    
    FastaWriter(folder_with_genbanks, 
                output_filepath = all_proteins)   
    
    #run makeblastdb
    print (f'Making database - {all_proteins_db}...')
    makeblastdb_subprocess(makeblastdb_exe_path, 
                           all_proteins, 
                           all_proteins_db)
    cds_count = get_cds_count(genbank_folder = folder_with_genbanks)
    
    #run blastp
    print (f'Running all vs all BLASTP for {len(get_gbk_files(folder_with_genbanks))} '\
               f'neighbourhoods in {folder_with_genbanks} - {cds_count} protein sequences, '\
                   f'{max_target_seqs} max_target_seqs...')  
    blastP_subprocess (e_value, 
                       all_proteins, 
                       blastp_exe_path, 
                       results_out_path, 
                       all_proteins_db, 
                       4,
                       max_target_seqs)
    print (f'Completed BLASTP - results at {results_out_path}...')
    
    return results_out_path

def get_best_hsp(alignment : Alignment, min_percent_identity : int) -> HSP:
    '''
    Get best high scoring pair (hsp) from a supplied BLASTP alignment, ranked 
    according to BLASTP score.

    Parameters
    ----------
    alignment : Alignment
        Biopython blast record alignment.
    min_percent_identity : int
        Minimum identity for a HSP to be considered.

    Returns
    -------
    HSP
        The best biopython blast record hsp.

    '''
    best_hsp = None
    for hsp in alignment.hsps:
        percent_identity = 100*(hsp.identities / hsp.align_length) 
        if percent_identity < min_percent_identity:
            continue
        if best_hsp == None:
            best_hsp = hsp
        else:
            if hsp.score > best_hsp.score:
                best_hsp = hsp
    return best_hsp

def read_xml(results_path : str) -> list:
    '''
    Read and parse BLASTP xml results.

    Parameters
    ----------
    results_path : str
        Path to xml file.

    Returns
    -------
    list
        List of biopython SeqRecords.

    '''
    with open(results_path, 'r') as result_handle:
        blast_records = list(NCBIXML.parse(result_handle))
    return blast_records

def results_to_hits(blast_records : list, min_percent_identity : int) -> dict:
    '''
    Parse xml results to dictionary of hits. 
                

    Parameters
    ----------
    blast_records : str
        XML data parsed to list of biopython SeqRecords.
    min_percent_identity : int
        Min alignment identity.

    Returns
    -------
    raw_results : dict
        Dictionary of hits, of the format:
            
            {query_file_A : 
                {query_protein_1 : 
                    {hit_file_A : 
                        {hit_protein_1A : best HSP (query_protein_1 vs hit_protein_1A),
                         hit_protein_2A : best HSP (query_protein_1 vs hit_protein_2A),
                         ...
                         },
                     hit_file_B : 
                         {hit_protein_1B : best HSP (query_protein_1 vs hit_protein_1B),
                          ...
                         }
                    },
                query_protein_2 : {..................................},
                ...
                },
            query_file_B : {.........................................},
            ...
            }

    '''
    raw_results = {}#
    for record in blast_records:
        query_scaffold, query_index = record.query.split('__')
        if query_scaffold not in raw_results.keys():
            raw_results[query_scaffold] = {}
        if query_index not in raw_results[query_scaffold].keys():
            raw_results[query_scaffold][query_index] = {}
        for hit in record.alignments:
            hit_scaffold, hit_index = hit.hit_def.split('__')
            best_hsp = get_best_hsp(hit, min_percent_identity)
            if hit_scaffold not in raw_results[query_scaffold][query_index].keys():
                raw_results[query_scaffold][query_index][hit_scaffold] = {}
            raw_results[query_scaffold][query_index][hit_scaffold][hit_index] = best_hsp
    return raw_results


def hits_to_best_hits(hit_matrix : dict) -> dict:
    '''
    Parse dictionary of hits to dictionary of best hits.

    Parameters
    ----------
    hit_matrix : dict
        Dictionary of hits, of the format:
            
            {query_file_A : 
                {query_protein_1 : 
                    {hit_file_A : 
                        {hit_protein_1A : best HSP (query_protein_1 vs hit_protein_1A),
                         hit_protein_2A : best HSP (query_protein_1 vs hit_protein_2A),
                         ...
                         },
                     hit_file_B : 
                         {hit_protein_1B : best HSP (query_protein_1 vs hit_protein_1B),
                          ...
                         }
                    },
                query_protein_2 : {..................................},
                ...
                },
            query_file_B : {.........................................},
            ...
            }.

    Returns
    -------
    best_hits : dict
        Dictionary of best hits, of the format:
            
            {query_file_A : 
                {query_protein_1 : 
                    {hit_file_A : id_of_best_query_protein_1_hit_in_hit_file_A,
                     hit_file_B : id_of_best_query_protein_1_hit_in_hit_file_B,
                     ...
                     },
                query_protein_2 : 
                    {hit_file_A : id_of_best_query_protein_2_hit_in_hit_file_A,
                     ...
                     },
                ...
                }
            query_file_B : 
                {.........................................
                },
            ..............................................
            }

    '''
    best_hits = {}
    for query_scaffold, query_proteins in hit_matrix.items():
        best_hits[query_scaffold] = {}
        for query_index, hit_scaffolds in query_proteins.items():
            best_hits[query_scaffold][query_index] = {}
            for hit_scaffold, hit_proteins in hit_scaffolds.items():
                best_hit = None
                best_score = 0
                for hit_index, hit_hsp in hit_proteins.items():
                    if hit_hsp == None:
                        continue
                    if hit_hsp.score > best_score:
                        best_hit = hit_index
                        best_score= hit_hsp.score
                if best_hit == None:
                    continue
                best_hits[query_scaffold][query_index][hit_scaffold] = best_hit
    return best_hits


def best_hits_to_rbh(best_hit_matrix : dict) -> dict:
    '''
    Parse dictionary of best hits to dictionary of reciprocal best hits

    Parameters
    ----------
    best_hit_matrix : dict
        Dictionary of hits, of the format:
            
            {query_file_A : 
                {query_protein_1 : 
                    {hit_file_A : id_of_best_query_protein_1_hit_in_hit_file_A,
                     hit_file_B : id_of_best_query_protein_1_hit_in_hit_file_B,
                     ...
                     },
                query_protein_2 : 
                    {hit_file_A : id_of_best_query_protein_2_hit_in_hit_file_A,
                     ...
                     },
                ...
                }
            query_file_B : 
                {.........................................
                },
            ..............................................
            }.

    Returns
    -------
    reciprocal_best_hits : dict
        Dictionary of reciprocal best hits, of the format (numbers are CDS):
            
            {'file1' : {'0' : {'file1' : '0',
                               'file2' : '0'},
                        '1' : {'file1' : '1'},
                        },

             'file2' : {'0' : {'file1' : '0',
                               'file2' : '0'},
                        '1' : {'file2' : '1'}
                        }
             }

    '''
    reciprocal_best_hits = {}
    for query_scaffold, query_proteins in best_hit_matrix.items():
        for query_index, best_hit_proteins in query_proteins.items():
            for hit_scaffold, best_hit_protein in best_hit_proteins.items():
                if query_scaffold not in best_hit_matrix[hit_scaffold][best_hit_protein].keys():
                    continue # one direction hit
                if best_hit_matrix[hit_scaffold][best_hit_protein][query_scaffold] == query_index:
                    if query_scaffold not in reciprocal_best_hits.keys():
                        reciprocal_best_hits[query_scaffold] = {}
                    if query_index not in reciprocal_best_hits[query_scaffold].keys():
                        reciprocal_best_hits[query_scaffold][query_index] = {}
                    reciprocal_best_hits[query_scaffold][query_index][hit_scaffold] = best_hit_protein
    return reciprocal_best_hits
    
def make_rbh_matrix(xml_path : str, min_percent_identity : int) -> dict:
    '''
    Convert xml file to reciprocal best hits dictionary.

    Parameters
    ----------
    xml_path : str
        Path to xml results.
    min_percent_identity : int
        Minimum BLASTP alignment percentage.

    Returns
    -------
    rbh_matrix : duct
        RBH matrix represented as a dictionary.

    '''
    blast_records = read_xml(xml_path)
    hit_matrix = results_to_hits(blast_records, 
                                 min_percent_identity)
    print ('Processed hits...')
    best_hit_matrix = hits_to_best_hits(hit_matrix)
    print ('Processed best hits')
    rbh_matrix = best_hits_to_rbh(best_hit_matrix)
    print ('Processed reciprocal best hits')
    return rbh_matrix