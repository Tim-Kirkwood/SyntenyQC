# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:45:42 2024

@author: u03132tk
"""
import argparse
from SyntenyQC.pipelines import collect, sieve
from SyntenyQC.helpers import get_gbk_files
import os

'''
This module contains command line parser functions to convert sys.argv commands 
into parameters supplied to either pipelines.collect() or pielines.sieve().
'''


def mixed_slashes(path : str) -> bool:
    '''
    Check if an input string (assumed to be a path) has both forward (/) 
    and backward (\\) slashes.
    
    Note - os.isdir and os.isfile dont care about mixed slashes, but if a user 
    inputs a mixed-slash path then that suggests the user may be doing something 
    wrong (so this function forces them to be careful). 


    Parameters
    ----------
    path : str
        Path to check.

    Raises
    ------
    ValueError
        No slashes are present - file/dir paths must have >= 1.

    Returns
    -------
    bool
        Mixed slashes are/are not present.

    '''
    is_forward_slash =  path.count('/') > 0
    is_backward_slash = path.count('\\') > 0
    if not (is_forward_slash or is_backward_slash):
        raise ValueError('Not a recognisable path (no slashes)')
    return is_forward_slash and is_backward_slash

def make_dirname(parent_dir : str, base_dir_name : str) -> str:
    '''
    Make a path to a new dir (base_dir_name) within parent_dir of the format 
    parent_dir/base_dir_name.  If a dir called base_dir_name is already present 
    in parent_dir, add an increment to the end of the path 
    - i.e. parent_dir/base_dir_name(increment).  

    Parameters
    ----------
    parent_dir : str
        Dir that will hold new base_dir_name dir.
    base_dir_name : str
        Name of new subdirectory in parent_dir.

    Raises
    ------
    ValueError
        parent_dir does not exist.

    Returns
    -------
    str
        Path of format parent_dir/base_dir_name or 
        parent_dir/base_dir_name(increment).
    '''
    if not os.path.isdir(parent_dir):
        raise ValueError(f'{parent_dir} is not a dir')
    temp_dir = base_dir_name
    count = 1
    while temp_dir in os.listdir(parent_dir):
        temp_dir = f'{base_dir_name}({count})'
        count += 1
    return f'{parent_dir}\\{temp_dir}'


def check_file_path_errors(long_var_name : str, path : str, suffixes : list) -> tuple:
    '''
    Return an error code and error message if path has mixed slashes, leads to 
    a file that does not exist, or has an unrecognised file type. 
    
    Main parser error codes:
        1 = incompatible inputs (e.g. input outside of allowed range)
        2 = non-existent files/folders

    Parameters
    ----------
    long_var_name : str
        Command line parameter associated with path (binary_path in this context).
    path : str
        Filepath to check.
    suffixes : list
        Accepted file types (.csv or .txt in csv format).

    Returns
    -------
    tuple
        (error code, message).

    '''
    if mixed_slashes(path):
        return (1, 
                f'--{long_var_name} path cannot contain forward and backward slashes.')
    if not os.path.isfile(path):
        return (2, 
               f'--{long_var_name} file does not exist.')
    suffix = path[path.rindex('.'):]
    if suffix not in suffixes:
        return (1, 
                f'--{long_var_name} must be in csv format (as a file with suffix {suffixes}).')
    return (None, None)

def check_dir_path_errors(long_var_name : str, path : str) -> tuple:
    '''
    Return an error code and error message if path has mixed slashes or does 
    not lead to an existing dir. 
    
    Main parser error codes:
        1 = incompatible inputs (e.g. input outside of allowed range)
        2 = non-existent files/folders

    Parameters
    ----------
    long_var_name : str
        Command line parameter associated with path (genbank_folder in this context).
    path : str
        Dirpath to check.

    Returns
    -------
    tuple
        (error code, message).

    '''
    if mixed_slashes(path):
        return (1, 
                f'--{long_var_name} path cannot contain forward and backward slashes.')
    if not os.path.isdir(path):
        return (2, 
               f'--{long_var_name} dir does not exist.')
    return (None, None)


def read_args(arg_list: list[str] | None = None):
    '''
    Read args supplied via arg_list or sys.argv.  Argparse will raise SystemExit if there is 
    an arg issue.

    Parameters
    ----------
    arg_list : list[str] | None, optional
        Arguments for argparse. The default is None.
        
    Raises
    ------
    SystemExit
        Argparse does not recognise arg in arg_list.

    Returns
    -------
    tuple
        (argument namespace, argparse.ArgumentParser).

    '''
    #https://pythontest.com/testing-argparse-apps/
    global_parser = argparse.ArgumentParser(prog="SyntenyQC")
    subparsers = global_parser.add_subparsers(title="subcommands", 
                                              help="Synteny quality control options",
                                              dest='command' #https://stackoverflow.com/a/9286586/11357695
                                              )
    collect_parser = subparsers.add_parser("collect", 
                                           description='''Write genbank files corresponding to cblaster neighbourhoods from a specified CSV-format binary file loacted at BINARY_PATH.  For each cblaster hit accession in the binary file:

1) A record is downloaded from NCBI using the accession.  NCBI requires a user EMAIL to search for this record programatically.  If WRITE_GENOMES is specified, this record is written to a local file according to FILENAMES (see final bulletpoint).
2) A neighbourhood of size NEIGHBOURHOOD_SIZE bp is defined, centered on the cblaster hits defined in the binary file for the target accession. 
3) (If STRICT_SPAN is specified:) If the accession's record is too small to contain a neighbourhood of the desired size, it is discarded.  For example, if an accession record is a 25kb contig and NEIGHBOURHOOD_SIZE is 50000, the record is discarded.
4) If FILENAMES is "organism", the nighbourhood is written to file called *organism*.gbk. If FILENAMES is "accession", the neighbourhood is written to *accession*.gbk. Synteny softwares such as clinker can use these filesnames to label synetny plot samples.
                                            
Once COLLECT has been run, a new folder with the same name as the binary file should be created in the directory that holds the binary file (i.e. the file "path/to/binary/file.txt" will generate the folder "path/to/binary/file"). This folder will have a subdirectory called "neighbourhood", containing all of the neighbourhood genbank files (i.e. "path/to/binary/file/neighbourhood"). If WRITE_GENOMES is specified, a second direcory ("genome") will also be present, containing the entire record associated with each cblaster accession (i.e. "path/to/binary/file/genome").  Finally, a log file will be present in the folder "path/to/binary/file", containing a summary of accessions whose neighbourhoods were discarded.''', 
                                           formatter_class=argparse.RawDescriptionHelpFormatter)
                                                
    collect_parser.add_argument("-bp", 
                                "--binary_path", 
                                type = str,
                                required=True,
                                metavar='\b',
                                help = '''Full filepath to the CSV-format cblaster binary file containing neighbourhoods 
                                that should be extracted''')
    collect_parser.add_argument("-ns", 
                                "--neighbourhood_size", 
                                type = int,
                                required=True,
                                metavar='\b',
                                help = '''Size (basepairs) of neighbourhood to be extracted (centered on middle of CBLASTER-defined 
                                neighbourhood)''')
    collect_parser.add_argument("-em", 
                                "--email", 
                                type = str,
                                required=True,
                                metavar='\b',
                                help = 'Email - required for NCBI entrez querying')
    collect_parser.add_argument("-fn", 
                                "--filenames", 
                                type = str,
                                choices=['organism', 'accession'],
                                default='organism',
                                metavar='\b',
                                help = '''If "organism", all collected files will be named according to organism.  
                                If "accession", all files will be named by NCBI accession. (default: %(default)s)''')
    collect_parser.add_argument("-sp", 
                                "--strict_span", 
                                action = 'store_true',
                                help = '''If set, will discard all neighbourhoods that are smaller than neighbourhood_size bp.  
                                For example, if you set a neighbourhood_size of 50000, a 50kb neighbourhood will be extracted 
                                from the NCBI record associateed with each cblaster hit.  If the record is too small for this 
                                to be done (i.e. the record is smaller then 50kb) it is discarded''') 
    collect_parser.add_argument("-wg", 
                                "--write_genomes", 
                                action = 'store_true',
                                help = '''If set, will write entire NCBI record containing a cblaster hit to file 
                                (as well as just the neighbourhood)''')

    sieve_parser = subparsers.add_parser("sieve", 
                                         description='''
Filter redundant genomic neighbourhoods based on neighbourhood similarity:
- First, an all-vs-all BLASTP is performed with user-specified BLASTP settings and the neighbourhoods in GENBANK_FOLDER.  
- Secondly, these are parsed to define reciprocal best hits between every pair of neighbourhoods.  
- Thirdly, these reciprocal best hits are used to derive a neighbourhood similarity network, where edges indicate two neighbourhood nodes that have a similarity > SIMILARITY_FILTER. Similarity = Number of RBHs / Number of proteins in smallest neighbourhood in pair.   
- Finally, this network is pruned to remove neighbourhoods that exceed the user's SIMILARITY_FILTER threshold. Nodes that remain are copied to the newly created folder "genbank_folder/sieve_results/genbank".
''',
                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    sieve_parser.add_argument("-gf", 
                              "--genbank_folder", 
                              type = str,
                              required=True,
                              metavar='\b',
                              help = 'Full path to folder containing neighbourhood genbank files requiring de-duplication')
    
    sieve_parser.add_argument("-ev", 
                              "--e_value", 
                              type = float, 
                              default = 10**-5,
                              metavar='\b',
                              help = 'BLASTP evalue threshold. (default: %(default)s)')
    sieve_parser.add_argument("-mi", 
                              "--min_percent_identity", 
                              type = int, 
                              default = 50,
                              metavar='\b',
                              help = 'BLASTP percent identity threshold. (default: %(default)s)')
    sieve_parser.add_argument("-mts", 
                              "--max_target_seqs", 
                              type = int, 
                              default = 200,
                              metavar='\b',
                              help = 'BLASTP -max_target_seqs. Maximum number of aligned sequences to keep. (default: %(default)s)')
    sieve_parser.add_argument("-mev", 
                              "--min_edge_view", 
                              type = float, 
                              #default = 0.5,
                              metavar='\b',
                              help = 'Minimum similarity between two neighbourhoods for an edge to be drawn betweeen them in the RBH graph.  Purely for visualisation of the graph HTML file - has no impact on the graph pruning results. (default: %(default)s)')
    
    sieve_parser.add_argument("-sf", 
                              "--similarity_filter", 
                              type = float,
                              required=True,
                              metavar='\b',
                              help = 'Similarity threshold above which two neighbourhoods are considered redundant')
    args = global_parser.parse_args(arg_list)
    return args, global_parser


        

def check_args(args, parser):
    '''
    Sanitise args and comapre args namespace from read_args to check whether each 
    argument is acceptable (i.e. points to existing files/folders and is within 
    an acceptable range). Call parser.exit() if a param fails a check (this in 
    turn raises a SystemExit).
    
    Error codes:
        1 = incompatible inputs (e.g. input outside of allowed range)
        2 = non-existent files/folders      
    
    Raises
    ------
    SystemExit
        A param fails it's check(s).
    '''
        
    if args.command == 'collect':
        
        #BINARY PATH - check and sanitise
        error_code, message = check_file_path_errors('binary_path', 
                                                     args.binary_path, 
                                                     ['.txt', '.csv'])
        if error_code != None:
            parser.exit(error_code, message)
        else:
            args.binary_path = args.binary_path.replace('/', '\\')
        
        #NEIGHBOURHOOD SIZE    
        if args.neighbourhood_size <=0:
            parser.exit(1, 
                        '--neighbourhood_size must be greater than 0')
        
        #EMAIL
        if '@' not in args.email:
            parser.exit(1, 
                        '--email is not correct - no @')
        
            
    elif args.command == 'sieve':
        
        #GENBANK FOLDER 
        #check/sanitise path
        error_code, message = check_dir_path_errors('genbank_folder', 
                                                    args.genbank_folder)
        if error_code != None:
            parser.exit(error_code, message)
        else:
            args.genbank_folder = args.genbank_folder.replace('/', '\\')
        #check genbank folder has >= 1 gbk file
        if get_gbk_files(args.genbank_folder) == []:
            parser.exit(2, 
                        f'No genbank (.gbk, .gb) files in {args.genbank_folder}')
        
        #SIMILARITY FILTER
        if not 0 < args.similarity_filter <= 1:
            parser.exit(1, 
                        '--similarity_filter must be between >0 and <=1.')
        
        #MIN EDGE VIEW
        if args.min_edge_view != None:
            if not 0 < args.min_edge_view <= 1:
                parser.exit(1, 
                            '--min_edge_view must be between >0 and <=1.')
            if not args.min_edge_view <= args.similarity_filter:
                parser.exit(1, 
                            '--min_edge_view must be <= similarity_filter.')
        else:
            args.min_edge_view = args.similarity_filter
        
        #BLASTP EVALUE
        if not 0 < args.e_value <= 1:
            parser.exit(1, 
                        '--e_value must be between >0 and <=1.')
        
        #BLASTP PERCENT IDENTITY
        if not 0 < args.min_percent_identity <= 100:
            parser.exit(1, 
                        '--min_percent_identity must be between >0 and <=100.')
        
        #BLASTP MAX TARGET SEQS
        if args.max_target_seqs <= 0:
            parser.exit(1, 
                        '--max_target_seqs must be between >0.')
        
    else:
        #DEFENSIVE PROGRAMMING - I dont think this logic block can be entered,
        #as argparse would fail without a recognised command, but this will 
        #remind me to update this function (and others) if I add more commands 
        parser.exit(1, f'Subcommand not recognised - {args.command}')



    
    
    
def main_cli(arg_list: list[str] | None = None):
    '''
    Parse and check args from the command line.  Perform some basic path construction 
    to get result folders.  Call selected pipeline with args from args namespace 
    object.

    Parameters
    ----------
    arg_list : list[str] | None, optional
        Arguments for argparse. The default is None.
        
    Raises
    ------
    SystemExit
        Successful pipeline completeion (error code 0) or command not recognised 
        (error code 1).

    '''
    args, parser = read_args(arg_list)
    check_args(args, parser) 
    if args.command == 'collect':
        
        #if binary_path = a/folder/binary.csv, make a new dir a/folder/binary 
        #to contain results 
        binary_folder = args.binary_path[0 : args.binary_path.rindex("\\")]
        binary_file = args.binary_path[args.binary_path.rindex("\\") +1 : args.binary_path.rindex(".")]
        results_dir = make_dirname(binary_folder, binary_file)
        os.makedirs(results_dir)
        
        #run pipelines.collect()
        results_path = collect(binary_path = args.binary_path, 
                               strict_span = args.strict_span, 
                               neighbourhood_size = args.neighbourhood_size, 
                               write_genomes = args.write_genomes, 
                               email = args.email, 
                               filenames = args.filenames,
                               results_dir = results_dir)
        
        #exit upon successuful completion 
        parser.exit(status = 0, 
                    message = f'Succsfully ran SyntenyQC -collect.  Results in {results_path}.')
        
        
    elif args.command == 'sieve':
        
        #if genbank_folder = a/folder, make new dir a/folder/sieve_results
        results_dir = make_dirname(args.genbank_folder,'sieve_results')
        
        #run pipelines.sieve()
        results_path = sieve(input_genbank_dir = args.genbank_folder, 
                             e_value = args.e_value, 
                             min_percent_identity = args.min_percent_identity, 
                             max_target_seqs = args.max_target_seqs,
                             similarity_filter = args.similarity_filter,
                             results_dir = results_dir,
                             min_edge_view = args.min_edge_view)
        
        #exit upon successuful completion 
        parser.exit(status = 0, 
                    message = f'Succsfully ran SyntenyQC -sieve.  Results in {results_path}.')
        
    else:
        
        #DEFENSIVE PROGRAMMING - I dont think this logic block can be entered,
        #as argparse would fail without a recognised command, but this will 
        #remind me to update this function (and others) if I add more commands 
        parser.exit(1, f'Subcommand not recognised - {args.command}')
    
