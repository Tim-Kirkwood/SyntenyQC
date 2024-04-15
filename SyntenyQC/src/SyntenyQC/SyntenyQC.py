# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:45:42 2024

@author: u03132tk
"""
import argparse
from collect import collect
from sieve import sieve
import os

#TODO lowercase
#add help messages and descriptions 

def read_args():
    global_parser = argparse.ArgumentParser(prog="ClusterSieve")
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

    #add these together as one tool - have the flags below as compulsory if -from_cblaster is true.
    #specify either -binary_path or -genbank_path to make sure the user is aware of what theyre doing
    #can check if set as the args are default None
    sieve_parser = subparsers.add_parser("sieve", 
                                         description='''
Filter redundant genomic neighbourhoods based on neighbourhood similarity:
- First, an all-vs-all BLASTP is performed with user-specified BLASTP settings and the neighbourhoods in GENBANK_FOLDER.  
- Secondly, these are parsed to define reciprocal best hits between every pair of neighbourhoods.  
- Thirdly, these reciprocal best hits are used to derive a neighbourhood similarity network, where edges indicate two neighbourhood nodes that have a similarity > SIMILARITY_FILTER. Similarity = (Number of RBHs - CONSERVED_NEIGHBOURHOOD_LEN) / (Number of proteins encoded in smallest neighbourhood - CONSERVED_NEIGHBOURHOOD_LEN).   
- Finally, this network is pruned to remove neighbourhoods that exceed the user's SIMILARITY_FILTER threshold
''',
                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    sieve_parser.add_argument("-gf", 
                              "--genbank_folder", 
                              type = str,
                              required=True,
                              metavar='\b',
                              help = 'Folder containing neighbourhood genbank files requiring de-duplication')
    
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
    sieve_parser.add_argument("-cnl", 
                              "--conserved_neighbourhood_len", 
                              type = int, 
                              default = 0,
                              metavar='\b',
                              help = 'Number of RBHs expected in every neighbourhood (e.g. if you have run a cblaster search requiring a minimum number of hits). (default: %(default)s)')
    sieve_parser.add_argument("-sf", 
                              "--similarity_filter", 
                              type = float,
                              required=True,
                              metavar='\b',
                              help = 'Similarity threshold above which two neighbourhoods are considered redundant')
                     
    args = global_parser.parse_args()
    return args, global_parser

def mixed_slashes(path):
    is_forward_slash =  path.count('/') > 0
    is_backward_slash = path.count('\\') > 0
    return is_forward_slash and is_backward_slash
        

def check_args(args, parser):
    '''
    Error codes:
        1 = incompatible inputs (e.g. input outside of allowed range)
        2 = non-existent files/folders
        3 = results folder(s) already exists 
    TODO:    
        - add blastp options 
        - see what cblaster results look like when not every gene is required
        
    '''
    if args.command == 'collect':
        if mixed_slashes(args.binary_path):
            parser.exit(1, 
                        '--binary_path cannot contain forward and backward slashes.')
        args.binary_path = args.binary_path.replace('/', '\\')
        if not os.path.isfile(args.binary_path):
            parser.exit(2, 
                        '--binary_path file does not exist.')
        binary_file = args.binary_path[args.binary_path.rindex('\\'):args.binary_path.rindex('.')]
        binary_folder = args.binary_path[0 : args.binary_path.rindex("\\")]
        if binary_file in os.listdir(binary_folder):
            parser.exit(3, 
                        f'{args.binary_folder} already contains the folder {binary_file} - please move or rename this folder to ensure file collections are kept seperate.')
            
        if args.neighbourhood_size <=0:
            parser.exit(1, 
                        '--neighbourhood_size must be greater than 0')
            
        if '@' not in args.email:
            parser.exit(1, 
                        '--email is not correct - no @')
        
            
    else:
        if mixed_slashes(args.genbank_folder):
            parser.exit(1, 
                        '--genbank_folder cannot contain forward and backward slashes.')
        args.genbank_folder = args.genbank_folder.replace('/', '\\')
        if not os.path.isdir(args.genbank_folder):
            parser.exit(2, 
                        '--genbank_folder directory does not exist.')
        if 'ClusterSieve' in os.listdir(args.genbank_folder):
            parser.exit(3, 
                        f'{args.genbank_folder} already contains the folder "ClusterSieve" - please move or rename this folder to ensure file collections are kept seperate.')
        if not 0 < args.similarity_filter <= 1:
            parser.exit(1, 
                        '--similarity_filter must be between >0 and <=1.')
        if not 0 < args.e_value <= 1:
            parser.exit(1, 
                        '--e_value must be between >0 and <=1.')
        if not 0 < args.min_percent_identity <= 100:
            parser.exit(1, 
                        '--min_percent_identity must be between >0 and <=100.')
        if not args.conserved_neighbourhood_len >= 0:
            parser.exit(1, 
                        '--conserved_neighbourhood_len must be >=0.')


def main_cli():
    #shift exceptions to parser.exit
    #dont use subcommands - should be single flow either from cblaster or a file for ease of use
    
    
    args, parser = read_args()
    check_args(args, parser) 
    #TODO check you have all expected parameter combinations etc
    if args.command == 'collect':
        collect(binary_path = args.binary_path, 
                strict_span = args.strict_span, 
                neighbourhood_size = args.neighbourhood_size, 
                write_genomes = args.write_genomes, 
                email = args.email, 
                filenames = args.filenames)
        parser.exit(status = 0, 
                    message = 'Succsfully ran ClusterSieve -collect')
    elif args.command == 'sieve':
        sieve(genbank_folder = args.genbank_folder, 
              e_value = args.e_value, 
              min_percent_identity = args.min_percent_identity, 
              conserved_neighbourhood_len = args.conserved_neighbourhood_len, 
              similarity_filter = args.similarity_filter)
        parser.exit(status = 0, 
                    message = 'Succsfully ran ClusterSieve -sieve')
        
    
if __name__ == '__main__':
    #think you need to make these path objects
    #ClusterCollect(binary_path = 'C:/Users/u03132tk/.spyder-py3/E121N603G591R92_binary.txt', 
                  # strict_span = True,# 
                  # neighbourhood_size = 50000, 
                  # write_genomes = True,
                  # email = 'tdjkirkwood@hotmail.com',
                  # filenames = 'organism'
                  # )
    main_cli()