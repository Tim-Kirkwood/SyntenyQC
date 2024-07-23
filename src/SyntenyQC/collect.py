# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 09:54:02 2024

@author: u03132tk
"""
import pandas as pd
import logging
import time 
from SyntenyQC.neighbourhood import Neighbourhood, write_results
from Bio import Entrez
from SyntenyQC.helpers import check_motif



def collect(binary_path, 
            strict_span,# 
            neighbourhood_size, 
            write_genomes, #
            email,
            filenames,
            results_dir):
    start = time.time()
    location_data = pd.read_csv(binary_path, sep = ',')
    assert len(location_data.columns)>=5
    Entrez.email = email
    #set up logging for issues with neighbourhood scraping and termini issues
    
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
    
    number_of_neighbourhoods = len(location_data['Scaffold'])
    logger.info(f"Extracting {number_of_neighbourhoods} gbks...")
    neighbourhood_data = zip (location_data['Scaffold'],
                             location_data['Start'],
                             location_data['End'])
    writtten = 0
    skip_termini = 0
    skip_scrape = 0
    skip_motif = 0
    for index, (accession, motif_start, motif_stop) in enumerate(neighbourhood_data):
        
        logger.info(f'accession {accession} #{index} of {number_of_neighbourhoods - 1}')  
        error_msg = check_motif(motif_start, motif_stop)
        if error_msg != '':
            logger.warn(error_msg) 
            raise ValueError (error_msg)
        neighbourhood = Neighbourhood(accession, 
                                      motif_start, 
                                      motif_stop, 
                                      neighbourhood_size, 
                                      strict_span)
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
        #TODO file handling
        write_results(results_folder=results_dir, 
                      neighbourhood=neighbourhood, 
                      filenames=filenames, 
                      scale = 'neighbourhood')
        if write_genomes:
            write_results(results_folder=results_dir, 
                          neighbourhood=neighbourhood, 
                          filenames=filenames, 
                          scale = 'genome')
        writtten += 1
    logger.info(f'Written {writtten} records to {results_dir} in {round(time.time() - start)} seconds - {skip_scrape} failed scraping, {skip_termini} had termini that exceeded the genome boundaries, {skip_motif} had a motif that was too long')
    return results_dir



