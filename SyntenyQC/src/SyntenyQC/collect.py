# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 09:54:02 2024

@author: u03132tk
"""
import pandas as pd
from Bio import Entrez, SeqIO
from http.client import IncompleteRead
import logging
import os
from urllib.request import HTTPError

class Neighbourhood:        
    def scrape_genome(self, number_of_attempts, accession):
        for attempt in range(number_of_attempts):
            print ('scraping')
            try:
                handle = Entrez.efetch(db = "nucleotide", 
                                        id = accession, 
                                        rettype = "gbwithparts ", 
                                        retmode = "text")  
            except HTTPError:
                print ('HTTPError - trying again')
                continue
            print ('reading')
            try:
                #self.accession = accession
                self.genome = SeqIO.read(handle, 
                                         "genbank")
            except IncompleteRead:
                print ('IncompleteRead - trying again')
                continue
            found = True
            break
        if not found:
            raise ValueError(f'{accession} efetch failed')
        if self.genome.id != accession:
            raise ValueError(f'''
                             record_id {self.genome.id} is not the same as accession {accession}.  
                             This error is unexpected - please raise an issue on the project github.
                             '''
                             )
    def define_neighbourhood(self, motif_start, motif_stop, neighbourhood_size, genome_length):
        if motif_start >= motif_stop:
            raise ValueError(f'Start ({motif_start}) is >= stop ({motif_stop})')
        #self.motif_start = motif_start
        #self.motif_stop = motif_stop
        motif_span = motif_stop-motif_start
        if motif_span > neighbourhood_size:
            raise ValueError(rf'{motif_span} larger than neighbourhood size {neighbourhood_size}')
        extension = (neighbourhood_size - motif_span)/2
        self.neighbourhood_start = motif_start - extension#extension
        self.neighbourhood_stop = motif_stop + extension#extension
        
    def check_neighbourhood(self, neighbourhood_start, neighbourhood_stop, genome_length, strict_span):
        if neighbourhood_stop > genome_length:
            if strict_span:
                raise ValueError('neighbourhood extends beyond genome terminus')
            else:
                self.neighbourhood_stop = genome_length
        if neighbourhood_start < 0:
            if strict_span:
                raise ValueError('neighbourhood extends beyond genome terminus')
            else:
                self.neighbourhood_start = 0
    def get_neighbourhood(self, start, stop, record):
        print (f'neighbourhood = {start} -> {stop}')
        cut_record = record[int(start) : int(stop)]
        cut_record.annotations = record.annotations
        self.neighbourhood = cut_record
    def write_record(self, record, path):
        with open(path, 'w') as handle:
            SeqIO.write (record, 
                         handle, 
                         "genbank")
    def build(self, 
              accession, 
              motif_start, 
              motif_stop, 
              neighbourhood_size, 
              strict_span):
        #keep error statements for downstream filtering of built neighbourhoods
        self.scrape_fail = False
        self.overlapping_termini = False
        try:
            self.scrape_genome(number_of_attempts = 5,
                               accession=accession)
        except ValueError:
            #TODO build CustomErrors - could theoretically get a value error that is unrelated to function raise
            self.scrape_fail = True
        if not self.scrape_fail:
            self.define_neighbourhood(motif_start= motif_start, 
                                     motif_stop=motif_stop,
                                     neighbourhood_size=neighbourhood_size,
                                     genome_length=len(self.genome)
                                     )
            try:
                self.check_neighbourhood(neighbourhood_start = self.neighbourhood_start,
                                        neighbourhood_stop = self.neighbourhood_stop, 
                                        genome_length = len(self.genome), 
                                        strict_span = strict_span
                                        )
            except ValueError:
                self.overlapping_termini = True
            if not self.overlapping_termini:
                self.get_neighbourhood(start = self.neighbourhood_start,
                                       stop = self.neighbourhood_stop,
                                       record=self.genome)
                

def write_results (results_folder, neighbourhood, filenames, scale):
    folder = f'{results_folder}\\{scale}'
    os.makedirs(folder, 
                exist_ok = True)
    file = make_filename(name_type = filenames, 
                         type_map = {'accession' : neighbourhood.genome.id, 
                                     'organism' : neighbourhood.genome.annotations["organism"]},
                         folder = folder)
    if scale == 'neighbourhood':
        neighbourhood.write_record(record = neighbourhood.neighbourhood,
                                  path = f'{folder}\\{file}')
    elif scale == 'genome':
        neighbourhood.write_record(record = neighbourhood.genome,
                                      path = f'{folder}\\{file}')
    else:
        raise ValueError(f"scale' must be 'neighbourhood' or 'genome', not {scale}")
def make_filename(name_type, type_map, folder):
    file = type_map[name_type]
    count = 0
    if f'{file}.gbk' in os.listdir(folder):
        file += f'_{count}'
        count += 1
    file += '.gbk'
    return file

def collect(binary_path, 
             strict_span,# 
             neighbourhood_size, 
             write_genomes, #
             email,
             filenames):
    
    suffix = binary_path[binary_path.rindex('.'):]
    assert suffix in ['.txt', '.csv'], suffix
    location_data = pd.read_csv(binary_path, sep = ',')
    assert len(location_data.columns)>=5
    Entrez.email = email
    
    #check binary file is a file in csv or txt format, and standardise slashes
    binary_folder = binary_path[0 : binary_path.rindex("\\")]
    binary_file = binary_path[binary_path.rindex("\\") +1 : binary_path.rindex(".")]
    results_dir = f'{binary_folder}\\{binary_file}'
    os.makedirs(results_dir)
    #set up logging for issues with neighbourhood scraping and termini issues
    log_file = f'{results_dir}\\log.txt'
    logging.basicConfig(filename=log_file, 
                        filemode='w')
    logger = logging.getLogger()

    
    number_of_neighbourhoods = len(location_data['Scaffold'])
    print (f"Extracting {number_of_neighbourhoods} gbks...")
    neighbourhood_data = zip (location_data['Scaffold'],
                             location_data['Start'],
                             location_data['End'])
    for index, (accession, motif_start, motif_stop) in enumerate(neighbourhood_data):
        print (f'accession {accession} #{index} of {number_of_neighbourhoods - 1}')  
        neighbourhood = Neighbourhood()
        neighbourhood.build(accession, 
                           motif_start, 
                           motif_stop, 
                           neighbourhood_size, 
                           strict_span)
        if neighbourhood.scrape_fail:
            logger.warning(f'scrape_fail - {accession}')
            continue
        elif neighbourhood.overlapping_termini:
            logger.warning(f'overlapping_termini - {accession}')
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
        
    return results_dir



