# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 12:39:27 2024

@author: u03132tk
"""
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from http.client import IncompleteRead
import os
from urllib.request import HTTPError
import logging 

'''
this module contains code for defining and writing neighbourhoods, used by 
pipelines.collect()
'''

#TODO make unique exceptions for specified fail conditons 
#TODO instead of specific error attr (e.g. scrape_error, motif_to_long...), have 
#     an error attr as list that you append error strings too (will cut down the 
#     code)

class Neighbourhood:       
    '''
    A neighbourhood centered on a motif identified in a query NCBI record by 
    cblaster
    '''
    def __init__(self, 
                 accession : str, 
                 motif_start : int, 
                 motif_stop : int, 
                 neighbourhood_size : int, 
                 strict_span : bool) -> None:
        '''
        Initialise a Neighbourhood object.

        Parameters
        ----------
        accession : str
            Accession used to query NCBI with entrez.
        motif_start : int
            Start of first cblaster hit in accession record.
        motif_stop : int
            End of last cblaster hit in accession record.
        neighbourhood_size : int
            Length of neighbourhood to extract from accession record.
        strict_span : bool
            If True, do not progress Neighbourhood if neighbourhood termini 
            exceed termini of record accession. Otherwise, set start and end
            of Neighbourhood to 0 (if start <0) and record length 
            (if stop > record length) respectively.
        '''
        
        #initialise with acession and no errors    
        self.accession = accession
        self.scrape_error = None
        self.overlapping_termini = False
        self.motif_to_long = False
        
        #try to download the genome - if there is a download error, record and 
        #stop building Neighbourhood 
        try:
            self.genome = self.scrape_genome(number_of_attempts = 5,
                                             accession = accession)
        except HTTPError:
            self.scrape_error = 'HTTP error'
        except IncompleteRead:
            self.scrape_error = 'IncompleteRead'
        except ValueError:
            self.scrape_error = 'ValueError'
        
        if  self.scrape_error == None:
            
            #Define neighbourhood.  If motif length (motif_start - motif_end) 
            #is greater than neighbourhood_size, record and stop building 
            #Neighbourhood 
            try:
                raw_start, raw_stop = self.define_neighbourhood(motif_start, 
                                                                motif_stop,
                                                                neighbourhood_size)
            except ValueError:
                self.motif_to_long = True
            if not self.motif_to_long:
                
                #Sanitise the start and stop values from define_neighbourhood()
                #If strict_span is True, do not progress Neighbourhood if 
                #neighbourhood termini exceed termini of record accession. 
                #Otherwise, set start and end of Neighbourhood to 0 (if start <0) 
                #and record length (if stop > record length) respectively
                genome_length = len(self.genome.seq)
                try:
                    self.neighbourhood_start, \
                        self.neighbourhood_stop = self.sanitise_neighbourhood(raw_start,
                                                                              raw_stop, 
                                                                              genome_length, 
                                                                              strict_span)
                except ValueError:
                    self.overlapping_termini = True
                if not self.overlapping_termini:
                    
                    #extract a biopython SeqRecord for the neighbourhood, 
                    #consisting of the scraped NCBI record from base neighbourhood_start 
                    #to base neighbourhood_end
                    self.neighbourhood = self.get_neighbourhood(self.neighbourhood_start,
                                                                self.neighbourhood_stop,
                                                                record = self.genome)
                    
    @staticmethod
    def scrape_genome(number_of_attempts : int, accession : str) -> SeqRecord:
        '''
        Scrape NCBI record associated with accession and parse into biopython 
        SeqRecord.

        Parameters
        ----------
        number_of_attempts : int
            Number of scrape attempts before giving up.
        accession : str
            Query accession.

        Raises
        ------
        HTTPError
            HTTPError associated with scraping.
        IncompleteRead
            IncompleteRead associated with reading.
        ValueError
            ValueError associated with reading.

        Returns
        -------
        SeqRecord
            NCBI record represented as SeqRecord object.

        '''
        
        #repeat until you run out or attempts or return a SeqRecord, then raise 
        #an error if no SeqRecord is found.
        for attempt in range(number_of_attempts):
            
            #Try to scrape the record from NCBI 
            try:
                handle = Entrez.efetch(db = "nucleotide", 
                                        id = accession, 
                                        rettype = "gbwithparts ", 
                                        retmode = "text")  
            
            #if scrape fails (e.g. accessiion is wrong)
            except HTTPError as e:
                if attempt == number_of_attempts - 1:
                    raise e
                else:
                    print ('Download error - trying again')
                    continue
            
            #Try to read the scraped file handle
            try:
                return SeqIO.read(handle, 
                                  "genbank")
            
            #if record cannot be parsed (e.g. accession was '' will rase a 
            #ValueError)
            except (IncompleteRead, ValueError) as e:
                if attempt == number_of_attempts - 1:
                    raise e
                else:
                    print ('Read error - trying again')
                    continue
    
    @staticmethod
    def define_neighbourhood(motif_start : int, motif_stop : int, 
                             neighbourhood_size : int) -> tuple:
        '''
        Given a motif_start/stop, define a neighbourhood (i) centered on the 
        midpoint of the motif and (ii) encompassing the entire motif.  If (ii) 
        is not possible, raise a ValueError. 

        Parameters
        ----------
        motif_start : int
            Start of first hit identified by cblaster in query record.
        motif_stop : int
            End of las hit identified by cblaster in query record..
        neighbourhood_size : int
            Size of neighbourhood to define (bp).

        Raises
        ------
        ValueError
            Motif length exceeds neighbourhood length.

        Returns
        -------
        tuple
            Neighbourhood start, neighbourhood end.

        '''
        motif_span = motif_stop - motif_start
        if motif_span > neighbourhood_size:
            raise ValueError()
        #round extension down to get an int - cannot have half a basepair
        extension = int((neighbourhood_size - motif_span)/2)
        neighbourhood_start = motif_start - extension
        neighbourhood_stop = motif_stop + extension
        return neighbourhood_start, neighbourhood_stop
    
    
    @staticmethod
    def sanitise_neighbourhood(raw_start : int, raw_stop : int, 
                               genome_length : int, strict_span : bool) -> tuple:
        '''
        Check the neiighbourhood start/end locci.  If they are impossible 
        (<0 or >record length respectively) then discard (if strict_span is True) 
        or correct them (to 0 and record legth respectively).

        Parameters
        ----------
        raw_start : int
            Defined neighbourhood start.
        raw_stop : int
            Defined neighbourhood end.
        genome_length : int
            Length of record.
        strict_span : bool
            If true, raise ValueError if raw_start or raw_end are incorrect.  
            Otherwise, correct raw_start and raw_end.

        Raises
        ------
        ValueError
            If strict_span is True and raw_start or raw_end are incorrect.

        Returns
        -------
        tuple
            either (raw_start, raw_end) or, if strict_span is False, 
            (corrected_start, corrected_stop).

        '''
        
        #START
        if raw_start < 0:
            if strict_span:
                raise ValueError()
            else:
               corrected_start = 0
        else:
            corrected_start = raw_start
        
        #STOP
        if raw_stop > genome_length:
            if strict_span:
                raise ValueError()
            else:
                corrected_stop = genome_length
        else:
            corrected_stop = raw_stop
            
        return corrected_start, corrected_stop
    
    @staticmethod 
    def get_neighbourhood(neighbourhood_start : int, neighbourhood_stop : int, 
                          record : SeqRecord) -> SeqRecord:
        '''
        Extract neighbourhood SeqRecord from parent SeqRecord.

        Parameters
        ----------
        neighbourhood_start : int
            Start of neighbourhood.
        neighbourhood_stop : int
            End of neighbourhood.
        record : SeqRecord
            Record in which the neighbourhood has been defined.

        Returns
        -------
        SeqRecord
            A biopython SeqRecord.

        '''
        #subsect record
        cut_record = record [neighbourhood_start : neighbourhood_stop]
        #copy over organism name etc
        cut_record.annotations = record.annotations
        return cut_record



def get_record(neighbourhood : Neighbourhood, scale : str) -> SeqRecord:
    '''
    Pull relevant SeqRecord attribute from Neighbourhood object.

    Parameters
    ----------
    neighbourhood : Neighbourhood
        A Neighbourhood with a neighbourhood SeqRecord and a genome (or contig) 
        SeqRecord.
    scale : str
        Which record should be returned (neighbourhood or genome).

    Raises
    ------
    ValueError
        Scale is not recognised.

    Returns
    -------
    SeqRecord
        User-specified SeqRecord.

    '''
    if scale == 'neighbourhood':
        return neighbourhood.neighbourhood
    elif scale == 'genome':
        return neighbourhood.genome
    else:
        raise ValueError(f"scale' must be 'neighbourhood' or 'genome', not {scale}")


def write_genbank_file(record : SeqRecord, path : str) -> None:
    '''
    Write SeqRecord to file.

    Parameters
    ----------
    record : SeqRecord
        SeqRecord to write.
    path : str
        Path to write SeqRecord.

    '''
    with open(path, 'w') as handle:
        SeqIO.write (record, 
                     handle, 
                     "genbank")

def make_filepath(name_type : str, type_map : dict, folder : str) -> str:
    '''
    Make a filepath for a Neighbourhood to be written in gbk format.  File is 
    named by either organism or accession of Neighbourhood's parent record, and 
    incremented if another genbank file of the same name has already been written. 

    Parameters
    ----------
    name_type : str
        Specifies what SeqRecord information should be used to name the files. 
        Must be a type_map key.
    type_map : dict
        Mapping of name_type to SeqRecord data 
        (e.g. {name_type : data, "organism" : "E. coli"})
    folder : str
        Folder that will contain the new file.

    Returns
    -------
    str
        Filepath.

    '''
    #name file as desired.  Note some organism names can contain slashes, so 
    #these are filtered out.  
    file = type_map[name_type]
    if '/' in file:
        file = file.replace('/', '_')
    if '\\' in file:
        file = file.replace('\\', '_')
    
    #Increment the file name if there is a file of the same name in folder 
    #(until that is no longer the case).
    temp_file = file
    count = 1
    while f'{temp_file}.gbk' in os.listdir(folder):
        temp_file = f'{file}_({count})'
        count += 1
    temp_file += '.gbk'
    
    return f'{folder}\\{temp_file}'

def write_results (results_folder : str, neighbourhood : Neighbourhood,
                   filenames : str, scale : str, logger_name : str) -> str:
    '''
    Make a scale-specific sub-folder in results_folder if one does not exist.  
    Make a filepath incorporating user-specified data.  Get a record of 
    user-specified scale.  Write record to filepath.

    Parameters
    ----------
    results_folder : str
        Folder that will contain the scale-specific subfolder.
    neighbourhood : Neighbourhood
        Neighbourhood with SeqRecord to be written.
    filenames : str
        What data should be used to name file (organism or accession).
    scale : str
        Write genome or neighbourhood SeqRecord to file.
    logger_name : str
        name of logger in which to record results

    Returns
    -------
    str
        Path to which SeqRecord was written.

    '''
    #make folder for scale if it doesnt exist
    folder = f'{results_folder}\\{scale}'
    os.makedirs(folder, 
                exist_ok = True)
    
    #make filename - note, you could get e.g. KeyError from a bad filenames value, 
    #but argument compatibility is confirmed at command line.
    path = make_filepath(name_type = filenames, 
                         type_map = {'accession' : neighbourhood.genome.id, 
                                     'organism' : neighbourhood.genome.annotations["organism"]},
                         folder = folder)
    
    #write record to filepath 
    record = get_record(neighbourhood, scale)
    write_genbank_file(record, path)
    
    #log results
    logger = logging.getLogger(logger_name)
    logger.info(f'written {neighbourhood.accession} {scale.upper()} to {path}')
    
    return path

def log_neighbourhood_details(neighbourhood : Neighbourhood, logger_name : str) -> str:
    '''
    Checks a Neoghbourhood object to see if it has been built successfully.  If
    it hasn't, a short description of the reason is returned for handling in 
    pipelines.sieve() and a full description is logged.

    Parameters
    ----------
    neighbourhood : Neighbourhood
        Neighbourhood being checked.
    logger_name : str
        name of logger in which to record results.

    Returns
    -------
    str
        Short description of logged message.

    '''
    logger = logging.getLogger(logger_name)
    if neighbourhood.scrape_error != None:
        logger.warning(f'scrape_fail - {neighbourhood.scrape_error} - {neighbourhood.accession}')
        return 'scrape_fail'
    elif neighbourhood.overlapping_termini:
        logger.warning(f'overlapping_termini - {neighbourhood.accession}')
        return 'overlapping_termini'
    elif neighbourhood.motif_to_long:
        logger.warning(f'motif is longer than specified neighbourhood - {neighbourhood.accession}')
        return 'long_motif'
    logger.info(f'neighbourhood = {neighbourhood.neighbourhood_start} -> '\
                    f'{neighbourhood.neighbourhood_stop}')
    return 'success'