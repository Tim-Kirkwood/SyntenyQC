# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:08:28 2024

@author: u03132tk
"""
from Bio import SeqIO
from Bio.Blast import NCBIXML
import os
import plotly.graph_objects as go
import networkx as nx
import shutil
import subprocess

def run_blast_process(cmd:list):
    '''this runs the BLAST+ BLASTP exe using Python via the subprocess module.  
    It uses error handlng to ease diagnostics - the return codes from BLAST+ 
    do not always match those of windows so (for example) 
    successful database creation can throw an error.'''
    try:
        return subprocess.run(cmd, check=True, capture_output=True, env={'BLASTDB_LMDB_MAP_SIZE':'1000000000'})
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

def makeblastdb_subprocess(makeblastdb_exe_path, input_fasta_path, db_out_path):#db_type_str param removed - redundant not used
    '''this runs the BLAST+ makeblastdb exe using Python via the subprocess module.  Input fasta is converted to BLAST db.  DB outpath is where you want the databse to be placed.
    NB - the DB outpath should have the full path up until where you want the files to be made, including the name of the files (path\to\dir\name).  
    Normally, name would have a suffix (e.g. .txt), but in this case do NOT include a file suffix.  
    Several files will be made with the same name and different suffixes, and adding a suffix does not appear to ipact databsae functionality, 
    but may lead the user to expect file functionality that does not exist - 
    ie adding a .txt suffix will not make the file behave like a .txt, as its final suffix wil be different.'''
    return run_blast_process([makeblastdb_exe_path, 
                              '-in', input_fasta_path, 
                              '-out', db_out_path, 
                              '-dbtype', 'prot'])
                              #'-parse_seqids'])#,,
                              #'-hash_index'])


def blastP_subprocess(evalue_user, query_user, blastp_exe_path, results_out_path, db_user, thread_num):
    '''this runs the biopython blastp wrapper, key details:  
        -xml output format 
        -parse deflines - this means you keep fasta heads, which is essential for coordinating rbh 
        -link to wrapper description - https://biopython.org/docs/1.75/api/Bio.Blast.Applications.html#Bio.Blast.Applications.NcbiblastpCommandline.  
    NB - the DB outpath should have the full path up until where you want the files to be made, including the name of the files (path\to\dir\name).  
    Normally, name would have a suffix (e.g. .txt), but in this case do NOT include a file suffix. 
    You are referenceing a collection of db files (.pdb, .phr, .pin, .pot, .psq, .ptf, .pto) not a single file.'''
    
    return run_blast_process([blastp_exe_path, 
                              '-out', fr'{results_out_path}', 
                              '-query', fr'{query_user}', 
                              '-db', fr'{db_user}', 
                              '-evalue', fr'{evalue_user}', 
                              '-outfmt', '5', 
                              '-num_threads', fr'{thread_num}',
                              #'-parse_deflines',#?
                              #https://microbiome.wordpress.com/research/orthologs/
                              #'-seg', 'yes',
                              #'-soft_masking', 'true',
                              #'-use_sw_tback',
                              #biopython defaults 
                              '-num_alignments', '200'])

def plot_graph(G, path):
    G = G.to_undirected()
    pos = nx.drawing.layout.spring_layout(G)
    nx.set_node_attributes(G, pos, 'pos')
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))
    node_degrees = []
    for adjacencies in G.adjacency():
        #print (adjacencies)
        node_degrees.append(len(adjacencies[1]))
    #TODO confirm that G.adjacency and G.nodes are same order
    node_text = []
    for degree, node in zip(node_degrees, G.nodes):
        node_text.append(f"{node} # of connections: {degree}")
    node_trace.marker.color = node_degrees
    node_trace.text = node_text
    fig = go.Figure(data=[edge_trace, node_trace],
                 layout=go.Layout(
                    title='<br>Network graph made with Python',
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    annotations=[ dict(
                        text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                        showarrow=False,
                        xref="paper", yref="paper",
                        x=0.005, y=-0.002 ) ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    fig.write_html(path)
    
def write_fasta(output_filepath, fasta_tuple):
    with open (output_filepath, 'w') as outfile:
        for index, (defline, seq) in enumerate(fasta_tuple):
            assert defline[0] == '>'
            outfile.write(defline + '\n')
            outfile.write(seq)
            if index < len(fasta_tuple) - 1:
                outfile.write('\n')

def build_fasta(genbank_folder):
    genbank_files = get_gbk_files(genbank_folder)
    fasta = []
    for file in genbank_files:
        cds_count = 0
        with open(f'{genbank_folder}\\{file}', 'r') as handle:
            record = SeqIO.read(handle, 
                                "genbank")
        record_id = file[0:file.rindex('.')]
        for feature in record.features:
            if feature.type == 'CDS':
                
                try:
                    protein_seq = ''.join(feature.qualifiers['translation'])
                except KeyError:
                    if 'pseudo' in feature.qualifiers.keys():
                        continue
                    else:
                        raise KeyError
                fasta += [(f'>{record_id}__{cds_count}', 
                           protein_seq)
                          ]
                cds_count += 1
    return fasta   
    
    

def build_neighbourhood_size_map(genbank_folder):
    genbank_files = get_gbk_files(genbank_folder)
    size_map = {}
    for file in genbank_files:
        with open(f'{genbank_folder}\\{file}', 'r') as handle:
            record = SeqIO.read(handle, 
                                "genbank")
        record_id = file[0:file.rindex('.')]
        size_map[record_id] = 0
        for feature in record.features:
            if feature.type == 'CDS':
                if 'pseudo' in feature.qualifiers.keys():
                    continue
                size_map[record_id] += 1
    return size_map

def find_blastp_exe(cmd, version):
    exe_path = shutil.which(cmd)
    if version not in exe_path:
        raise ValueError(f'incorrect {cmd} version stored in path - {exe_path}')
    elif exe_path is None:
        raise ValueError('could not find BLAST+ (version {version}) {cmd} in PATH - please check you installed BLAST+ correctly')
    return exe_path

def run_all_vs_all_blast(folder_with_genbanks, 
                         blast_file_folder, 
                          
                         
                         e_value):
    
    #folder = f'{folder_with_genbanks}\\ClusterSieve_results'
    all_proteins = f'{blast_file_folder}\\all_proteins.txt' 
    all_proteins_db = f'{blast_file_folder}\\all_proteins'
    results_out_path =  f'{blast_file_folder}\\results.xml' 
    makeblastdb_exe_path = find_blastp_exe(cmd = 'makeblastdb',
                                           version = 'blast-2.10.1+')
    blastp_exe_path = find_blastp_exe(cmd = 'blastp',
                                      version = 'blast-2.10.1+')
    #query_user = f'{folder}/query.txt'  
    fasta = build_fasta(genbank_folder=folder_with_genbanks)
    write_fasta(output_filepath = all_proteins, 
                fasta_tuple = fasta)
    print (f'Running all vs all BLASTP - {len(fasta)} protein sequences...')    
    
    # query_user = 'C:/Users/u03132tk/.spyder-py3/boundary_dataset/NZ_CP042324.1.region011.fasta'
    
    makeblastdb_subprocess(makeblastdb_exe_path, 
                            all_proteins, 
                            all_proteins_db)
    print ('Made database...')
    blastP_subprocess (e_value, 
                        all_proteins, 
                        blastp_exe_path, 
                        results_out_path, 
                        all_proteins_db, 
                        4)
    print ('Completed BLASTP...')
    return results_out_path

def make_hit_matrix(results_path, min_percent_identity):
    raw_results = {}#
    with open(results_path, 'r') as result_handle:
        blast_records = list(NCBIXML.parse(result_handle))
        for record in blast_records:
            query_scaffold, query_index = record.query.split('__')
            if query_scaffold not in raw_results.keys():
                raw_results[query_scaffold] = {}
            if query_index not in raw_results[query_scaffold].keys():
                raw_results[query_scaffold][query_index] = {}
            for hit in record.alignments:
                hit_scaffold, hit_index = hit.hit_def.split('__')
                best_hsp = None
                for hsp in hit.hsps:
                    percent_identity = 100*(hsp.identities / hsp.align_length) 
                    #query_len = protein_length_map[f'{query_scaffold}_{query_index}']
                    #query_coverage = hsp.align_length/query_len
                    if percent_identity < min_percent_identity:
                        continue
                    if best_hsp == None:
                        best_hsp = hsp
                    else:
                        if hsp.score > best_hsp.score:
                            best_hsp = hsp
                
                if hit_scaffold not in raw_results[query_scaffold][query_index].keys():
                    raw_results[query_scaffold][query_index][hit_scaffold] = {}
                raw_results[query_scaffold][query_index][hit_scaffold][hit_index] = best_hsp
    return raw_results


def make_best_hit_matrix(hit_matrix):
    best_hits = {}
    for query_scaffold, query_proteins in hit_matrix.items():
        best_hits[query_scaffold] = {}
        for query_index, hit_proteins in query_proteins.items():
            best_hits[query_scaffold][query_index] = {}
            for hit_scaffold, hit_proteins in hit_proteins.items():
                best_hit = None
                best_score = 0
                for hit_index, hit_hsp in hit_proteins.items():
                    if hit_hsp == None:
                        continue
                    if hit_hsp.score > best_score:
                        best_hit = hit_index
                        best_score= hit_hsp.score
                if best_hit == None:
                    #no hsps that meet coverage lims
                    continue
                best_hits[query_scaffold][query_index][hit_scaffold] = best_hit
    return best_hits


def make_reciprocal_best_hits_matrix(best_hit_matrix):
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

def make_graph(reciprocal_best_hit_matrix, 
               neighbourhood_size_map, 
               conserved_neighbourhood_len,
               similarity_filter):
    #do not draw an edge between two neighbourhood nodes that share fewer RBHs than specified by conserved_neighbourhood_len
    #do not include neighbourhoods with fewer proteins than conserved_neighbourhood_len
    small_neighbourhoods = []
    edge_map = {}
    for record_accession in reciprocal_best_hit_matrix.keys():
        edge_map[record_accession] = {}
        for other_record_accession in reciprocal_best_hit_matrix.keys():
            if record_accession != other_record_accession:
                number_of_rbh = 0
                smallest_record_protein_count = min(neighbourhood_size_map[record_accession], 
                                                    neighbourhood_size_map[other_record_accession])
                
                for query_protein, rbh_scaffolds in reciprocal_best_hit_matrix[other_record_accession].items():
                    if record_accession in rbh_scaffolds.keys():
                        number_of_rbh += 1
                corrected_rbh_number = number_of_rbh - conserved_neighbourhood_len
                if corrected_rbh_number < 0:
                    continue
                corrected_neighbourhood_len = smallest_record_protein_count - conserved_neighbourhood_len
                if corrected_neighbourhood_len <0:
                    if neighbourhood_size_map[record_accession] - conserved_neighbourhood_len < 0:
                        small_neighbourhoods += [record_accession]
                    if neighbourhood_size_map[other_record_accession]  - conserved_neighbourhood_len < 0:
                        small_neighbourhoods += [other_record_accession]
                    continue
                similarity_score = corrected_rbh_number/corrected_neighbourhood_len
                if similarity_score > similarity_filter:
                    edge_map[record_accession][other_record_accession] = {'weight' : similarity_score} 
    G = nx.Graph(edge_map)
    all_possible_nodes = set(reciprocal_best_hit_matrix.keys()) - set(small_neighbourhoods)
    G.add_nodes_from(all_possible_nodes)
    return G

def prune_graph(graph):
    nodes = graph.nodes 
    while True:
        temp_G = graph.subgraph(nodes)
        temp_nodes = []
        degrees = []
        for node, degree in temp_G.degree:
            temp_nodes += [node]
            degrees += [degree]
        max_degree = max(degrees)
        if max_degree == 0:
            break
        delete_index = degrees.index(max_degree)#ONLY TAKE ONE NODE, NOT ALL NODES THAT HAVE SAME DEGREE
        nodes = [n for i, n in enumerate(temp_nodes) if i != delete_index]
    return nodes
def get_gbk_files(folder):
    files = os.listdir(folder)
    suffixes = ['.gbk','.gb'] 
    gbk_files = []
    for file in files:
        if os.path.isfile(f'{folder}\\{file}'):        
            if file[file.rindex('.') : ] in suffixes:
                gbk_files += [file]
    return gbk_files

def write_nodes(nodes, genbank_folder, results_folder):
    #'check old and new paths arent the same and that 
    
    gbk_files = get_gbk_files(genbank_folder)
    written_nodes = []
    for file in gbk_files:
        #with open(f'{genbank_folder}\\{file}', 'r') as query_handle:
        #    record = SeqIO.read (query_handle, 
        #                            "genbank")
        record_id = file[0:file.rindex('.')]
        if record_id in nodes:
            this_path = f'{genbank_folder}\\{file}'
            new_path = f'{results_folder}\\{file}'
            shutil.copy(src = this_path, 
                        dst = new_path)
            #note you are working from a directory so every file has to be unique
            #chekcing/updating filenames will occur in the clustercollect stage
            # with open(f'{results_folder}\\{file}', 'w') as result_handle:
            #     SeqIO.write (record, 
            #                  result_handle, 
            #                  "genbank")
            written_nodes += [record_id]
    if set(written_nodes) != set(nodes):
        raise ValueError('Not all nodes are found in the genbank folder - this should not be possible...')

def sieve(genbank_folder, 
                  
                 e_value, 
                 min_percent_identity, 
                 conserved_neighbourhood_len, 
                 similarity_filter):
    #check slashes are standardised  and do os.make_dirs
    folder_with_genbanks = genbank_folder.replace('/', '\\')
    current_dir = os.path.dirname(os.path.abspath(__file__))
    blast_file_folder = f'{current_dir}\\BLASTP_files'
    plot_path = f'{folder_with_genbanks}\\RBH_graph.html'
    results_folder = f'{folder_with_genbanks}\\ClusterSieve'
    
    #TODO - other blastp params
    results_path = run_all_vs_all_blast(folder_with_genbanks, 
                                        blast_file_folder, 
                                         
                                        e_value)
    #results_folder = f'{genbank_folder}\\{ClusterSieve_results}' 
    hit_matrix = make_hit_matrix(results_path, 
                                 min_percent_identity)
    print ('Processed hits...')
    best_hit_matrix = make_best_hit_matrix(hit_matrix)
    print ('Processed best hits')
    reciprocal_best_hit_matrix = make_reciprocal_best_hits_matrix(best_hit_matrix)
    print ('Processed reciprocal best hits')
    neighbourhood_size_map = build_neighbourhood_size_map(folder_with_genbanks)
    G = make_graph(reciprocal_best_hit_matrix, 
                   neighbourhood_size_map, 
                   conserved_neighbourhood_len,
                   similarity_filter)
    #sort file - put in genbank folder
    plot_graph(G, plot_path)
    print (f'Made RBH graph - written to {plot_path}') 
    nodes = prune_graph(G)
    print ('Pruned graph')
    os.makedirs(results_folder, 
                exist_ok=True)
    write_nodes(nodes, 
                folder_with_genbanks, 
                results_folder)
    