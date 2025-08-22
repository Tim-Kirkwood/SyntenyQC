# This app has been tested on Mac, Windows and Linux.

# SyntenyQC
## Motivation: 
Synteny plots are widely used for the comparison of genomic neighbourhoods.  Whilst synteny plots are often included as part of larger software suites (e.g. the `antiSMASH` ClusterBlast module), various low-code, stand-alone tools are now available that allow users to source candidate neighbourhoods and build their own synteny plots.  However, a gap remains between: 

**(i) tools that source these candidate neighbourhoods** (e.g. `cblaster`, which can find hundreds of candidates), 

**(ii) tools that build the synteny plots** (e.g. `clinker`, which struggles as the number of neighbourhoods exceeds 30-50) and 

**(iii) the synteny plots themselves**, which become much harder to analyse/present as the number of neighbourhoods they include increases.

## Description: 
`SyntenyQC` is a python app for the curation of neighbourhoods immediately prior to synteny plot creation. `SyntenyQC collect` supports the systematic definition and annotation of candidate neighbourhoods based on a direct integration to `cblaster`.  `SytenyQC sieve` offers a flexible method for objectively removing redundant neighbourhoods (sourced using `cblaster` or any other tool) prior to synteny plot creation.  This is in some cases an absolute requirement (e.g. `cblaster` called via the `CAGECAT` webserver places a limit of 50 neighbourhoods).  

## Installation 
**What do I do?**

```
conda create --name syntenyqc_env pip python=3.12.9
conda activate syntenyqc_env
pip install SyntenyQC
```
 
**Why do I do it?**

You should install SyntenyQC within a virtual environment to make sure it doesn't interfere with any other software you have installed ([read more here](https://stackoverflow.com/q/41972261/11357695)).  There are various options for working with virtual environments, but I use `conda` - see their tutorial [here](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html#).  As `SyntenyQC` is uploaded to the Python Package Index ([PyPI](https://pypi.org/)) and not Anaconda (yet - also, how is [`conda` different from Anaconda](https://stackoverflow.com/questions/30034840/what-are-the-differences-between-conda-and-anaconda) and [whats the relationship between `pip` and PyPI?](https://stackoverflow.com/questions/74307171/does-pip-only-use-pypi-or-does-it-use-other-domains-to-find-packages)), we will set up an environment using `conda`, install `pip` in that environment, and then install `SyntenyQC` using `pip` (note, Python can be any version between 3.10.0 and 3.12.9 inclusive).  

**Software you need to install manually**

`SyntenyQC` depends on [DIAMOND](https://github.com/bbuchfink/diamond), which must be [installed](https://github.com/bbuchfink/diamond/wiki/2.-Installation) by the user (tested with v2.1.12.166 - but should work with other versions unless there are parameter changes). If this is installed correctly, you should be able to see help messages after typing `diamond help` (with no "-" or "--") in the command line. 

Note - after you download the diamond executable file (.exe), you will probably need to add it to your [path](https://stackoverflow.com/a/56929848/11357695) - this allows your computer to understand what you mean when you type `diamond help` into the command line.  When you add to DIAMOND to your path, `diamond help` becomes equivalent to `path/to/diamond.exe help`. This is [easy to do](https://www.eukhost.com/kb/how-to-add-to-the-path-on-windows-10-and-windows-11/) - on Windows, just go to **Start**, search for **Edit the system environment variables**, click **Environment Variables**, under **User variables** click **Path** and then **Edit**, then finally add the **path to the folder with the exe** (not the .exe filepath) to the Path variable.     

## Tests
Tests are performed using [pytest](https://pypi.org/project/pytest/), but are not distributed with `SyntenyQC`.  To run tests:

1) Install pytest
2) Clone the `SyntenyQC` github repository.
3) Update `path/to/cloned/repository/tests/email.txt` with your email (if left blank, two webscraping tests will fail).
4) Navigate to the cloned repository via command line (on Windows, type `cd path/to/cloned/repository`).
5) Type `pytest` in the command line.
 
## Usage
### General help:
```
>SyntenyQC -h
usage: SyntenyQC [-h] {collect,sieve} ...

options:
  -h, --help       show this help message and exit

subcommands:
  {collect,sieve}  Synteny quality control options
```
### Collect subcommand:
```
>SyntenyQC collect -h
usage: SyntenyQC collect [-h] -bp -ns -em [-fn] [-sp] [-wg]

Write genbank files corresponding to cblaster neighbourhoods from a specified CSV-format binary file located at
BINARY_PATH.  For each cblaster hit accession in the binary file:

1) A record is downloaded from NCBI using the accession.  NCBI requires a user EMAIL to search for this record
   programatically.  If WRITE_GENOMES is specified, this record is written to a local file according to FILENAMES
   (see final bulletpoint).
2) A neighbourhood of size NEIGHBOURHOOD_SIZE bp is defined, centered on the cblaster hits defined in the binary
   file for the target accession.
3) (If STRICT_SPAN is specified:) If the accession's record is too small to contain a neighbourhood of the
   desired size, it is discarded.  For example, if an accession record is a 25kb contig and NEIGHBOURHOOD_SIZE
   is 50000, the record is discarded.
4) If FILENAMES is "organism", the nighbourhood is written to file called *organism*.gbk. If FILENAMES is
   "accession", the neighbourhood is written to *accession*.gbk. Synteny softwares such as clinker can use these
   filenames to label synetny plot samples.

Once COLLECT has been run, a new folder with the same name as the binary file will be seen in the directory
that holds the binary file (i.e. the file "path/to/binary/file.txt" will generate the folder "path/to/binary/file").
This folder will have a subdirectory called "neighbourhood", containing all of the neighbourhood genbank files
(i.e. "path/to/binary/file/neighbourhood"). If WRITE_GENOMES is specified, a second direcory ("genome") will also
be present, containing the entire record associated with each cblaster accession (i.e. "path/to/binary/file/genome").
Finally, a log file will be present in the folder "path/to/binary/file", containing all run details.

options:
  -h, --help            show this help message and exit
  -bp, --binary_path
                        Full filepath to the CSV-format cblaster binary file containing neighbourhoods that should
                        be extracted
  -n, --neighbourhood_size
                        Size (basepairs) of neighbourhood to be extracted (centered on middle of CBLASTER-defined
                        neighbourhood)
  -em, --email          Email - required for NCBI entrez querying
  -fn, --filenames
                        If "organism", all collected files will be named according to organism. If "accession", all
                        files will be named by NCBI accession. (default:
                        organism)
  -sp, --strict_span
                        If set, will discard all neighbourhoods that are smaller than neighbourhood_size bp. For
                        example, if you set a neighbourhood_size of 50000, a 50kb neighbourhood will be extracted
                        from the NCBI record associateed with each cblaster hit. If the record is too small for this
                        to be done (i.e. the record is smaller then 50kb) it is discarded.
  -wg, --write_genomes
                        If set, will write entire NCBI record containing a cblaster hit to file (as well as just the
                        neighbourhood)
```
### Sieve subcommand:
```
>syntenyqc sieve -h
usage: SyntenyQC sieve [-h] -gf [-ev] [-mts] [-mev] [-sf] [-am] [-dmts] [-ex] [-qc] [-sc] [-id] [-ks]

Filter redundant genomic neighbourhoods based on neighbourhood similarity:
- First, an all-vs-all BLASTP is performed with user-specified BLASTP settings and the neighbourhoods in GENBANK_FOLDER.
- Secondly, these are parsed to define reciprocal best hits between every pair of neighbourhoods.
- Thirdly, these reciprocal best hits are used to derive a neighbourhood similarity network, where edges indicate two
  neighbourhood nodes that have a similarity > SIMILARITY_FILTER. Similarity = Number of RBHs / Number of proteins in
  smallest neighbourhood in pair.
- Finally, this network is pruned to remove neighbourhoods that exceed the user's SIMILARITY_FILTER threshold. Nodes
  that remain are copied to the newly created folder "genbank_folder/sieve_results/genbank".

options:
  -h, --help            show this help message and exit
  -gf, --genbank_folder
                        Full path to folder containing neighbourhood genbank files requiring de-duplication
  -ev, --e_value    BLASTP evalue threshold. (default: 1e-05)
  -mts, --max_target_seqs
                        BLASTP -max_target_seqs. Maximum number of aligned sequences to keep. (default: 200)
  -mev, --min_edge_view
                        Minimum similarity between two neighbourhoods for an edge to be drawn betweeen them in the
                        RBH graph. Purely for visualisation of the graph HTML file - has no impact on the graph pruning
                        results. (default: --similarity_filter)
  -sf, --similarity_filter
                        Similarity threshold above which two neighbourhoods are considered redundant (default: 0.7)
  -am, --alignment_mode
                        Alignment mode used by DIAMOND (choices: fast, mid-sensitive, sensitive, more-sensitive,
                        very-sensitive, ultra-sensitive). Without using any sensitivity option, the default mode
                        will run which is designed for finding hits of >60 percent identity and short read alignment.
                        Its sensitivity is between --fast and --mid-sensitive. See here
                        https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes
  -dmts, --dynamic_max_target_seqs
                        If set, --max_target_seqs will be automatically defined as the numer of genbank files
                        within --genbank_folder or --max_target_seqs, whichever is larger
  -ex, --expand         If set, DO NOT gzip compress DIAMOND results file (will increase disk space requirments)
  -qc, --query_cover
                        Report only alignments above the given percentage of query cover. Note that using this option
                        reduces performance.
  -sc, --subject_cover
                        Report only alignments above the given percentage of subject cover. Note that using this option
                        reduces performance.
  -id, --identity   Report only alignments above the given percentage of sequence identity. Note that using this option
                    reduces performance.
  -ks, --keep_pseudo    if set, will count pseudo entries (or missing sequences) when counting the number of proteins
                        within a given neighbourhood for the inter-neighbourhood similarity score.
```
### Sieve pruning algorithm:
```
Data: RBH graph
Result: Nodes from pruned RBH graph
Procedure:
    while max(node degrees in RBH graph) > 0:
        delete nodes = []
        for node in RBH graph:
            if node degree == max(node degrees in RBH graph):
                delete nodes + node
        delete node = random node from delete nodes
        RBH graph = RBH graph - delete node 
    return nodes in RBH graph 
```

# Example use:
## (1) Preliminaries
Take the [actinorhodin BGC from MIBIG](https://mibig.secondarymetabolites.org/repository/BGC0000194/index.html#r1c1), and analyse with `cblaster` via the command line or `CAGECAT` web server. This will generate a set of files, one of which will be a 'binary file', stored at `folder/with/binary.csv`. To show that even stringent searches can generate many hits, all core biosynthetic genes (`-r CAC44200.1;CAC44201.1`) and >= 11 hits (`-mh 11`) were required for each hit neighbourhood.

**⚠️WARNING: To use a binary file in `SyntenyQC collect`, `cblaster` must be run with a ',' binary delimeter and no intermediate genes (`-bde ','` and no `-ig` flag)⚠️**  
## (2) Collect neighbourhoods  
### Starting directory structure:
```
folder/with/binary.csv
```
        
### Command (neighbourhood size 42566kb = 2 x BGC length):
     
```
#Note - you should add your email
SyntenyQC collect -bp path/to/BGC0000194_binary.txt -ns 42566 -em my_email@domain.com -fn organism -sp -wg
```
   
### Finishing directory structure: 
```
folder/with/binary/neighbourhood/organism1.gbk, organism2.gbk...organism157.gbk
                  /genomes      /organism1.gbk, organism2.gbk...organism157.gbk         ###only if -wg!
                  /log.txt
```
#### 🔴 157 neighbourhoods is a lot for a synteny plot 🔴


## (3) Sieve neighbourhoods
### Starting directory structure 
From `SyntenyQC collect` in this example, but any folder with at least one genbank file can be used:
```
folder/with/binary/neighbourhood/organism1.gbk, ...
```
### Command:
```
SyntenyQC sieve -gf folder/with/binary/neighbourhood
```
### Finishing directory structure: 
```
folder/with/binary/neighbourhood/organism1.gbk, ...
                                /sieve_results/blastp        /results.xml, db.txt, db.pin, ...   #call all be deleted
                                              /genbank       /full_adjacency.csv, sieved_adjacency.csv, organism1.gbk, ...organism38.gbk   #use as e.g. clinker input
                                              /visualisations/RBH_graph.html, RBH_histogram.html #see what is being pruned
                                              /log.txt
```
#### :green_heart: 38 neighbourhoods is OK for a synteny plot :green_heart:
## Notes 
- The above example was performed on an earlier `SyntenyQC` version (1.0) using `BLASTP` from `BLAST+` instead of `DIAMOND`.  Everything is the same command/filepath wise (there are some new parameters, but they have defaults so the above commands would work).  However, the number of filtered neighbourhoods may be different if you try and replicate this with `SyntenyQC` version 2.  See our Application Note for neighbourhod count data using `SyntenyQC` version 2.
- Given filenames are purely to show number of files - `neighbourhood/sieve_results/genbank/organism1.gbk` is one of the genbanks in the `neighbourhood` folder, but may be different to     `neighbourhood/organism1.gbk`.
- `RBH_graph` is an interactive html picture of the similarity graph created by `SyntenyQC sieve` (before pruning), only showing edges with a similarity > `min_edge_view`. 
 Edges are black (< `similairty_filter`) or red (>= `similarity_filter`).  The associated adjacency matricies for this graph before and after pruning are also saved in CSV format.  Note, unless `--min_edge_view` is set to < `--similarity_filter`, the edges in `sieved_adjacency.csv` will all be 0 (only similarities > `--min_edge_view` are included in the initial graph, so if `--min_edge_view` = `--similarity_filter` then all edges will be removed by `SyntenyQC Sieve`.).  
- `RBH_histogram` shows the distribution of edge weights.
- Most neighbourhoods that meet the user-defined similarity threshold `-sf` will be removed in a single sieve run.  However, the `sieve -mts` setting can impact final results.  If a protein has homologs in 251 neighboughoods and `-mts` is 250, then one of the homologs will be missed by BLASTP, and the host neighbourhood may appear less similar to the neighbourhood with the query.  Whilst a high `-mts` setting could be used for `sieve` runs pruning many neighbourhoods, this will generate large blast files that may take up a lot of space and slow down the run.  Thus, if users wish to remove the (typically 2-3) redundant neighbourhoods remaining after a single `sieve` call, they can run sieve again on the pruned results (`syntenyqc sieve -g path/to/sieve_results/genbank -sf 0.7`).  See our paper, Supplementary Methods S2, for a suggested search strategy that optimises disk space and run time for large searches.   

## References
### `cblaster`

**Paper:** Cameron L M Gilchrist, Thomas J Booth, Bram van Wersch, Liana van Grieken, Marnix H Medema, Yit-Heng Chooi, cblaster: a remote search tool for rapid identification and visualization of homologous gene clusters, Bioinformatics Advances, Volume 1, Issue 1, 2021, vbab016, https://doi.org/10.1093/bioadv/vbab016 

**Docs:** https://cblaster.readthedocs.io/en/latest/ 

### `clinker` 

**Paper:** Cameron L M Gilchrist, Yit-Heng Chooi, clinker & clustermap.js: automatic generation of gene cluster comparison figures, Bioinformatics, Volume 37, Issue 16, August 2021, Pages 2473–2475, https://doi.org/10.1093/bioinformatics/btab007 

**Docs:** https://github.com/gamcil/clinker 

### `CAGECAT`

**Paper:** van den Belt, M., Gilchrist, C., Booth, T.J. et al. CAGECAT: The CompArative GEne Cluster Analysis Toolbox for rapid search and visualisation of homologous gene clusters. BMC Bioinformatics 24, 181 (2023). https://doi.org/10.1186/s12859-023-05311-2

**Website:** https://cagecat.bioinformatics.nl/ 

### `antiSMASH` 

**Paper (ClusterBlast was introduced in version 1):** Medema MH, Blin K, Cimermancic P, de Jager V, Zakrzewski P, Fischbach MA, Weber T, Takano E, Breitling R. antiSMASH: rapid identification, annotation and analysis of secondary metabolite biosynthesis gene clusters in bacterial and fungal genome sequences. Nucleic Acids Res. 2011 Jul;39(Web Server issue):W339-46. doi: 10.1093/nar/gkr466. Epub 2011 Jun 14. PMID: 21672958; PMCID: PMC3125804

**Website (latest version):** https://antismash.secondarymetabolites.org/#!/start 
