# Taxonomic classification of environmental phages and viruses from novel lineages

Viral database: 
- Currently the database is focused on viruses from the human gut

Use cases:
- Determine novelty: does my viral genome represent a novel species? genus? family?
- Ecological analysis: compare viral phylogenetic groups across samples
- Comparative genomes: retrieve other viral genomes from the same phylogenetic group
- Infer host and lifestyle: impute characteristics of the virus based on it's phylogenetic group
- Update the UHGV: cluster unclassified viral genomes into de novo vOTUs

## Installation

Install program using git and pip (add `--user` if you don't have root access):  
`pip install git+https://github.com/snayfach/UHGV-toolkit.git`

Install external dependencies using conda:  
`conda install -c bioconda prodigal-gv diamond blast -y`

View available modules:  
`uhgv-tools -h`

Download and unpack the latest database:   
`uhgv-tools download_database`

View command line usage for `classify` module:  
`uhgv-tools classify -h`

>usage: uhgv-class.py [-h] -i PATH -o PATH -d PATH [-t THREADS] [-c]
>
>options:
>  -h, --help  show this help message and exit
>
>required arguments:
>  -i PATH     Path to nucleotide seqs<br>
>  -o PATH     Path to output directory<br>
>  -d PATH     Path to database directory<br>
>  -t THREADS  Number of threads to run program with (1)<br>
>  --continue  Continue where program left off<br>
>  --quiet     Suppress logging messages<br>

Download a test dataset:  
`wget XXXX`

Classify contigs from the test dataset:  
`uhgv-tools classify -i test_sequences.fna -o output -d uhgv-db-v0.1/`

Main results are found in `output/classify_summary.tsv`:

- genome : 
- classification : 
- classification_method : 
- ani_reference : 
- ani_identity : 
- ani_query_af : 
- ani_target_af : 
- ani_taxonomy : 
- aai_reference : 
- aai_shared_genes : 
- aai_identity : 
- aai_score : 
- aai_taxonomy : 
