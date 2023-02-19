# Taxonomic classification of human gut viruses

The code and database described here will allow you to obtain a taxonomic label for your gut virus based on the UHGV taxonomy. This is useful to determine novelty relative to database, identify characteristics of the nearest viral group, and allow you to identify other phylogenetically related viruses in the database.

<img src="img/classify_workflow.png" width="900">

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

## Example usage

Download the test dataset consisting of 5 phage genomes from [Nishijima et al.](https://www.nature.com/articles/s41467-022-32832-w)  

If you've cloned the repo, these are found in `UHGV/example/viral_sequences.fna`

Otherwise, download using wget:  
`wget https://raw.githubusercontent.com/snayfach/UHGV/main/example/viral_sequences.fna?token=GHSAT0AAAAAAB5YRYUZ2FVVNDY5NIHRTC44Y7SNMYA -O viral_sequences.fna`

Classify sequences, replacing `</path/to/uhgv-db>` as appropriate:   
`uhgv-tools classify -i viral_sequences.fna -o output -d </path/to/uhgv-db> -t 10`

Expected logging messages:

> UHGV-tools v0.0.1: classify<br>
> [1/10] Reading input sequences<br>
> [2/10] Reading database sequences<br>
> [3/10] Estimating ANI with blastn<br>
> [4/10] Identifying genes using prodigal-gv<br>
> [5/10] Performing self alignment<br>
> [6/10] Aligning proteins to database<br>
> [7/10] Calculating amino acid similarity scores<br>
> [8/10] Finding top database hits<br>
> [9/10] Performing phylogenetic assignment<br>
> [10/10] Writing output file(s)<br>

There are two main output files:

- `output/classify_summary.tsv`: information related to classification 
- `output/taxon_info.tsv`: details about the classified taxa (ex: lifestyle, genome size, host)

Here are field definitions for `classify_summary.tsv`:

- genome : 
- classification : 
- classification\_method : 
- ani\_reference : 
- ani\_identity : 
- ani\_query_af : 
- ani\_target_af : 
- ani\_taxonomy : 
- aai\_reference : 
- aai\_shared_genes : 
- aai\_identity : 
- aai\_score : 
- aai\_taxonomy : 

Here are field definitions for `taxon_info.tsv`:
