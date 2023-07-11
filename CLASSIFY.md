# Taxonomic classification of human gut viruses

The code and database described here will allow you to obtain a taxonomic label for your gut virus based on the UHGV taxonomy. This is useful to determine novelty relative to database, identify characteristics of the nearest viral group, and allow you to identify other phylogenetically related viruses in the database.

<img src="img/classify_workflow.png" width="900">

## Installation

Install program using git and pip (add `--user` if you don't have root access):  
`pip install git+https://github.com/snayfach/UHGV.git`

Install external dependencies using conda:  
`conda install -c bioconda prodigal-gv diamond blast -y`

View available modules:  
`uhgv-tools -h`

Download and unpack the latest database:   
`uhgv-tools download .`

> UHGV-tools: download<br>
> [1/5] Checking latest version of database...<br>
> [2/5] Downloading 'uhgv-db'...<br>
> [3/5] Extracting 'uhgv-db'...<br>
> [4/5] Building BLASTN database...<br>
> [5/5] Building DIAMOND database...<br>
> Run time: 121.6 seconds<br>
> Peak mem: 1.58 GB<br>

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

Download a test dataset of 5 phages from [Nishijima et al.](https://www.nature.com/articles/s41467-022-32832-w) using wget:  
`wget https://raw.githubusercontent.com/snayfach/UHGV/main/example/viral_sequences.fna -O viral_sequences.fna`

Classify sequences, replacing `</path/to/uhgv-db>` as appropriate:   
`uhgv-tools classify -i viral_sequences.fna -o output -d </path/to/uhgv-db> -t 10`

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

Here are field definitions and example values for `classify_summary.tsv`:

| Field      | Description | Example |
| ----------- | ----------- |----------- |
| genome_id | user genome identifier | 0008\_k141_99927 |
| genome_length | length in bp | 96989 |
| genome\_num_genes | count of CDS | 106 |
| taxon_id | UHGV taxon identifier | vSUBGEN-22354 |
| class_method | nucleotide or protein based classification | protein |
| class_rank | lowest classified rank | subgenus |
| ani_reference | nearest reference based on ANI | UHGV-0030436 |
| ani_identity | nucleotide identity | 93.65 |
| ani\_query_af | % of query covered | 86.6 |
| ani\_target_af | % of target covered | 83.58 |
| ani_taxonomy | taxonomy of reference genome | vFAM-00050;vSUBFAM-00057;vGENUS-00180;vSUBGEN-22354;vOTU-000988 |
| aai_reference | nearest reference based on AAI | UHGV-0030436 |
| aai\_shared_genes | number of proteins aligned | 93 |
| aai_identity | amino acid identity | 89.33 |
| aai_score | normalized, cumulative bitscore | 82.57 |
| aai_taxonomy | taxonomy of reference genome | vFAM-00050;vSUBFAM-00057;vGENUS-00180;vSUBGEN-22354;vOTU-000988 |

Here field definitions and example values for `taxon_info.tsv`:

| Field      | Description | Example |
| ----------- | ----------- |----------- |
| genome_id | user genome identifier | 0008\_k141_99927 |
| taxon_id | UHGV taxon identifier | vSUBGEN-22354 |
| taxon\_lineage | UHGV taxon lineage | vFAM-00050;vSUBFAM-00057;vGENUS-00180;vSUBGEN-22354 |
| host\_lineage | Consensus GTDB host lineage | d\_\_Bacteria;p\_\_Bacteroidota;c\_\_Bacteroidia;o\_\_Bacteroidales;f\_\_Bacteroidaceae;g\_\_Prevotella (100.0) |
| ictv\_lineage | Consensus ICTV taxon lineage | r\_\_Duplodnaviria;k\_\_Heunggongvirae;p\_\_Uroviricota;c\_\_Caudoviricetes;o\_\_Crassvirales;f\_\_Beta-crassviridae (100.0) |
| lifestyle | Consensus virus lifestyle | virulent (100.0) |
| genome\_length_median | median genome length of viruses in lineage | 100566.0 |
| genome\_length_iqr | interquartile range of genome length | 100566.0 - 100566.0 |

                 
        
