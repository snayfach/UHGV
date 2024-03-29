# Unified Human Gut Virome Catalog (UHGV)

The UHGV is a comprehensive genomic resource of viruses from the human microbiome. Genomes were derived from [12 independent data sources](#data-sources) and annotated using a [uniform bioinformatics pipeline](#bioinformatics-pipeline):

<img src="img/data_workflow.png" width="900">


## Table of contents
1. [Methods](#methods)
   * [Data sources](#data-sources)
   * [Bioinformatics pipeline](#bioinformatics-pipeline)
2. [Data availability](#data-availability)
   * [Recommended files](#recommended-files)
   * [All available files](#all-available-files)
3. [Bioinformatics tools that use the UHGV](#code-availability) 
   * [Contig-level taxonomic classification](CLASSIFY.md)
   * [Read-level abundance profiling](#read-level-abundance-profiling-with-phanta)
   * [Genome visualization](#genome-visualization)
      

## Methods

### Data sources

We constructed the UHGV by integrating gut virome collections from a number of recent studies: 

- Metagenomic Gut Virus Compendium (MGV): https://doi.org/10.1038/s41564-021-00928-6
- Gut Phage Database (GPD): https://doi.org/10.1016/j.cell.2021.01.029
- Metagenomic Mobile Genetic Elements Database (mMGE): https://doi.org/10.1093/nar/gkaa869
- IMG Virus Resource v4 (IMG/VR): https://doi.org/10.1093/nar/gkac1037
- Hadza Hunter Gatherer Phage Catalog (Hadza): https://doi.org/10.1101/2022.03.30.486478
- Cenote Human Virome Database (CHVD): https://doi.org/10.1073/pnas.2023202118
- Human Virome Database (HuVirDB): https://doi.org/10.1016/j.chom.2019.08.008
- Gut Virome Database (GVD): https://doi.org/10.1016/j.chom.2020.08.003
- Atlas of Infant Gut DNA Virus Diversity (COPSAC): https://doi.org/10.1101/2021.07.02.450849
- Circular Gut Phages from NCBI (Benler et al.): https://doi.org/10.1186/s40168-021-01017-w
- Danish Enteric Virome Catalogue (DEVoC): https://doi.org/10.1128/mSystems.00382-21
- Stability of the human gut virome and effect of gluten-free diet (GFD): https://doi.org/10.1016/j.celrep.2021.109132

### Bioinformatics pipeline

Sequences from these studies were combined and run through the following bioinformatics pipeline:
- [geNomad](https://portal.nersc.gov/genomad/), [viralVerify](https://github.com/ablab/viralVerify), and [CheckV](https://bitbucket.org/berkeleylab/checkv) were used to remove sequences from cellular organisms and plasmids, as necessary
- [CheckV](https://bitbucket.org/berkeleylab/checkv) was used to trim remaining bacterial DNA from virus ends, estimate completeness, and identify closed genomes. Sequences >10Kb or >50% complete were retained and classified as either complete, high-quality (>90% complete), medium-quality (50-90% complete), or low-quality (<50% complete)
- [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) was used to calculate the average nucleotide identity between viruses using a [custom script](https://bitbucket.org/berkeleylab/checkv/src/master/scripts/anicalc.py)
- [DIAMOND](https://github.com/bbuchfink/diamond/) was used to blast proteins between viral genomes. Pairwise alignments were used to calculate a genome-wide protein-based similarity metric. 
- [MCL](http://micans.org) was used to cluster genomes into viral operational taxonomic units (vOTUs) at approximately the species, subgenus, genus, subfamily, and family-level ranks using a combination of genome-wide ANI for the species level and genome-wide proteomic similarity for higher ranks
- A representative genome was selected for each species level vOTU based on: presence of terminal repeats, completeness, and ratio of viral:non-viral genes
- [ICTV](https://ictv.global/vmr) taxonomy was inferred using a best-genome-hit approach to phage genomes from [INPHARED](https://github.com/RyanCook94/inphared) and using taxon-specific marker genes from [geNomad](https://portal.nersc.gov/genomad/) 
- [CRISPR](https://github.com/snayfach/MGV/tree/master/crispr_spacers) spacer matching and kmer matching with [PHIST](https://github.com/refresh-bio/PHIST) were used to connect viruses and host genomes. A voting procedure was used to then identify the host taxon at the lowest taxonomic rank comprising at least 70% of connections
- [HumGut](https://arken.nmbu.no/~larssn/humgut/) genomes and MAGs from a [Hadza](https://www.biorxiv.org/content/10.1101/2022.03.30.486478v2) hunter-gatherer population were used for host prediction and read mapping (HumGut contains all genomes from the [UHGG v1.0](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/) combined with NCBI genomes detected in gut metagenomes)
- [GTDB r207](https://gtdb.ecogenomic.org/) and [GTDB-tk](https://github.com/Ecogenomics/GTDBTk) were used to assign taxonomy to all prokaryotic genomes
- [BACPHLIP](https://github.com/adamhockenberry/bacphlip) was used for prediction of phage lifestyle together with integrases from the [PHROG ](https://phrogs.lmge.uca.fr/) database and prophage information from geNomad. Note: BACPHLIP tends to over classify viral genome fragments as lytic
- [Prodigal-gv](https://github.com/apcamargo/prodigal-gv) was used to identify protein-coding genes and alternative genetic codes
- [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper), [PHROGs](https://phrogs.lmge.uca.fr/), [KOfam](https://www.genome.jp/ftp/db/kofam/), [Pfam](http://pfam.xfam.org/), [UniRef_90](https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/), [PADLOC](https://github.com/padlocbio/padloc), and the [AcrCatalog](http://acrcatalog.pythonanywhere.com/) were used for phage gene functional annotation
- [PhaNNs](https://github.com/Adrian-Cantu/PhANNs) were used to infer phage structural genes
- [DGRscan](https://github.com/YuzhenYe/DGRscan) was used to identify diversity-generating retroelements on viruses containing reverse transcriptases
- [Bowtie2](https://github.com/BenLangmead/bowtie2) was used to align short reads from 1798 whole-metagenomes and 673 viral-enriched metagenomes against the UHGV and database of prokaryotic genomes. [ViromeQC](https://github.com/SegataLab/viromeqc) was used to select human gut viromes. [CoverM](https://github.com/wwood/CoverM) was used to estimate the breadth of coverage and we applied a 50% threshold for classifying virus presence-absence

For additional details, please refer to our manuscript: (in preparation).

## Data availability

The entire resource is freely available at: https://portal.nersc.gov/UHGV

We provide genomes for three quality tiers: 
- Full: >50% complete or >10Kbp, high-confidence & uncertain viral predictions
- Medium-quality: >50% complete, high-confidence viral predictions
- High-quality : >90% complete, high-confidence viral predictions 

Additionally, we provide data for:
- vOTU representatives
- All genomes in each vOTU

### Recommended files
For most analyses, we recommend using these files:
- [High-quality representative genomes](https://portal.nersc.gov/UHGV/genome_catalogs/votus_hq_plus.fna.gz)
- [Metadata for all species level vOTUs](https://portal.nersc.gov/UHGV/metadata/votus_full_metadata.tsv)

### All available files:

- metadata/

   - uhgv_full_metadata.tsv : detailed information on each of the 874,104 UHGV genome sequences
   - votus_full_metadata.tsv : detailed information on each of the 168,570 species level viral clusters
   - votus_metadata_extended.tsv: additional information on each vOTU
   - host_metadata.tsv : taxonomy and other info for prokaroytic genomes (completeness, contamination, n50)

- genome_catalogs/

   - uhgv_full.[fna|faa].gz : sequences for all genomes >10kb or >50% completeness 
   - uhgv_mq_plus.[fna|faa].gz : sequences for all genomes with >50% completeness 
   - uhgv_hq_plus.[fna|faa].gz : sequences for all genomes with >90% completeness 
   - votus_full.[fna|faa].gz : sequences for for vOTU representatives >10kb or >50% completeness
   - votus_mq_plus.[fna|faa].gz : sequences for for vOTU representatives with >50% completeness 
   - votus_hq_plus.[fna|faa].gz : sequences for vOTU representatives with >90% completeness 

- votu_reps/

   - [genome_id].fna : DNA sequence FASTA file of the genome assembly of the species representative
   - [genome_id].faa : protein sequence FASTA file of the species representative
   - [genome_id].gff : genome GFF file with various sequence annotations
   - [genome_id]_emapper.tsv : eggNOG-mapper annotations of the protein-coding sequences
   - [genome_id]_annotations.tsv : tab-delimited file containing diverse protein-coding annotations (PHROG, Pfam, UniRef90, eggNOG-mapper, PhANNs, KEGG)

- host_predictions/ 

   - crispr_spacers.fna : 5,318,089 CRISPR spacers from UHGG (3,143,456), NCBI (1,568,807), and Hadza genomes (605,826)
   - host_genomes_info.tsv : GTDB r207 taxonomy for genomes from the UHGG (286,387), NCBI (123,500), and Hadza genomes (54,779)
   - host_assignment_crispr.tsv : detailed information for host prediction with CRISPR spacers
   - host_assignment_kmers.tsv : detailed information for host prediction with PHIST kmer matching

- annotations/
   - functional annotation matrices: vOTUs x functions (PHROG, Pfam, KOfam, PADLOC)

- read_mapping/ 

   - metagenomes_prok_vir_counts_matrix.tsv.gz : coverM mapping statistics for viruses and bacteria across bulk metagenomes 
   - viromes_prok_vir_counts_matrix.tsv.gz : coverM mapping statistics for viruses and bacteria across viral-enriched metagenomes 
   - sample_metadata.tsv: human sample metadata (country, lifestyle, age, gender, bmi, study)
   - fastq_summary.tsv: information on sequencing reads (sra, bulk/virome metagenome, viromeQC enrichment, read counts)
   - study_metadata.tsv: information on individual studies for read mapping
   
   - bowtie2_indexes/
     - prokaryote_reps.fna.gz: FASTA of prokaryotic genomes used for read mapping
     - prokaryote_metadata_table.tsv.gz: prok genome metadata
     - prokaryote_reps.1.bt*: bowtie2 indexes

## Code availability

### Contig-level taxonomic classification with the UHGV toolkit

- Code to assign viral genomes to taxonomic groups from the UHGV
- View the [README](CLASSIFY.md) for download and usage instructions.

### Read-level abundance profiling with Phanta

- Phanta (https://github.com/bhattlab/phanta) is a fast and accurate virus-inclusive profiler of human gut metagenomes based on the classification of short reads with Kraken2. 
- Follow the instructions to install the software at the [Phanta Github page](https://github.com/bhattlab/phanta#quick-start)
- Download a custom-built UHGV database for Phanta:
  - HQ plus: `wget http://ab_phanta.os.scg.stanford.edu/Phanta_DBs/humgut_uhgv_hqplus_v1.tar.gz`
  - MQ plus: `wget http://ab_phanta.os.scg.stanford.edu/Phanta_DBs/humgut_uhgv_mqplus_v1.tar.gz`
  - These databases are similar to Phanta's default database as described in Phanta's manuscript but replacing the viral portion of Phanta’s default DB with UHGV.
- Phanta can be executed based on the instructions on its GitHub page.


### Genome visualization

- Species level genomes can be visualized using [Geneious](https://www.geneious.com/) or other tools that accept GFF3 format.
- Example:
   * Identify a species of interest: UHGV-0014815
   * Download a GFF file for species of interest: https://portal.nersc.gov/UHGV/votu_reps/UHGV-001/UHGV-0014815/UHGV-0014815.gff)
   * Geneious > Import GFF
   * Menu > Sequence > Circularize
