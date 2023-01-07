# Unified Human Gut Virome Catalog (UHGV)

The UHGV is a comprehensive genomic resource of viruses from the human microbiome. Genomes were derived from [12 independent data sources](#data-sources) and annotated using a [uniform bioinformatics pipeline](#bioinformatics-pipeline):

<img src="data_workflow.png" width="900">


## Table of contents
1. [Methods](#methods)
   * [Data sources](#data-sources)
   * [Bioinformatics pipeline](#bioinformatics-pipeline)
2. [Data availability](#data-availability)
   * [Recommended files](#recommended-files)
   * [All available files](#all-available-files)
3. [Code availability](#code-availability) 
   * [Phylogenetic placement](#phylogenetic-placement)
   * [Read mapping](#read-mapping)
      * [Using Bowtie2](#bowtie2)
      * [Using Phanta](#phanta)
   * [Genome visualization](#genome-visualization)
4. [Keeping the UHGV up-to-date](#updating-the-uhgv)
      

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
- [geNomad](https://portal.nersc.gov/genomad/) was used to confirm a viral origin and to excise provirues from bacterial chromosomes (as necessary)
- [CheckV](https://bitbucket.org/berkeleylab/checkv) was used to trim remaining bacterial DNA from virus ends, estimate completeness, and identify closed genomes. Sequences >10Kb or >50% complete were retained and classified as either complete, high-quality (>90% complete), medium-quality (50-90% complete), or low-quality (<50% complete)
- [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) was used to calculate the average nucleotide identity between viruses using a [custom script](https://bitbucket.org/berkeleylab/checkv/src/master/scripts/anicalc.py)
- Viral genomes were clustered into viral operational taxonomic units (vOTUs) at approximately the species, subgenus, genus, subfamily, and family-level ranks using a combination of genome-wide ANI for the species level and genome-wide proteomic similarity for higher ranks
- A representative genome was selected for each species level vOTU based on: presence of terminal repeats, completeness, and ratio of viral:non-viral genes
- [ICTV](https://ictv.global/vmr) taxonomy was inferred using a best-genome-hit approach to phage genomes from [INPHARED](https://github.com/RyanCook94/inphared) and using taxon-specific marker genes from [geNomad](https://portal.nersc.gov/genomad/) 
- [CRISPR](https://github.com/snayfach/MGV/tree/master/crispr_spacers) spacer matching and kmer matching with [PHIST](https://github.com/refresh-bio/PHIST) were used to connect viruses and host genomes. A voting procedure was used to then identify the host taxon at the lowest taxonomic rank comprising at least 70% of connections
- [HumGut](https://arken.nmbu.no/~larssn/humgut/) genomes and MAGs from a [Hadza](https://www.biorxiv.org/content/10.1101/2022.03.30.486478v2) hunter gatherer population were used for host prediction and read mapping (HumGut contains all genomes from the [UHGG v1.0](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/) combined with NCBI genomes detected in gut metagenomes)
- [GTDB r207](https://gtdb.ecogenomic.org/) and [GTDB-tk](https://github.com/Ecogenomics/GTDBTk) were used to assign taxonomy to all prokaryotic genomes
- [Placeholder] was used to estimate the host range of each virus based on the average phylogenetic distance between connected host genomes
- [BACPHLIP](https://github.com/adamhockenberry/bacphlip) was used for prediction of phage lifestyle together with integrases from the [PHROG ](https://phrogs.lmge.uca.fr/) database and prophage information from geNomad. Note: BACPHLIP tends to over classify viral genome fragments as lytic
- [Prodigal-gv](https://github.com/apcamargo/prodigal-gv) was used to identify protein coding genes and alternative genetic codes
- [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper), [PHROGs](https://phrogs.lmge.uca.fr/), [KOfam](https://www.genome.jp/ftp/db/kofam/), [Pfam](http://pfam.xfam.org/), [UniRef_90](https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/), [PADLOC](https://github.com/padlocbio/padloc), and the [AcrCatalog](http://acrcatalog.pythonanywhere.com/) were used for phage gene functional annotation
- [PhaNNs](https://github.com/Adrian-Cantu/PhANNs) were used to infer phage structural genes
- [DGRscan](https://github.com/YuzhenYe/DGRscan) was used to identify diversity generating retroelements on viruses containing reverse transcriptases
- [Bowtie2](https://github.com/BenLangmead/bowtie2) was used to align short reads from 1798 whole-metagenomes and 673 viral-enriched metagenomes against the UHGV and database of prokaryotic genomes. [ViromeQC](https://github.com/SegataLab/viromeqc) was used to select human gut viromes. [CoverM](https://github.com/wwood/CoverM) was used to estimate breadth of coverage and we applied a 50% threshold for classifing virus presence-absence

For additional details, please refer to our manuscript: (in preparation).

## Data availability

The entire resource is freely available at: https://portal.nersc.gov/UHGV

#### Recommended files:

- [Representative genomes with >50% completeness](https://portal.nersc.gov/UHGV/genome_catalogs/votus_mq_plus.fna.gz)
- [Metadata for all species level vOTUs](https://portal.nersc.gov/UHGV/metadata/votus_full_metadata.tsv)

#### All available files:

- metadata/

   - uhgv_full_metadata.tsv : detailed information on each of the 884,377 UHGV genome sequences
   - votus_full_metadata.tsv : detailed information on each of the 171,338 species level viral clusters

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
   - [genome_id]_emapper.tsv : eggNOG-mapper annotations of the protein coding sequences
   - [genome_id]_annotations.tsv : tab-delimited file containg diverse protein coding annotations (PHROG, Pfam, UniRef90, eggNOG-mapper, PhANNs, KEGG)

- host_predictions/ 

   - crispr_spacers.fna : 5,318,089 CRISPR spacers from UHGG (3,143,456), NCBI (1,568,807), and Hadza genomes (605,826)
   - host_genomes_info.tsv : GTDB r207 taxonomy for genomes from the UHGG (286,387), NCBI (123,500), and Hadza genomes (54,779)
   - host_assignment_crispr.tsv : detailed information for host prediction with CRISPR spacers
   - host_assignment_kmers.tsv : detailed information for host prediction with PHIST kmer matching

- read_mapping/ 

   - metagenomes_prok_vir_counts_matrix.tsv.gz : mapping statistics for viruses and bacteria across bulk metagenomes 
   - viromes_prok_vir_counts_matrix.tsv.gz : mapping statistics for viruses and bacteria across viral-enriched metagenomes 
   - study_sample_metadata.xlsx : information on bulk metagenomes and viral-enriched metagenomes  
   - viral_bowtie2_index : 
   - bacterial_bowtie2_index : 

## Code availability

### Phylogenetic placement

Available at https://github.com/snayfach/classiPhi

Use cases:

- Determine novelty of new genome: does my viral genome represent a novel species? genus? family?
- Ecological analysis: compare viral phylogenetic groups across samples
- Comparative genomics: retrieve other viral genomes from the same phylogenetic group
- Infer host and lifestyle: impute characteristics of the virus based on it's phylogenetic group
- Update the UHGV: cluster unclassified viral genomes into de novo vOTUs

### Read-mapping 

#### Bowtie2

(in preparation)

#### Phanta

Phanta uses Kraken2 to efficiently quantity the presence of viruses and prokyotes
(in preparation)

### Genome visualization

Species level genomes can be visualized using [Geneious](https://www.geneious.com/) or other tools that accept GFF3 format.

Example:
* Identify a species of interest: UHGV-0014815
* Download a GFF file for species of interest: https://portal.nersc.gov/UHGV/votu_reps/UHGV-001/UHGV-0014815/UHGV-0014815.gff)
* Geneious > Import GFF
* Menu > Sequence > Circularize


## Updating the UHGV

The human gut virome harbors immense diversity that may not be fully captured by the UHGV, particularly below the genus rank, or for understudied human populations.

Please email snayfach@gmail.com if you'd like to include your sequences in the next version of the data resource. Please provide the following:
- FASTA file of predicted viruses
- INDSC accession numbers (SRA/ENA/DDBJ) corresponding to sequencing reads used to assemble viruses
- Tool used for bioinformatic virus prediction
- Whether the sample was derived from bulk or viral metagenome
- Brief description of the study
- Publication DOI

If your sequences are unpublished, you can update the UHGV using [classiPhi](https://github.com/snayfach/ClassiPhi). The pipeline will compare your viruses against the UHGV, classify them to their appropriate rank, and cluster unclassified viruses into **de novo** viral clusters.
