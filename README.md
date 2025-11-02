# Unified Human Gut Virome Catalog (UHGV)

The **UHGV** is a comprehensive genomic resource of viruses from the human gut microbiome. Genomes were derived from [12 independent data sources](#data-sources) and annotated using a [uniform bioinformatics pipeline](#bioinformatics-pipeline):


<img src="img/data_workflow.png" width="900">


## Table of Contents
1. [Methods](#methods)
   - [Data sources](#data-sources)
   - [Bioinformatics pipeline](#bioinformatics-pipeline)
2. [Data Availability](#data-availability)
   - [Recommended files](#recommended-files)
   - [All available files](#all-available-files)
3. [Tools Using the UHGV](#code-availability)
   - [Genome Taxonomy Classification](#genome-taxonomy-classification)
   - [Read-Level Abundance Profiling](#read-level-abundance-profiling)
   - [Genome Visualization](#genome-visualization)
4. [Citation](#citation)


## Methods

### Data sources

The UHGV integrates gut virome collections from recent studies:

1. [Metagenomic Gut Virus Compendium (MGV)](https://doi.org/10.1038/s41564-021-00928-6)
2. [Gut Phage Database (GPD)](https://doi.org/10.1016/j.cell.2021.01.029)
3. [Metagenomic Mobile Genetic Elements Database (mMGE)](https://doi.org/10.1093/nar/gkaa869)
4. [IMG Virus Resource v4 (IMG/VR)](https://doi.org/10.1093/nar/gkac1037)
5. [Hadza Hunter Gatherer Phage Catalog (Hadza)](https://doi.org/10.1101/2022.03.30.486478)
6. [Cenote Human Virome Database (CHVD)](https://doi.org/10.1073/pnas.2023202118)
7. [Human Virome Database (HuVirDB)](https://doi.org/10.1016/j.chom.2019.08.008)
8. [Gut Virome Database (GVD)](https://doi.org/10.1016/j.chom.2020.08.003)
9. [Atlas of Infant Gut DNA Virus Diversity (COPSAC)](https://doi.org/10.1101/2021.07.02.450849)
10. [Circular Gut Phages from NCBI (Benler et al.)](https://doi.org/10.1186/s40168-021-01017-w)
11. [Danish Enteric Virome Catalogue (DEVoC)](https://doi.org/10.1128/mSystems.00382-21)
12. [Stability of the human gut virome and effect of gluten-free diet (GFD)](https://doi.org/10.1016/j.celrep.2021.109132)


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

The UHGV resource is freely available at: [https://portal.nersc.gov/UHGV](https://portal.nersc.gov/UHGV)

We provide genomes at three quality tiers:

| Tier | Criteria |
|------|----------|
| Full | >50% complete or >10 Kbp; high-confidence & uncertain viral predictions |
| Medium-quality | >50% complete; high-confidence viral predictions |
| High-quality | >90% complete; high-confidence viral predictions |

These data are provided for either vOTU representatives or all genomes in each vOTU.

### Recommended files

| File | Description | Link |
|------|-------------|------|
| `votus_hq_plus.fna.gz` | High-quality representative genomes | [Download](https://portal.nersc.gov/UHGV/genome_catalogs/votus_hq_plus.fna.gz) |
| `votus_metadata.tsv` | Metadata for all species-level vOTUs | [Download](https://portal.nersc.gov/UHGV/metadata/votus_metadata.tsv) |


### All available files:

**metadata/**
- `uhgv_metadata.tsv`: information for each of the 873,995 UHGV genomes
- `votus_metadata.tsv`: information for 168,536 species-level viral clusters
- `votus_metadata_extended.tsv`: additional vOTU details
- `host_metadata.tsv`: taxonomy, completeness, contamination, N50 for prokaryotic genomes
- `source_biosample_metadata.tsv`: information for the samples from which virus genomes were obtained

**genome_catalogs/**
- `uhgv_full.[fna|faa].gz`: all genomes >10 kb or >50% complete
- `uhgv_mq_plus.[fna|faa].gz`: genomes >50% complete
- `uhgv_hq_plus.[fna|faa].gz`: genomes >90% complete
- `votus_full.[fna|faa].gz`: vOTU representatives >10 kb or >50% complete
- `votus_mq_plus.[fna|faa].gz`: vOTU representatives >50% complete
- `votus_hq_plus.[fna|faa].gz`: vOTU representatives >90% complete
- `prokaryote_reps.fna.gz`: genomic sequences of gut prokaryotes

**votu_reps/**
- `[genome_id].fna`: DNA sequence
- `[genome_id].faa`: protein sequence
- `[genome_id].gff`: genome annotations
- `[genome_id]_emapper.tsv`: eggNOG-mapper annotations
- `[genome_id]_annotations.tsv`: PHROG, Pfam, UniRef90, eggNOG-mapper, PhANNs, KEGG annotations

> Only available for genomes >50% complete with confident virus prediction

**host_predictions/**
- `crispr_spacers.fna`: 5,318,089 CRISPR spacers
- `host_genomes_info.tsv`: GTDB r207 taxonomy for UHGG, NCBI, Hadza genomes
- `host_assignment_crispr.tsv`: host predictions via CRISPR
- `host_assignment_kmers.tsv`: host predictions via PHIST

**annotations/**
- Functional annotation matrices (vOTUs × functions: PHROG, Pfam, KOfam, PADLOC)

**read_mapping/**
- `metagenomes_prok_vir_counts_matrix.tsv.gz`: CoverM statistics for bulk metagenomes
- `viromes_prok_vir_counts_matrix.tsv.gz`: CoverM statistics for viral-enriched metagenomes
- `relative_abundance.tsv`: Per-sample relative abundances of viruses and hosts derived from read mapping data
- `sample_metadata.tsv`: sample metadata (country, lifestyle, age, gender, BMI, study)
- `fastq_summary.tsv`: sequencing reads info
- `study_metadata.tsv`: per-study metadata

- **bowtie2_indexes/**
  - `prokaryote_reps.fna.gz`: prokaryotic genome FASTA
  - `prokaryote_metadata_table.tsv.gz`: prok genome metadata
  - `prokaryote_reps.1.bt*`: Bowtie2 indexes

## Code availability

### Genome Taxonomy Classification
**UHGV-classifier**: command-line tool for classifying genomes using UHGV.
- [GitHub & installation](https://github.com/snayfach/UHGV-classifier)

### Read-level Abundance Profiling
**Phanta**: virus-inclusive profiler for human gut metagenomes.
- [GitHub & installation](https://github.com/bhattlab/phanta#quick-start)
- UHGV databases:
  - HQ plus: `wget http://ab_phanta.os.scg.stanford.edu/Phanta_DBs/humgut_uhgv_hqplus_v1.tar.gz`
  - MQ plus: `wget http://ab_phanta.os.scg.stanford.edu/Phanta_DBs/humgut_uhgv_mqplus_v1.tar.gz`
> Databases replace the viral portion of Phanta’s default DB with UHGV sequences.

### Genome Visualization
- Use [Geneious](https://www.geneious.com/) or any GFF3-compatible tool.
- Example workflow for a species (`UHGV-0014815`):
  1. Download GFF: `https://portal.nersc.gov/UHGV/votu_reps/UHGV-001/UHGV-0014815/UHGV-0014815.gff`
  2. Import into Geneious
  3. Menu → Sequence → Circularize
> Can also be applied with other GFF3 visualization software.

## Citation

If you use the UHGV in your research, please cite both the database and the underlying publication:

**Data resource:**  
Nayfach, S., & Camargo, A. (2025). Unified Human Gut Virome (UHGV) (1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17402089

**Publication:**
Coming soon
