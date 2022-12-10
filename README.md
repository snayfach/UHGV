# Unified Human Gut Virome Catalog (UHGV)

Catalog of human gut viruses integrated from 12 databases

(Manuscript in preparation)

## Data availability

The entire resource is freely available at our HTTP site: https://portal.nersc.gov/cfs/m342/UHGV/

## Methods

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

Sequences from these studies were combined and run through the following bioinformatics pipeline:
- [geNomad](https://portal.nersc.gov/genomad/) was used to confirm a viral origin and to excise provirues from bacterial chromosomes (as necessary)
- [CheckV](https://bitbucket.org/berkeleylab/checkv) was used to trim remaining bacterial DNA from virus ends, estimate completeness, and identify closed genomes. Sequences >10Kb or >50% complete were retained and classified as either complete, high-quality (>90% complete), medium-quality (50-90% complete), or low-quality (<50% complete)
- Viral genomes were clustered into viral operational taxonomic units (vOTUs) at approximately the species, subgenus, genus, subfamily, and family-level ranks using a combination of genome-wide ANI for the species level and genome-wide proteomic similarity for higher ranks
- A representative genome was selected for each species level vOTU based on: presence of terminal repeats, completeness, and ratio of viral:non-viral genes
- [ICTV](https://ictv.global/vmr) taxonomy was inferred using a best-genome-hit approach to phage genomes from [INPHARED](https://github.com/RyanCook94/inphared) and using taxon-specific marker genes from [geNomad](https://portal.nersc.gov/genomad/) 
- [CRISPR](https://github.com/snayfach/MGV/tree/master/crispr_spacers) spacer matching and kmer matching with [PHIST](https://github.com/refresh-bio/PHIST) were used to connect viruses and host genomes
- [HumGut](https://arken.nmbu.no/~larssn/humgut/) genomes and MAGs from a [Hadza](https://www.biorxiv.org/content/10.1101/2022.03.30.486478v2) hunter gatherer population were used for host prediction and read mapping. Note: HumGut contains all genomes from the [UHGG v1.0](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/) combined with NCBI genomes detected in gut metagenomes
- [GTDB r207](https://gtdb.ecogenomic.org/) and [GTDB-tk](https://github.com/Ecogenomics/GTDBTk) were used to assign taxonomy to all prokaryotic genomes
- [BACPHLIP](https://github.com/adamhockenberry/bacphlip) was used for prediction of phage lifestyle together with integrases from the [PHROG ](https://phrogs.lmge.uca.fr/) database and prophage information from geNomad
- [Prodigal-gv](https://github.com/apcamargo/prodigal-gv) was used to identify protein coding genes and alternative genetic codes
- [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper), [PHROGs](https://phrogs.lmge.uca.fr/), [KOfam](https://www.genome.jp/ftp/db/kofam/), [Pfam](http://pfam.xfam.org/), and [UniRef_90](https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/) were used for phage gene functional annotation

For additional details, please refer to our manuscript: (in preparation).

## The following files are available for download


**- metadata/**

- uhgv_full_884377.tsv : detailed information on each of the 884,377 UHGV genome sequences
- votus_full_171338.tsv : detailed information on each of the 171338 species level viral clusters

(in preparation) cluster info 

**- genomes/** nucleotide sequences for different subsets of the database

- uhgv_hq+_212886.fna.gz : all genomes with >90% completeness 
- uhgv_mq+_434701.fna.gz : all genomes with >50% completeness  
- uhgv_full_884377.fna.gz : all genomes with >50% completeness or >10Kb
- votus_hq+_58864.fna.gz : representative genomes for vOTUs with >=1 genome w/ >90% completeness
- votus_mq+_102355.fna.gz : representative genomes for vOTUs with >=1 genome w/ >50% completeness
- votus_full_171338.fna.gz : representative genomes for vOTUs with >=1 genome w/ >50% completeness or >10Kb

**- genes/**

(in preparation) will include GFF file for each species level representative 

**- host_prediction/**

- crispr_spacers.fna : 5,318,089 CRISPR spacers from UHGG (3,143,456), NCBI (1,568,807), and Hadza genomes (605,826)
- host_genomes_info.tsv : GTDB r207 taxonomy for genomes from the UHGG (286,387), NCBI (123,500), and Hadza genomes (54,779)
- host_assignment_crispr.tsv : detailed information for host prediction with CRISPR spacers
- host_assignment_kmers.tsv : detailed information for host prediction with PHIST kmer matching


## Code to compare new viruses to the UHGV

(in preparation)

## Code to perform assembly-free detection of UHGV viral taxa

(in preparation)
