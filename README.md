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
- Viral genomes were clustered into operational taxonomic units at approximately the species, subgenus, genus, subfamily, and family-level ranks using a combination of genome-wide ANI (for species level vOTUs) and genome-wide proteomic similarity (for higher ranks). 
- A representative genome was selected for each vOTU based on: terminal repeats, completeness, and viral genes
- Taxonomy was inferred using a best-hit approach to phage genomes from the ICTV and using taxon-specific marker genes from geNomad. 
- Prokaryotic hosts were predicted using a large database of bacterial and archael genomes that included the [UHGG v2](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/) and [>48K MAGs from Hadza hunter gatherers](https://www.biorxiv.org/content/10.1101/2022.03.30.486478v2). Prediction methods included [CRISPR spacer matching](https://github.com/snayfach/MGV/tree/master/crispr_spacers) and kmer matching with [PHIST](https://github.com/refresh-bio/PHIST)
- Phage lifestyle was inferred using [BACPHLIP](https://github.com/adamhockenberry/bacphlip), integrases from the [PHROG database](https://phrogs.lmge.uca.fr/), and excision information from geNomad
- Protein coding genes and alternative genetic codes were predicted using [Prodigal-gv](https://github.com/apcamargo/prodigal-gv)
- Proteins were functionally annotated using the eggNOG-mapper, PHROGs, KOfam, Pfam, and UniProt

For additional details, please refer to our manuscript: (in preparation).

## The following files are available for download


**- metadata/**

- uhgv_full_884377.tsv : detailed information on each of the 884,377 UHGV genome sequences
- votus_full_171338.tsv : detailed information on each of the 171338 species level viral clusters

**- genomes/** nucleotide sequences for different subsets of the database

- uhgv_hq+_212886.fna.gz : all genomes with >90% completeness 
- uhgv_mq+_434701.fna.gz : all genomes with >50% completeness  
- uhgv_full_884377.fna.gz : all genomes with >50% completeness or >10Kb
- votus_hq+_58864.fna.gz : representative genomes for vOTUs with >=1 genome w/ >90% completeness
- votus_mq+_102355.fna.gz : representative genomes for vOTUs with >=1 genome w/ >50% completeness
- votus_full_171338.fna.gz : representative genomes for vOTUs with >=1 genome w/ >50% completeness or >10Kb
 
**- host_prediction/**



<b>High quality MAGs (N=24345)</b>   
download from your browser: [download link](https://bit.ly/HGM_hq_24345_fna)  
download via wget: `wget http://bit.ly/HGM_hq_24345_fna -O HGM_v1.0_hq_24345_fna.tar.bz2`  

* High quality MAGs meet the following criterea:
	*  \>=90% estimated completeness
	*  <=5% estimated contamination
	*  <=500 contigs
	*  \>=5kb average contig length
	*  \>=10kb contig N50
	*  \>=5x read depth over >=90% of contigs

<b>High and medium quality MAGs (N=60664)</b>   
download from your browser: [download link](https://bit.ly/HGM_all_60664_fna)  
download via wget: `wget http://bit.ly/HGM_all_60664_fna -O HGM_v1.0_all_60664_fna.tar.bz2` 


## Code to compare new viruses to the UHGV

## Code to perform assembly-free detection of UHGV viral taxa
