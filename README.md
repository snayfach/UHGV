# Unified Human Gut Virome Catalog (UHGV)

We constructed the UHGV by integrated gut virome collections from a number of recent studies: 

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

## Bioinformatics pipeline

Viral genomes were downloaded from the sources listed above. [geNomad](https://portal.nersc.gov/genomad/) was used to confirm a viral origin and to excise provirues from bacterial chromosomes, as necessary. [CheckV](https://bitbucket.org/berkeleylab/checkv) was used to trim any remaining bacterial DNA from virus ends, estimate completeness, and identify closed genomes. Sequences >10Kb or >50% complete were retained and classified as either complete, high-quality (>90% complete), medium-quality (50-90% complete), or low-quality (<50% complete). Viral genomes were clustered into operational taxonomic units at approximately the species, subgenus, genus, subfamily, and family-level ranks using a combination of genome-wide ANI (for species level vOTUs) and genome-wide proteomic similarity (for higher ranks). Taxonomy was inferred using a best-hit approach to phage genomes from the ICTV and using taxon-specific marker genes from geNomad. Hosts were predicted using a large database of bacterial and archael genomes that included the [UHGG v2](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/) and [>48K MAGs from Hadza hunter gatherers](https://www.biorxiv.org/content/10.1101/2022.03.30.486478v2). Prediction methods included [CRISPR spacer matching](https://github.com/snayfach/MGV/tree/master/crispr_spacers) and kmer matching with [PHIST](https://github.com/refresh-bio/PHIST). Phage lifestyle was inferred using [BACPHLIP](https://github.com/adamhockenberry/bacphlip), integrases from the [PHROG database](https://phrogs.lmge.uca.fr/), and excision information from geNomad. Protein coding genes and alternative genetic codes were predicted using [Prodigal-gv](https://github.com/apcamargo/prodigal-gv). For additional detail, please refer to our manuscript (in preparation).

## Data Resource

The data resource can be accessed without restrictions at: https://portal.nersc.gov/cfs/m342/UHGV/

Details on individual files are listed below


## Compare new viruses to the UHGV

## Assembly-free detection of UHGV viral taxa from metagenomic reads


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

* All MAGs are estimated to be >=50% complete and <=10% contaminated

## Integrated genomes from the gut and other environments dataset (IGG)

Note: this set includes non-gut genomes

The 60,664 genomes from the HGM dataset were integrated together with 145,917 reference genomes from PATRIC and IMG, which include 16,525 publicly available MAGs from other studies <b> as well as genomes from other non-gut environments </b>. All 206,581 genomes met the MIMAG medium quality draft genome standard of >=50% completeness and <=10% contamination. Genomes were clustered into 23,790 species-level OTUs based on 95% genome-wide average nucleotide identity. 

<b>Representative genomes for all species (N=23,790)</b>  
download from your browser: [download link](https://bit.ly/IGG_all_23790_fna)  
download via wget: `wget http://bit.ly/IGG_all_23790_fna -O IGG_v1.0_all_23790_fna.tar.bz2`

* The representative genome is the highest quality genome that is most closely related to other genomes in the same species
* The selection was based on: completeness, contamination, N50, and % DNA identity  and % DNA alignment to other genomes in the same species
* High-quality genomes are picked as representatives when available

<b>Representative genomes for species with a high-quality genome (N=16,136)</b>  
download via from your browser: [download link](https://bit.ly/IGG_hq_16136_fna)  
download via wget: `wget http://bit.ly/IGG_hq_16136_fna -O IGG_v1.0_hq_16136_fna.tar.bz2`

* This is a more conservative, higher-quality set of representative genomes 

<b>Representative genomes for human gut species (N=4,558)</b>  
download from your browser: [download link](https://bit.ly/IGG_gut_4558_fna)  
download via wget: `wget http://bit.ly/IGG_gut_4558_fna -O IGG_v1.0_gut_4558_fna.tar.bz2`  

* Gut species were defined on the basis of 3 criteria:  
	1) Contain a genome from an organism isolated from the gut (N=955), or  
	2) Contain a MAG from the HGM dataset (N=2,962), or  
	3) Were detected in a human gut metagenome based on read-mapping to conserved species-specific marker genes using [IGGsearch](https://github.com/snayfach/IGGsearch) (N=4,347)
* 2,058 (45%) gut species are new while 2,500 (55%) contain a reference genome

<b>Representative genomes for human gut species with a high-quality genome (N=2,935)</b>  
download from your browser: [download link](https://bit.ly/IGG_gut_2935_fna)  
download via wget: `wget http://bit.ly/IGG_gut_2935_fna -O IGG_v1.0_gut_2935_fna.tar.bz2` 

* 684 (23%) gut species are new while 2,251 (77%) contain a reference genome


## Phylogenome trees

Phylogenetic trees were constructed for all Bacterial and Archaeal species in the IGGdb using concatenated alignments of conserved, single-copy marker gene families from the PhyEco database (N=88 for Bacteria and 100 for Archaea). Protein-based multiple sequence alignment was performed using FAMSA v1.2.5, which is designed for fast and accurate alignment of thousands of sequences. Gene alignments were concated together, columns with >15% gaps were dropped, and seuqneces with >70% gaps and were removed (N=39). FastTree2 was used to build a maximum likelihood phylogeny.

<b>All Bacterial species (N=22,515)</b>  
[download alignments](https://bit.ly/IGG_bact_22515_msa)  
[download tree](https://bit.ly/IGG_bact_22515_tre)

<b>All Archaeal species (N=1,236)</b>  
[download alignments](https://bit.ly/IGG_arch_1236_msa)   
[download tree](https://bit.ly/IGG_arch_1236_tre)


## Metadata for the HGM and IGG datasets

<b>Reference genomes and MAGs from HGM dataset (N=206,581)</b>  
[download link](https://bit.ly/IGG_genome_info_206581)

* Contains rich metadata on all reference genomes and MAGs, including quality information and habitat
* Completeness and contamination estimates were made using [CheckM](https://github.com/Ecogenomics/CheckM)

<b>All species from the IGGdb (N=23,790)</b>  
[download link](https://bit.ly/IGG_species_info_23790)  

* Contains rich metadata on each species, including taxonomic information
* Species were taxonomically annotated based on the [Genome Taxonomy Database](https://github.com/Ecogenomics/GTDBTk)
* Additionaly, species were assigned to an operational taxonomy based on phylogenetic clustering using rank-specific distance cutoffs
