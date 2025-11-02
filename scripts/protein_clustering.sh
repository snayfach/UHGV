#!/usr/bin/env bash

mmseqs cluster -s 7.0 -c 0.8 --cov-mode 0 -e 1e-3 --cluster-mode 0 --cluster-steps 3 \
  --cluster-reassign 1 --kmer-per-seq 50 uhgv_db/uhgv_db cluster_db/cluster_db tmp
rm -rf tmp
mmseqs result2msa uhgv_db/uhgv_db uhgv_db/uhgv_db cluster_db/cluster_db msa_db/msa_db \
  --msa-format-mode 2
mmseqs msa2profile --msa-type 2 --match-mode 1 --match-ratio 0.5 msa_db/msa_db profile_db/profile_db
mmseqs profile2consensus profile_db/profile_db consensus_db/consensus_db
mmseqs search --cov-mode 0 -c 0.9 -s 8.0 -e 1e-4 --add-self-matches 1 -a 1 profile_db/profile_db \
  consensus_db/consensus_db profile_consensus_search_db/profile_consensus_search_db tmp
rm -rf tmp
mmseqs clust profile_db/profile_db profile_consensus_search_db/profile_consensus_search_db \
  profile_cluster_db/profile_cluster_db
mmseqs mergeclusters uhgv_db/uhgv_db merged_cluster_db/merged_cluster_db cluster_db/cluster_db \
  profile_cluster_db/profile_cluster_db
mmseqs createtsv uhgv_db/uhgv_db uhgv_db/uhgv_db merged_cluster_db/merged_cluster_db  \
  uhgv_protein_clusters.tsv
