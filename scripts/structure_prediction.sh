#!/usr/bin/env bash

mmseqs createseqfiledb --min-sequences 20 uhgv_db/uhgv_db merged_cluster_db/merged_cluster_db \
  unaligned_db/unaligned_db
mmseqs apply unaligned_db/unaligned_db aligned_db/aligned_db -- famsa STDIN STDOUT -t 1
mmseqs unpackdb --unpack-suffix ".faa" --unpack-name-mode 0 aligned_db/aligned_db aligned_unpacked_db
ls -1 aligned_unpacked_db | rush --eta 'hhconsensus -v 0 -M first -id 97 -i aligned_unpacked_db/{} \
  -o aligned_unpacked_a3m_db/{.}.a3m'
ls -1 aligned_unpacked_a3m_db | rush --eta 'awk -v FILEN=$(basename -s .a3m {}) '\''{if (NR==1) {next} else {print}}'\'' aligned_unpacked_a3m_db/{} | sponge aligned_unpacked_a3m_db/{}'
ls -1 aligned_unpacked_a3m_db | rush --eta '[ $(rg -c "^>" aligned_unpacked_a3m_db/{}) -lt 20 ] && rm aligned_unpacked_a3m_db/{} || :'
mv aligned_unpacked_a3m_db MSAs/UHGV

mkdir uhgv_hhsuite_db
ffindex_build -s uhgv_hhsuite_db/uhgv_hhsuite_db_a3m.ff{data,index} MSAs/UHGV/*.a3m
ffindex_apply uhgv_hhsuite_db/uhgv_hhsuite_db_a3m.ff{data,index} \
  -i uhgv_hhsuite_db/uhgv_hhsuite_db_hhm.ffindex \
  -d uhgv_hhsuite_db/uhgv_hhsuite_db_hhm.ffdata -- hhmake -i stdin -o stdout -v 0
cstranslate -f -x 0.3 -c 4 -I a3m -i uhgv_hhsuite_db/uhgv_hhsuite_db_a3m \
  -o uhgv_hhsuite_db/uhgv_hhsuite_db_cs219

mkdir -p enriched/UHGV
fd a3m MSAs/UHGV | rush -j 12 --eta './run_hhblits.sh {}'