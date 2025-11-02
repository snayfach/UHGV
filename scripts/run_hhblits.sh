#!/usr/bin/env bash

INPUT=${1}
OUTPUT=$(echo ${1} | sed 's/MSAs/enriched/g')

hhblits \
  -d UniRef30_2023_02/UniRef30_2023_02 \
  -d uhgv_hhsuite_db/uhgv_hhsuite_db \
  -i $INPUT -o /tmp/tmp.hhr -oa3m $OUTPUT -id 92.5 \
  -n 3 -cpu 4 -e 0.001 -v 0
