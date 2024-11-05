#! /bin/bash

set -e

## Sites found in plasmids

python make_alignments.py \
    --att-sites results/att_sites.json \
    --output-dir results/alignment \

python make_consensus_sites.py \
    --input-dir results/alignment \
    --output-file results/consensus_sites.tsv

## Sites found in plasmids + combinatorial sites

python make_alignments.py \
    --att-sites results/att_sites_combinatorial.json \
    --output-dir results/alignment_combinatorial \

python make_consensus_sites.py \
    --input-dir results/alignment_combinatorial \
    --output-file results/consensus_sites_combinatorial.tsv

## Make consensus alignment

python make_consensus_alignment.py