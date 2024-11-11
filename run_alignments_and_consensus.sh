#! /bin/bash

set -e

## Sites found in plasmids
echo -e "\033[32mMaking alignments for plasmid sites\033[0m"
python make_alignments.py \
    --att-sites results/att_sites.json \
    --output-dir results/alignment \

echo -e "\033[32mMaking consensus sites for plasmid sites\033[0m"
python make_consensus_sites.py \
    --input-dir results/alignment \
    --output-file results/consensus_sites.tsv

echo -e "\033[32mMaking logos for plasmid sites\033[0m"
python make_logos.py --alignment-files results/alignment/*.clu

## Sites found in plasmids + combinatorial sites
echo -e "\033[32mMaking alignments for combinatorial sites\033[0m"
python make_alignments.py \
    --att-sites results/att_sites_combinatorial.json \
    --output-dir results/alignment_combinatorial \

echo -e "\033[32mMaking consensus sites for combinatorial sites\033[0m"
python make_consensus_sites.py \
    --input-dir results/alignment_combinatorial \
    --output-file results/consensus_sites_combinatorial.tsv

echo -e "\033[32mMaking logos for combinatorial sites\033[0m"
python make_logos.py --alignment-files results/alignment_combinatorial/*.clu

## Make consensus alignment (not used)

python make_consensus_alignment.py