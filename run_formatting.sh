#! /bin/bash
set -e

python make_plasmid_site_dict.py
python make_plasmid_summary.py
python make_feature_dict.py
python make_unique_sites.py
python make_combinatorial_att_sites.py

