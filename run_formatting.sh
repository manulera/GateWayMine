#! /bin/bash
set -e

echo -e "\033[32mRunning make_plasmid_site_dict.py\033[0m"
python make_plasmid_site_dict.py

echo -e "\033[32mRunning make_plasmid_summary.py\033[0m"
python make_plasmid_summary.py

echo -e "\033[32mRunning make_feature_dict.py\033[0m"
python make_feature_dict.py

echo -e "\033[32mRunning make_unique_sites.py\033[0m"
python make_unique_sites.py

echo -e "\033[32mRunning make_combinatorial_att_sites.py\033[0m"
python make_combinatorial_att_sites.py

