#! /bin/bash
set -e

# No need to repeat this
# python get_snapgene_files.py --folder /Applications/SnapGene.app/Contents/Resources/Plasmids

python get_all_gateway_plasmids.py
python get_addgene_kits_info.py
python get_addgene_article_refs.py



python get_other_plasmids.py
python get_addgene_kit_plasmids.py

