"""
Make a dictionary of att sites for each plasmid, where the keys are the paths of plasmids, and the values
are the att sites found in each plasmid, e.g:

    "data/snapgene_plasmids/pDEST15.dna": {
        "attR1": [
            "ACAAGTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATG"
        ],
        "attR2": [
            "ACCACTTTGTACAAGAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCACTATG"
        ]
    },
"""

import argparse
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import json
import re
import glob
import warnings
from tqdm import tqdm


def main(input_folder, output_file):
    # Get all *.dna *.gb *.gbk in /data and subfolders
    files = list()
    for extension in ["dna", "gb", "gbk"]:
        files.extend(glob.glob(os.path.join(input_folder, "**", f"*.{extension}")))

    out_dict = dict()

    for file in tqdm(files, desc="Processing plasmid files"):
        plasmid_dict = dict()

        file_format = "snapgene" if file.split(".")[-1] == "dna" else "genbank"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plasmid_record = SeqIO.read(file, file_format)

        # Find att sites
        for feature in plasmid_record.features:
            feature: SeqFeature
            if "label" not in feature.qualifiers:
                continue
            label = feature.qualifiers["label"][0]
            if re.match(r"att[PLBR]\d+", label):
                if label not in plasmid_dict:
                    plasmid_dict[label] = list()
                site_seq = str(feature.location.extract(plasmid_record).seq)
                plasmid_dict[label].append(site_seq)

        out_dict[file] = plasmid_dict

    with open(output_file, "w") as handle:
        json.dump(out_dict, handle, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract att sites from sequence files"
    )
    parser.add_argument(
        "--input-folder",
        help="Path to the folder containing the sequence files",
        default="data",
    )
    parser.add_argument(
        "--output-file",
        help="Path to the output file",
        default="results/plasmid_site_dict.json",
    )
    args = parser.parse_args()

    main(args.input_folder, args.output_file)
