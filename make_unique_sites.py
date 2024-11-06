"""
From a dictionary plasmid -> sites, extract all unique sites of each
type, and create a dictionary site_type -> sequences.

E.g. 

{
    "attR1": [
        "ACAACTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATGCAGTCACTATG",
        "ACAACTTTGTACAAAAAAGTTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATGCAGTCACTATG",
        "ACAAGTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATCCAGTCA",
        ...
        ]
}
"""

import json


def main(input_file: str, output_file: str):
    with open(input_file, "r") as f:
        plasmid_site_dict = json.load(f)

    unique_site_dict = dict()
    for plasmid, sites in plasmid_site_dict.items():
        for site_type, site_seqs in sites.items():
            # Remove those that contain non-ACGT characters
            site_seqs = [
                seq for seq in site_seqs if all(base in "ACGT" for base in seq)
            ]
            if site_type not in unique_site_dict:
                unique_site_dict[site_type] = set()
            unique_site_dict[site_type].update(site_seqs)

    # Sort sequences alphabetically
    for site_type in unique_site_dict:
        unique_site_dict[site_type] = sorted(unique_site_dict[site_type])

    with open(output_file, "w") as f:
        json.dump(unique_site_dict, f, indent=4)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_file",
        type=str,
        help="Path to plasmid site dictionary",
        default="results/plasmid_site_dict.json",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        help="Path to output file",
        default="results/att_sites.json",
    )
    args = parser.parse_args()
    main(args.input_file, args.output_file)
