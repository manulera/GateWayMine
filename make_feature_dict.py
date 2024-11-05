"""
Make a dictionary of features for each plasmid
"""

import json
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


def main(plasmid_summary_file, plasmid_site_dict_file, output_file):
    with open(plasmid_summary_file) as f:
        plasmid_summary = json.load(f)

    with open(plasmid_site_dict_file) as f:
        plasmid_site_dict = json.load(f)

    for plasmid in plasmid_summary:
        # If no att sites, skip and remove from plasmid summary
        if len(plasmid_site_dict[plasmid["file"]]) == 0:
            plasmid_summary.remove(plasmid)
            continue

        if plasmid["source"] == "snapgene":
            with open(plasmid["file"], "br") as f:
                record = SeqIO.read(f, "snapgene")
        elif plasmid["source"] == "addgene":
            with open(plasmid["file"], "r") as f:
                record = SeqIO.read(f, "genbank")
        features: list[SeqFeature] = record.features
        plasmid["att_sites"] = list(sorted(plasmid_site_dict[plasmid["file"]].keys()))
        plasmid["features"] = set()
        for feature in features:
            if "label" in feature.qualifiers:
                plasmid["features"].add(feature.qualifiers["label"][0])
            elif "gene" in feature.qualifiers:
                plasmid["features"].add(feature.qualifiers["gene"][0])

        plasmid["features"] = list(sorted(plasmid["features"]))

        # Remove file from plasmid summary
        plasmid.pop("file")

    with open(output_file, "w") as f:
        json.dump(plasmid_summary, f, indent=4)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__,
    )
    parser.add_argument(
        "--plasmid-summary",
        help="Path to the plasmid summary file",
        default="results/plasmid_summary.json",
    )
    parser.add_argument(
        "--plasmid-site-dict",
        help="Path to the plasmid site dictionary",
        default="results/plasmid_site_dict.json",
    )
    parser.add_argument(
        "--output-file",
        help="Path to the output file",
        default="results/plasmid_features.json",
    )
    args = parser.parse_args()

    main(args.plasmid_summary, args.plasmid_site_dict, args.output_file)
