"""
Make a dictionary of features for each plasmid, including:
- plasmid_name
- att_sites
- features (name extracted from gene or label)
- source (snapgene or addgene)
- file (path to the plasmid file in this repository)
- if it's addgene:
    - addgene_id
    - references, if available as links
    - kit (name and url), if it belongs to a kit

This file is then used in the GatewayMine web app.

E.g.
 {
        "source": "snapgene",
        "plasmid_name": "pDEST15",
        "att_sites": [
            "attR1",
            "attR2"
        ],
        "features": [
            "AmpR",
            "AmpR promoter",
            ...
        ]
    },
{
        "source": "addgene",
        "plasmid_name": "pDONR223_C1orf150_p.G2E",
        "sequence-type": "addgene-full",
        "addgene_id": "81309",
        "references": [],
        "kit": {
            "name": "Broad Target Accelerator Plasmid Collections",
            "url": "https://www.addgene.org/1000000103/"
        },
        "att_sites": [
            "attL1",
            "attL2"
        ],
        "features": [
            "L4440",
            "M13 Forward",
            ...
        ]
    },
"""

import json
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import warnings
from tqdm import tqdm


def main(plasmid_summary_file, plasmid_site_dict_file, output_file):
    with open(plasmid_summary_file) as f:
        plasmid_summary = json.load(f)

    with open(plasmid_site_dict_file) as f:
        plasmid_site_dict = json.load(f)

    for plasmid in tqdm(plasmid_summary, desc="Extracting plasmid features"):
        # If no att sites, skip and remove from plasmid summary
        if len(plasmid_site_dict[plasmid["file"]]) == 0:
            plasmid_summary.remove(plasmid)
            continue

        if plasmid["source"] == "snapgene":
            with open(plasmid["file"], "br") as f:
                record = SeqIO.read(f, "snapgene")
        elif plasmid["source"] == "addgene":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
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
