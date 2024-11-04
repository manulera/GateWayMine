import argparse
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import json
import re
import glob


def main(input_folder, output_folder):
    # Get all *.dna *.gb *.gbk in /data and subfolders
    files = list()
    for extension in ["dna", "gb", "gbk"]:
        files.extend(glob.glob(os.path.join(input_folder, "**", f"*.{extension}")))

    out_dict = dict()

    for file in files:
        plasmid_dict = dict()
        plasmid_folder = os.path.split(os.path.split(file)[0])[1]
        plasmid_name = os.path.basename(file).split(".")[0]

        file_format = "snapgene" if file.split(".")[-1] == "dna" else "genbank"
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

        out_dict[plasmid_folder + "/" + plasmid_name] = plasmid_dict

    with open(os.path.join(output_folder, "plasmid_site_dict.json"), "w") as handle:
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
        "--output-folder",
        help="Path to the folder to store the output files",
        default="results",
    )
    args = parser.parse_args()

    main(args.input_folder, args.output_folder)
