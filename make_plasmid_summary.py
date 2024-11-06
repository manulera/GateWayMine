"""
Make a dictionary of features for each plasmid
"""

import glob
import os
import json


def main(
    input_folder, output_file, all_gateway_plasmids, addgene_kits, addgene_articles
):
    addgene_id2name = dict()
    with open(all_gateway_plasmids) as f:
        for line in f:
            addgene_id, name = line.strip().split("\t")[:2]
            addgene_id2name[addgene_id] = name

    addgene_id2references = dict()
    with open(addgene_articles) as f:
        for line in f:
            ls = line.strip().split("\t")
            addgene_id = ls[0]
            references = ls[1:]
            addgene_id2references[addgene_id] = references

    addgene_id2kit = dict()
    with open(addgene_kits) as f:
        kit_dict = json.load(f)
        for kit_name, kit_info in kit_dict.items():
            # Copy of kit info without addgene_ids and names
            short_kit_info = {
                k: v
                for k, v in kit_info.items()
                if k not in ["addgene_ids", "plasmid_names"]
            }
            for i, addgene_id in enumerate(kit_info["addgene_ids"]):
                addgene_id2kit[addgene_id] = short_kit_info
                addgene_id2name[addgene_id] = kit_info["plasmid_names"][i]

    files = list()
    for extension in ["dna", "gb", "gbk"]:
        files.extend(glob.glob(os.path.join(input_folder, "**", f"*.{extension}")))

    output = list()
    for file in files:
        basename = os.path.basename(file)
        plasmid_dict = dict()
        if basename.endswith(".dna"):
            plasmid_dict["source"] = "snapgene"
            plasmid_dict["plasmid_name"] = os.path.basename(file).split(".")[0]
        elif "addgene-full" in basename or "depositor-full" in basename:
            addgene_id = basename.split(".")[0]
            plasmid_dict["source"] = "addgene"
            plasmid_dict["plasmid_name"] = addgene_id2name[addgene_id]
            plasmid_dict["sequence-type"] = basename.split(".")[1]
            plasmid_dict["addgene_id"] = addgene_id
            if addgene_id in addgene_id2references:
                plasmid_dict["references"] = addgene_id2references[addgene_id]
            else:
                plasmid_dict["references"] = []
            if addgene_id in addgene_id2kit:
                plasmid_dict["kit"] = addgene_id2kit[addgene_id]
            else:
                plasmid_dict["kit"] = None
        else:
            raise ValueError(f"Unknown file type: {basename}")
        plasmid_dict["file"] = file
        output.append(plasmid_dict)

    with open(output_file, "w") as f:
        json.dump(output, f, indent=2)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__,
    )
    parser.add_argument(
        "--input-folder",
        help="Path to the folder containing the sequence files",
        default="data",
    )
    parser.add_argument(
        "--all-gateway-plasmids",
        help="Path to the file containing all gateway plasmids",
        default="data/all_gateway_plasmids.tsv",
    )
    parser.add_argument(
        "--addgene-kits",
        help="Path to the json file containing addgene kits",
        default="data/addgene_kit_plasmids.json",
    )
    parser.add_argument(
        "--addgene-articles",
        help="Path to the tsv file containing addgene articles",
        default="data/addgene_article_refs.tsv",
    )
    parser.add_argument(
        "--output-file",
        help="Path to the output file",
        default="results/plasmid_summary.json",
    )
    args = parser.parse_args()

    main(
        args.input_folder,
        args.output_file,
        args.all_gateway_plasmids,
        args.addgene_kits,
        args.addgene_articles,
    )
