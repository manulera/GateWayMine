"""
Download plasmids from AddGene listed in a table file, into a given output directory.

Only new plasmids (not present in the output folder) are downloaded. The plasmids are
identified by their AddGene ids, which are the first column in the table file. If a plasmid
without a sequence is encountered, its id is added to a separate file "no_sequence.txt"
in the output directory, to prevent attempts to download it in subsequent runs.
"""

from get_addgene_kit_plasmids import get_addgene_plasmid
import argparse
import asyncio
import re
import os


async def main(input_file: str, output_dir: str):

    addgene_ids_to_download = []
    with open(input_file, "r") as no_sequence_file:
        for line in no_sequence_file:
            addgene_ids_to_download.append(line.strip().split("\t")[0])

    # Validate the ids
    for addgene_id in addgene_ids_to_download:
        if not re.match(r"\d+", addgene_id):
            raise ValueError(f"Invalid addgene id: {addgene_id}")

    existing_files = os.listdir(output_dir)
    existing_addgene_ids = []
    for no_sequence_file in existing_files:
        addgene_id = re.search(r"(\d+)\..+\.gb", no_sequence_file)
        if addgene_id:
            existing_addgene_ids.append(addgene_id.group(1))

    no_sequence_plasmids = []
    with open(os.path.join(output_dir, "no_sequence.txt"), "r") as no_sequence_file:
        no_sequence_plasmids = no_sequence_file.read().splitlines()

    with open(os.path.join(output_dir, "no_sequence.txt"), "a") as no_sequence_file:
        for addgene_id in addgene_ids_to_download:
            if addgene_id in existing_addgene_ids or addgene_id in no_sequence_plasmids:
                print("Skipping", addgene_id)
                continue
            print("Downloading", addgene_id)
            if not await get_addgene_plasmid(addgene_id, output_dir):
                no_sequence_file.write(addgene_id + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input_file", type=str, default="data/all_gateway_plasmids.tsv"
    )
    parser.add_argument("--output-dir", type=str, default="data/addgene_plasmids")
    args = parser.parse_args()
    asyncio.run(main(args.input_file, args.output_dir))
