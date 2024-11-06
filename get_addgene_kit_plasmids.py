"""
Download plasmids from AddGene kits, using the addgene_kit_plasmids.json dictionary.
"""

import httpx
import os
import asyncio
from bs4 import BeautifulSoup
import argparse
import json
import re


async def get_addgene_plasmid(addgene_id: str, output_dir: str):
    url = f"https://www.addgene.org/{addgene_id}/sequences/"
    async with httpx.AsyncClient() as client:
        resp = await client.get(url)
    if resp.status_code == 404:
        raise ValueError(f"The requested plasmid does not exist on Addgene: {url}")
    soup = BeautifulSoup(resp.content, "html.parser")

    sequence_file_url_dict = dict()
    for _type in [
        "depositor-full",
        "depositor-partial",
        "addgene-full",
        "addgene-partial",
    ]:
        sequence_file_url_dict[_type] = []
        if soup.find(id=_type) is not None:
            sequence_file_url_dict[_type] = [
                a.get("href")
                for a in soup.find(id=_type).findAll(class_="genbank-file-download")
            ]

    for _type in ["addgene-full", "depositor-full"]:
        if len(sequence_file_url_dict[_type]) > 0:
            for seq_url in sequence_file_url_dict[_type]:
                output_name = os.path.join(output_dir, f"{addgene_id}.{_type}.gb")
                # Get file content from url
                async with httpx.AsyncClient() as client:
                    resp = await client.get(seq_url)
                with open(output_name, "wb") as handle:
                    handle.write(resp.content)
                    return True
    return False


async def main(input_file: str, output_dir: str):

    with open(input_file, "r") as no_sequence_file:
        addgene_plasmids = json.load(no_sequence_file)

    addgene_ids_to_download = []
    for kit_name, kit in addgene_plasmids.items():
        addgene_ids_to_download.extend(kit["addgene_ids"])

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
    parser = argparse.ArgumentParser(
        description=__doc__,
    )
    parser.add_argument(
        "--input_file", type=str, default="data/addgene_kit_plasmids.json"
    )
    parser.add_argument("--output-dir", type=str, default="data/addgene_plasmids")
    args = parser.parse_args()
    asyncio.run(main(args.input_file, args.output_dir))
