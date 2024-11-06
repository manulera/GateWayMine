"""
Reads all att sites from the att_sites.json file and creates fasta files for alignments.
"""

import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import glob


def prepare_inputs(input_file: str, output_dir: str):
    with open(input_file, "r") as file:
        att_sites = json.load(file)

    all_records = []
    for att_site, sequences in att_sites.items():
        for i, sequence in enumerate(sequences):
            all_records.append(
                SeqRecord(
                    Seq(sequence),
                    id=f"{att_site}_{i+1}",
                    name=f"{att_site}_{i+1}",
                    description="",
                )
            )

    for i in range(1, 6):
        record_subset = [
            record for record in all_records if record.id.split("_")[0][-1] == str(i)
        ]
        SeqIO.write(record_subset, f"{output_dir}/attX{i}.fasta", "fasta")

    for site_type in "BPLR":
        record_subset = [
            record for record in all_records if record.id.split("_")[0][-2] == site_type
        ]
        SeqIO.write(record_subset, f"{output_dir}/att{site_type}x.fasta", "fasta")

    for site_type in "BPLR":
        for i in range(1, 6):
            record_subset = [
                record
                for record in all_records
                if record.id.split("_")[0][-2] == site_type
                and record.id.split("_")[0][-1] == str(i)
            ]
            SeqIO.write(record_subset, f"{output_dir}/att{site_type}{i}.fasta", "fasta")

    SeqIO.write(all_records, f"{output_dir}/all.fasta", "fasta")


def run_clustalo(output_dir: str, clustalo_bin: str):
    for input_file in glob.glob(f"{output_dir}/*.fasta"):
        output_file = input_file.replace(".fasta", ".clu")
        command = (
            f"{clustalo_bin} -i {input_file} --force --outfmt=clu -o {output_file}"
        )
        os.system(command)


def main(
    input_file: str,
    output_dir: str,
    clustalo_bin: str,
):
    os.makedirs(output_dir, exist_ok=True)
    prepare_inputs(input_file, output_dir)
    run_clustalo(output_dir, clustalo_bin)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--att-sites", type=str, default="results/att_sites.json")
    parser.add_argument("--output-dir", type=str, default="results/alignment")
    parser.add_argument("--clustalo-bin", type=str, default="./clustalo")
    args = parser.parse_args()
    main(
        args.att_sites,
        args.output_dir,
        args.clustalo_bin,
    )
