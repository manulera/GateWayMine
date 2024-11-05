"""
Reads all alignment files in the input directory and creates a consensus sequence for each alignment file.
"""

import os
import glob
import re
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data.IUPACData import ambiguous_dna_values

# We create a dictinary that maps a set of bases to the ambibuous base
# that represents them.
ambiguous_base_dict = {}
for ambiguous, bases in ambiguous_dna_values.items():
    ambiguous_base_dict["".join(sorted(bases))] = ambiguous


def get_consensus_base(all_letters: str) -> str:
    if "-" in all_letters:
        return "-"

    key = "".join(sorted(set(*[all_letters])))
    return ambiguous_base_dict[key]


def get_consensus_sequence(alignment: MultipleSeqAlignment):
    consensus = ""
    for column in range(alignment.alignment.shape[1]):
        all_letters = alignment[:, column]
        consensus += get_consensus_base(all_letters)

    # Remove flanking - characters
    consensus = re.sub(r"^-+", "", consensus)
    consensus = re.sub(r"-+$", "", consensus)

    # Remove flanking N characters
    consensus = re.sub(r"^N+", "", consensus)
    consensus = re.sub(r"N+$", "", consensus)

    # There should be no - left in the consensus
    if "-" in consensus:
        raise ValueError(f"Consensus sequence contains '-': {consensus}")

    return consensus


def main(input_dir: str, output_file: str):

    out_dict = {}
    for alignment_file in glob.glob(f"{input_dir}/*.clu"):
        site_name = os.path.basename(alignment_file).split(".")[0]
        alignment = AlignIO.read(alignment_file, "clustal")
        try:
            consensus = get_consensus_sequence(alignment)
        except ValueError as e:
            raise ValueError(f"Error getting consensus for {site_name}: {e}")

        out_dict[site_name] = consensus

    # Sort the sites alphabetically
    out_dict = dict(sorted(out_dict.items()))

    with open(output_file, "w") as f:
        for site_name, consensus in out_dict.items():
            f.write(f"{site_name}\t{consensus}\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=str, default="results/alignment")
    parser.add_argument(
        "--output-file", type=str, default="results/consensus_sites.tsv"
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_file)
