"""
Reads all alignment files in the input directory and creates a consensus
sequence for each alignment file. In addition, it makes a "merged" consensus
sequences using all the sites of a certain type (e.g. all attB sites), excluding
the overlap sequence, common to all att sites of a given number
(e.g. all attX1 sites contain `twtGTACAAAaaa` as the overlap sequence).
"""

import os
import glob
import re
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data.IUPACData import ambiguous_dna_values
from make_combinatorial_att_sites import overlap_dict, compute_regex_site

# We create a dictinary that maps a set of bases to the ambibuous base
# that represents them.
ambiguous_base_dict = {}
for ambiguous, bases in ambiguous_dna_values.items():
    ambiguous_base_dict["".join(sorted(bases))] = ambiguous


def validate_consensus_sequence(
    consensus: str, alignment: MultipleSeqAlignment, site_name: str
):
    pattern = compute_regex_site(consensus)
    for seq in alignment:
        if seq.id.startswith(site_name) and re.search(pattern, str(seq.seq)) is None:
            raise ValueError(
                f"{site_name}: Consensus sequence {consensus} does not match sequence from {seq.id}: {seq.seq}"
            )
        elif site_name.endswith("x") and re.search(pattern, str(seq.seq)) is None:
            raise ValueError(
                f"{site_name}: Consensus sequence {consensus} does not match sequence from {seq.id}: {seq.seq}"
            )


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

        validate_consensus_sequence(consensus, alignment, site_name)
        out_dict[site_name] = consensus

        # These are alignments that contain all sites of a certain type
        # (e.g. attB1, attB2, etc.)
        # We make a consensus where the only differing sequence
        # is the overlap sequence
        if site_name.endswith("x"):
            # This is the consensus of the overlap sequence present in all sites
            # (capital letters in overlap_dict)
            # twtGTACAAAaaa
            # twtGTACAAGaaa
            # twtGTATAATaaa
            # twtGTATAGAaaa
            # twtGTATACAaaa
            #    GTAYAVD
            core_consensus = "GTAYAVD"
            # Should be present only once in the consensus
            if len(re.findall(core_consensus, consensus)) != 1:
                raise ValueError(f"Expected 1 core sequence for {site_name}")
            for num in range(1, 6):
                core = overlap_dict[str(num)][3:-3]
                merged_consensus = consensus.replace(core_consensus, core)
                merged_site = f"merged_{site_name[:-1]}{num}"

                validate_consensus_sequence(
                    merged_consensus, alignment, f"{site_name[:-1]}{num}"
                )
                out_dict[merged_site] = merged_consensus

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
