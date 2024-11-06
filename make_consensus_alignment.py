"""
Align the consensus sequences (unused).
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from make_alignments import run_clustalo


def prepare_inputs(file_name: str, output_dir: str, combinatorial: bool):
    seq_dict = {}
    with open(file_name, "r") as file:
        for line in file:
            name, seq = line.strip().split("\t")
            seq_dict[name] = seq

    all_records = []
    for name, seq in seq_dict.items():
        if combinatorial:
            name = f"{name}_comb"
        all_records.append(SeqRecord(Seq(seq), id=name, name=name, description=""))

    if combinatorial:
        output_file = "consensus_sites_combinatorial.fasta"
    else:
        output_file = "consensus_sites.fasta"

    SeqIO.write(all_records, f"{output_dir}/{output_file}", "fasta")


def main(
    consensus_sites: str,
    consensus_sites_combinatorial: str,
    output_dir: str,
    clustalo_bin: str,
):
    prepare_inputs(consensus_sites, output_dir, combinatorial=False)
    prepare_inputs(consensus_sites_combinatorial, output_dir, combinatorial=True)
    run_clustalo(output_dir, clustalo_bin)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--consensus-sites", type=str, default="results/consensus_sites.tsv"
    )
    parser.add_argument(
        "--consensus-sites-combinatorial",
        type=str,
        default="results/consensus_sites_combinatorial.tsv",
    )
    parser.add_argument("--output-dir", type=str, default="results/aligment_consensus")
    parser.add_argument("--clustalo-bin", type=str, default="./clustalo")
    args = parser.parse_args()

    main(
        args.consensus_sites,
        args.consensus_sites_combinatorial,
        args.output_dir,
        args.clustalo_bin,
    )
