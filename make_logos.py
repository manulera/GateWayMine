"""
Make sequence logos from alignment files.
"""

import glob
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import logomaker
import pandas as pd
import matplotlib.pyplot as plt


def make_logo(alignment: MultipleSeqAlignment, output_file: str):
    freqs = alignment.alignment.frequencies
    matrix = pd.DataFrame.from_dict(freqs, orient="index").T
    matrix = matrix[["A", "C", "G", "T"]]
    fig, ax = plt.subplots(1, 1, figsize=[matrix.shape[0] / 7, 2])
    # Remove axis and labels
    ax.axis("off")
    logomaker.Logo(matrix, ax=ax)
    plt.savefig(output_file)


def main(alignment_files: list[str]):
    for alignment_file in alignment_files:
        alignment = AlignIO.read(alignment_file, "clustal")
        output_file = alignment_file[:-4] + ".svg"
        make_logo(alignment, output_file)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--alignment-files",
        nargs="+",
        help="List of alignment files to process",
        default=glob.glob("results/alignment/*.clu"),
    )
    args = parser.parse_args()
    main(args.alignment_files)
