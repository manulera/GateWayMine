import os
import re
import glob
from Bio import SeqIO

# Read all sequences:
sequences = list()
for file in sorted(glob.glob("att*.fa")):
    sequences += list(SeqIO.parse(file, "fasta"))

groups = {
    "attP": "attP",
    "attB": "attB",
    "attR": "attR",
    "all": ".*",
    "group1": "att.1",
    "group2": "att.2",
    "group3": "att.3",
    "group4": "att.4",
    "no_B": "^(?!.*B).*$",
}

for label, pattern in groups.items():
    seq_subset = [seq for seq in sequences if re.search(pattern, seq.id)]

    with open(f"alignments/input.fa", "w") as handle:
        for seq in seq_subset:
            SeqIO.write(seq, handle, "fasta")

    os.system(
        f"./clustalo --use-kimura -i alignments/input.fa --force --outfmt=clu -o alignments/{label}.clu"
    )

    if label in ["all", "no_B"]:
        # Sort the sequences by group:
        seq_subset = sorted(seq_subset, key=lambda x: x.id[4])
        with open(f"alignments/input.fa", "w") as handle:
            for seq in seq_subset:
                SeqIO.write(seq, handle, "fasta")
        os.system(
            f"./clustalo --use-kimura -i alignments/input.fa --force --outfmt=clu -o alignments/{label}_sorted.clu"
        )


# # Manually amend the attB2 alignment in the 'all' and 'all_sorted' alignments
# for file in [
#     "alignments/all.clu",
#     "alignments/all_sorted.clu",
#     "alignments/group2.clu",
# ]:
#     os.system(
#         f"sed -i '' 's/---ACCCAGCTTTCTTGTACAAAGTGG---/ACCCAGCTTTCTTGTACAAAGTGG------/g' {file}"
#     )
