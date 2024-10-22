"""
Takes a series of plasmids as sequence files and annotates them with the Gateway sequence sites.

Usage:
python annotate_plasmid.py <plasmid_file> <plasmid_file> ...

The output is written to a file with the same name as the input file, but with a .annotated.gb extension.
"""

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
from pydna.utils import shift_location
from Bio.SeqFeature import SeqFeature, SimpleLocation

features = list()

for seq_file in glob.glob("*.fa"):
    features += list(SeqIO.parse(seq_file, "fasta"))

for f in features:
    f.seq = f.seq.upper()


def annotate_plasmid(plasmid: SeqRecord):
    circ_seq: Seq = (plasmid.seq + plasmid.seq).upper()
    for feature in features:
        for rvs in [False, True]:
            query = feature.seq if not rvs else feature.seq.reverse_complement()
            for match in list(circ_seq.search([query])):
                pos, _ = match
                # Skip if the match is in the second half of query sequence (since it's a circular sequence)
                if pos >= len(plasmid.seq):
                    continue
                strand = 1 if not rvs else -1
                location = shift_location(
                    SimpleLocation(pos, pos + len(feature.seq), strand),
                    0,
                    len(plasmid.seq),
                )
                plasmid.features.append(
                    SeqFeature(
                        location,
                        type="protein_bind",
                        qualifiers={"label": feature.id},
                    )
                )
    return plasmid


def main(plasmid_files):
    for plasmid_file in plasmid_files:
        plasmid = SeqIO.read(plasmid_file, "snapgene")
        annotated_plasmid = annotate_plasmid(plasmid)
        SeqIO.write(annotated_plasmid, plasmid_file + ".annotated.gb", "genbank")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python annotate_plasmid.py <plasmid_file> <plasmid_file> ...")
        sys.exit(1)

    plasmid_files = sys.argv[1:]
    main(plasmid_files)
