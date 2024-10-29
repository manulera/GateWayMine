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
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import reverse_complement
import re

for normal_base in "ACGT":
    del ambiguous_dna_values[normal_base]


def compute_regex_site(site: str):
    for k, v in ambiguous_dna_values.items():
        if len(v) > 1:
            site = site.replace(k, f"[{''.join(v)}]")
    return site


gateway_features = dict()
with open("gateway_features.txt") as f:
    for line in f:
        line = line.strip().split("\t")
        gateway_features[line[0]] = line[1]

features = list()

for seq_file in glob.glob("*.fa"):
    features += list(SeqIO.parse(seq_file, "fasta"))

for f in features:
    f.seq = f.seq.upper()


def annotate_plasmid(plasmid: SeqRecord):
    circ_seq: Seq = (plasmid.seq + plasmid.seq).upper()
    for feature_name, feature_seq in gateway_features.items():
        for rvs in [False, True]:
            query = feature_seq if not rvs else reverse_complement(feature_seq)
            for match in re.finditer(compute_regex_site(query), str(circ_seq)):
                print(feature_name)
                pos = match.start()
                match_len = match.end() - match.start()

                # Skip if the match is in the second half of query sequence (since it's a circular sequence)
                if pos >= len(plasmid.seq):
                    continue
                strand = 1 if not rvs else -1
                location = shift_location(
                    SimpleLocation(pos, pos + match_len, strand),
                    0,
                    len(plasmid.seq),
                )
                plasmid.features.append(
                    SeqFeature(
                        location,
                        type="protein_bind",
                        qualifiers={"label": feature_name},
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
