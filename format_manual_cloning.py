import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob


def main(folder_path):
    out_dict = dict()
    for file in glob.glob(os.path.join(folder_path, "*.dna")):
        with open(file, "rb") as handle:
            record = SeqIO.read(handle, "snapgene")
            file_name = file.split("/")[-1].split(".")[0]

            # The name contains the attPx_y, where y is the version
            # of the attPx site
            attP_sites = dict()
            for attPsite in file_name.split("-")[1:]:
                # We make an x: y dictionary
                splitted = attPsite.split("_")
                if len(splitted) == 2:
                    attP_sites[splitted[0][-1]] = splitted[1]
                else:
                    attP_sites[splitted[0][-1]] = ""

            for feature in record.features:
                label = feature.qualifiers["label"][0]
                if feature.type != "primer_bind" and label.startswith("att"):
                    # The x in attLx
                    attL_x = label[-1]
                    if attP_sites[attL_x]:
                        attL_label = label + "_" + attP_sites[attL_x]
                    else:
                        attL_label = label
                    attL_seq = str(feature.location.extract(record).seq)
                    out_dict[attL_label] = attL_seq

    for att_type in ["L", "B"]:
        with open(f"att{att_type}_manual.fa", "w") as handle:
            for label, seq in sorted(out_dict.items()):
                if label.startswith(f"att{att_type}"):
                    SeqIO.write(
                        SeqRecord(Seq(seq), id=label, description=""),
                        handle,
                        "fasta",
                    )


if __name__ == "__main__":
    main("manual_cloning")
