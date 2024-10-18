import argparse
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
import json


def main(folder_path, align):
    donor_plasmid_folder = os.path.join(folder_path, "Gateway Donor Vectors")
    destination_plasmid_folder = os.path.join(
        folder_path, "Gateway Destination Vectors"
    )

    # Verify that folders exist
    for folder in [folder_path, donor_plasmid_folder, destination_plasmid_folder]:
        if not os.path.isdir(folder):
            raise ValueError(f"The specified folder '{folder}' does not exist.")

    for folder in [donor_plasmid_folder, destination_plasmid_folder]:
        sites: dict[str, set[str]] = dict()
        for plasmid_file in glob.glob(os.path.join(folder, "*.dna")):
            plasmid_record = SeqIO.read(plasmid_file, "snapgene")
            for feature in plasmid_record.features:
                feature: SeqFeature
                label = feature.qualifiers["label"][0]
                if label.startswith("att"):
                    if label not in sites:
                        sites[label] = set()
                    sites[label].add(str(feature.location.extract(plasmid_record).seq))

        # Transform the sets into sorted lists
        sites = {label: sorted(list(site)) for label, site in sorted(sites.items())}
        name = "destination" if "Destination" in os.path.basename(folder) else "donor"
        with open(f"{name}_sites.json", "w") as f:
            json.dump(sites, f, indent=4)

        if not align:
            continue
        for label, sites in sites.items():
            if len(sites) > 1:
                print(f"Aligning {label}")
                # Create a temp fasta file with the sites
                with open(f"alignments/input.fa", "w") as handle:
                    for i, site in enumerate(sites):
                        SeqIO.write(
                            SeqRecord(Seq(site), id=f"{label}_{i}"), handle, "fasta"
                        )

                # Run clustalo
                os.system(
                    f"./clustalo -i alignments/input.fa --force --outfmt=clu -o alignments/{label}_aligned.clu"
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a folder")
    parser.add_argument("--folder", required=True, help="Path to the folder to process")
    parser.add_argument("--align", action="store_true", help="Run alignments")
    args = parser.parse_args()

    main(args.folder, args.align)
