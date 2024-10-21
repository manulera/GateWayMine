import argparse
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob


def main(folder_path):
    donor_plasmid_folder = os.path.join(folder_path, "Gateway Donor Vectors")
    destination_plasmid_folder = os.path.join(
        folder_path, "Gateway Destination Vectors"
    )

    # Verify that folders exist
    for folder in [folder_path, donor_plasmid_folder, destination_plasmid_folder]:
        if not os.path.isdir(folder):
            raise ValueError(f"The specified folder '{folder}' does not exist.")

    plasmid2sites: dict[str, list[str, str]] = dict()
    site_seq2site_name: dict[str, str] = dict()
    for folder in [donor_plasmid_folder, destination_plasmid_folder]:
        sites: dict[str, list[str]] = dict()
        for plasmid_file in glob.glob(os.path.join(folder, "*.dna")):
            plasmid_record = SeqIO.read(plasmid_file, "snapgene")
            plasmid_name = plasmid_file.split("/")[-1].split(".")[0]
            plasmid2sites[plasmid_name] = list()
            for feature in plasmid_record.features:
                feature: SeqFeature
                label = feature.qualifiers["label"][0]
                if label.startswith("att"):
                    if label not in sites:
                        sites[label] = list()
                    site_seq = str(feature.location.extract(plasmid_record).seq)
                    plasmid2sites[plasmid_name].append(site_seq)
                    sites[label].append(site_seq)

        # Get first key from the dictionary (should start by attP or attR)
        label = next(iter(sites))
        # Use it to name the output file
        output_file = f"{label[:4]}.fa"
        with open(output_file, "w") as handle:
            # Transform the list into sorted list, where the order is given by how many
            # times the site appears in the list
            sites2 = {
                label: sorted(set(site), key=lambda x: -sites[label].count(x))
                for label, site in sorted(sites.items())
            }
            sites = sites2
            for site in sites:
                multi_site = len(sites[site]) > 1
                for i, seq in enumerate(sites[site]):
                    site_name = f"{site}_{i+1}" if multi_site else site
                    site_seq2site_name[seq] = site_name
                    SeqIO.write(
                        SeqRecord(Seq(seq), id=site_name, description=""),
                        handle,
                        "fasta",
                    )
    with open("plasmid2sites.tsv", "w") as handle:
        for plasmid, sites in plasmid2sites.items():
            handle.write(
                f"{plasmid}\t{','.join(site_seq2site_name[site] for site in sites)}\n"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a folder")
    parser.add_argument("--folder", required=True, help="Path to the folder to process")
    args = parser.parse_args()

    main(args.folder)
