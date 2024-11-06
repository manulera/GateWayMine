"""
Gather Gateway vectors from the SnapGene installation folder
"""

import argparse
import os
import shutil
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

    # Copy all *dna files into the /data folder
    for folder in [donor_plasmid_folder, destination_plasmid_folder]:
        print(os.path.join(folder, "*.dna"))
        for file in glob.glob(os.path.join(folder, "*.dna")):
            print(file)
            shutil.copy(
                file, os.path.join("data", "snapgene_plasmids", os.path.basename(file))
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--snapgene-dir",
        required=True,
        help="Path to the folder to process (in mac: /Applications/SnapGene.app/Contents/Resources/Plasmids)",
    )
    parser.add_argument(
        "--output-dir",
        default="data/snapgene_plasmids",
        help="Path to the output directory",
    )
    args = parser.parse_args()

    main(args.snapgene_dir, args.output_dir)
