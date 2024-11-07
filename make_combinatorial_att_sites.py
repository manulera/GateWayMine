"""
Recombine all att sites found in plasmids via BP or LR reaction to to generate
even more att sites.
"""

import json
import re
from Bio.Data.IUPACData import ambiguous_dna_values as _ambiguous_dna_values
import copy

overlap_dict = {
    "1": "twtGTACAAAaaa",
    "2": "twtGTACAAGaaa",
    "3": "twtGTATAATaaa",
    "4": "twtGTATAGAaaa",
    "5": "twtGTATACAaaa",
}

ambiguous_only_dna_values = {**_ambiguous_dna_values}
for normal_base in "ACGT":
    del ambiguous_only_dna_values[normal_base]


def compute_regex_site(site: str) -> str:
    """
    Compute the regex for a given site containing ambiguous bases.
    """
    upper_site = site.upper()
    for k, v in ambiguous_only_dna_values.items():
        if len(v) > 1:
            upper_site = upper_site.replace(k, f"[{''.join(v)}]")

    # Make case insensitive
    upper_site = f"(?i){upper_site}"
    return upper_site


overlap_regex = {k: compute_regex_site(v) for k, v in overlap_dict.items()}


def main(att_sites_file, output_file, output_file_combinatorial_only):
    with open(att_sites_file) as f:
        att_sites = json.load(f)

    # All sequences should containt the corresponding overlap
    # regex only once
    for site in att_sites:
        site_number = site[-1]
        pattern = overlap_regex[site_number]
        for seq in att_sites[site]:
            if len(re.findall(pattern, seq)) < 1:
                raise ValueError(
                    f"Site {site} with {seq} does not contain overlap regex"
                )
            elif len(re.findall(pattern, seq)) > 1:
                raise ValueError(
                    f"Site {site} with {seq} contains overlap regex more than once"
                )
    # All possible reactions (inputs, output)
    reactions = (
        ("PB", "L"),
        ("BP", "R"),
        ("RL", "B"),
        ("LR", "P"),
    )

    computed_att_sites = copy.deepcopy(att_sites)

    for input_types, output_type in reactions:
        left_input, right_input = input_types

        for i in range(1, 6):
            left_site = f"att{left_input}{i}"
            right_site = f"att{right_input}{i}"
            output_site = f"att{output_type}{i}"
            print(f"{left_site} + {right_site} -> {output_site}")

            left_seqs = att_sites[left_site]
            right_seqs = att_sites[right_site]
            pattern = overlap_regex[str(i)]
            for left_seq in left_seqs:
                for right_seq in right_seqs:
                    end_left = re.search(pattern, left_seq).start()
                    start_right = re.search(pattern, right_seq).start()
                    computed_seq = left_seq[:end_left] + right_seq[start_right:]
                    computed_att_sites[output_site].append(computed_seq)

    for site in computed_att_sites:
        computed_att_sites[site] = list(sorted(set(computed_att_sites[site])))

    with open(output_file, "w") as f:
        json.dump(computed_att_sites, f, indent=4)

    # Keep only the new ones
    for computed_site in computed_att_sites:
        computed_att_sites[computed_site] = list(
            sorted(
                set(computed_att_sites[computed_site]) - set(att_sites[computed_site])
            )
        )

    with open(output_file_combinatorial_only, "w") as f:
        json.dump(computed_att_sites, f, indent=4)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--att-sites-file", type=str, default="results/att_sites.json")
    parser.add_argument(
        "--output-file", type=str, default="results/att_sites_combinatorial.json"
    )
    parser.add_argument(
        "--output-file-combinatorial-only",
        type=str,
        default="results/att_sites_combinatorial_only.json",
    )
    args = parser.parse_args()
    main(
        args.att_sites_file,
        args.output_file,
        args.output_file_combinatorial_only,
    )
