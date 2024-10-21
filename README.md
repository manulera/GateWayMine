# Gateway Sequence Sites

This repository contains the code and data for extracting Gateway sequence sites from donor and destination plasmids from the SnapGene plasmid collection.

These files can be found inside the SnapGene installation folder, in Mac the path is `/Applications/SnapGene.app/Contents/Resources/Plasmids`. To access it go to `Applications`, then right click SnapGene and select "Show Package Contents".

## Installing dependencies

### With pip

```bash
python -m venv .venv
source .venv/bin/activate
pip install biopython
```

### With poetry

```bash
poetry install
poetry shell # Activate the virtual environment
```

## Running the script

```bash
# without alignment
python main.py --folder <path/to/plasmids>
# with alignment
python main.py --folder <path/to/plasmids> --align
```

This will create two files in the current directory: `donor_sites.json` and `destination_sites.json`, which contain dictionaries with the Gateway sequence sites for each plasmid.

If you want to align the sites, you can use the `--align` flag. This will run `clustalo` for each site and create an aligment file in the `alignments` folder. You need to have the `clustalo` binary in the project directory. You can download it from [here](http://www.clustal.org/omega/), and rename it to `clustalo`.

## Contributing

This could be extended to support other plasmid sources, but for now it only supports the SnapGene plasmid collection. Feel free to submit a PR.

## Conclusion

Marked with `x` are the positions that constitute the `O` overlap of the attP and attRsites

```
                                   xxxxxxx
attP1_1      GCTTTTTTATAATGCCAACTTTGTACAAAAAAGTTGAACGAGAAACGTAAAATGATATAA
attP1_2      GCTTTTTTATAATGCCAACTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAA
attP1_3      GCTTTTTTATAATGCCAAGTTTGTACAAAAAAGCAGAACGAGAAACGTAAAATGATATAA
attP1_4      GCTTTTTTATAATGCCAACTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAA
attP2_1      GCTTTCTTATAATGCCAACTTTGTACAAGAAAGCTGAACGAGAAACGTAAAATGATATAA
attP2_2      GCTTTCTTATAATGCCCACTTTGTACAAGAAAGCTGAACGAGAAACGTAAAATGATATAA
attP2_3      GCTTTTTTATAATGCCAACTTTGTACAAGAAAGTTGAACGAGAAACGTAAAATGATATAA
attP3_1      GCTTTTTTATAATGCCAACTTTGTATAATAAAGTTGAACGAGAAACGTAAAATGATATAA
attP3_2      GCTTTTTTATAATGCCAACTTTGTATAATAAAGTTGAACGAGAAACGTAAAATGATATAA
attP4_1      GCTTTCTTATAATGCCAACTTTGTATAGAAAAGTTGAACGAGAAACGTAAAATGATATAA
attP4_2      GCTTTTTTATAATGCCAACTTTGTATAGAAAAGTTGAACGAGAAACGTAAAATGATATAA
attP5_1      GCTTTTTTATAATGCCAACTTTGTATACAAAAGTTGAACGAGAAACGTAAAATGATATAA
attR1_1      --------------ACAAGTTTGTACAAAAAAGCTGAACGAGAAACGTAAAATGATATAA
attR2_1      --------------ACCACTTTGTACAAGAAAGCTGAACGAGAAACGTAAAATGATATAA
attR3_1      ---------------CAACTATGTATAATAAAGTTGAACGAGAAACGTAAAATGATATAA
attR4_1      ---------------CAACTTTGTATAGAAAAGTTGAACGAGAAACGTAAAATGATATAA
attR4_2      --------------TCAACTTTGTATAGAAAAGTTGAACGAGAAACGTAAAATGATATAA
                            * * * **** *  ****  *************************
```