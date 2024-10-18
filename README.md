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
python main.py --folder <path/to/plasmids>
```

This will create two files in the current directory: `donor_sites.json` and `destination_sites.json`, which contain dictionaries with the Gateway sequence sites for each plasmid.

## Contributing

This could be extended to support other plasmid sources, but for now it only supports the SnapGene plasmid collection. Feel free to submit a PR.

