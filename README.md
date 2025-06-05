# PDB EXTRACTION

This repository contains tools for parsing PDB files and analyzing intermolecular contacts.

## Scripts

### `pdb_contact_analysis.py`

A utility that reads all PDB files in `Inputs/PDB-Inputs/` (or a user specified directory) and calculates residues contacting a user defined selection.

**Features**

- Supports selections of residues, ranges or entire chains (e.g. `"A:100"`, `"A:1-10,B:"`, `"LIG"`).
- Calculates minimum atom distance to the selection using Bio.Python.
- Annotates each residue with domain (`N-term`, `LBD`, `DNA`) and whether the interaction is `cis` or `trans` relative to the selected chain(s).
- Processes multiple PDBs (states) and reports shared vs. state specific contacts.

**Basic usage**

```bash
python pdb_contact_analysis.py --selection "A:100" --cutoff 4.0
```

For interactive use in Jupyter, open `pdb_contact_analysis.py` and adjust the
variables `PDB_DIR`, `SELECTION`, and `CUTOFF` near the top of the file before
running `run_analysis()`.

Use `--domain` or `--interaction` to filter results. The default PDB directory is
`C:\Users\shirk\Documents\Github\AclR\PDB-Extraction\PDB-Inputs\`.
