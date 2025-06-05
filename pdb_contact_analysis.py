import os
from collections import defaultdict
from typing import List, Tuple, Dict

import pandas as pd
from Bio.PDB import PDBParser, is_aa

from Bio.PDB.NeighborSearch import NeighborSearch

# Default parameters for notebook use. Adjust these before running the analysis
# Double backslash at the end keeps the trailing separator on Windows.
PDB_DIR = r"C:\Users\shirk\Documents\Github\AclR\PDB-Extraction\PDB-Inputs\\"
SELECTION = "A:100"
CUTOFF = 4.0
FILTER_DOMAIN = None  # set to 'N-term', 'LBD', or 'DNA' to filter
FILTER_INTERACTION = None  # set to 'cis' or 'trans' to filter


DOMAIN_MAP = {
    'N-term': range(1, 28),
    'LBD': range(28, 190),
    'DNA': range(190, 253),
}


def assign_domain(resnum: int) -> str:
    for name, rng in DOMAIN_MAP.items():
        if resnum in rng:
            return name
    return 'Unknown'


def parse_selection(model, selection: str) -> Tuple[List, List[str]]:
    """Return atoms for the selection and list of chains involved."""
    atoms = []
    chains_involved = set()

    tokens = [tok.strip() for tok in selection.split(',') if tok.strip()]
    for token in tokens:
        if ':' in token:
            chain_id, res_part = token.split(':', 1)
            if chain_id not in model:
                continue
            chain = model[chain_id]
            chains_involved.add(chain_id)
            if not res_part:
                for res in chain:
                    atoms.extend(res.get_atoms())
            else:
                for seg in res_part.split(';'):
                    seg = seg.strip()
                    if '-' in seg:
                        start, end = seg.split('-', 1)
                        start, end = int(start), int(end)
                        for res in chain:
                            if start <= res.id[1] <= end:
                                atoms.extend(res.get_atoms())
                    elif seg.isdigit():
                        try:
                            res = chain[(' ', int(seg), ' ')]
                            atoms.extend(res.get_atoms())
                        except KeyError:
                            pass
                    else:  # resname or ligand id
                        for res in chain:
                            if res.resname == seg or str(res.id[1]) == seg:
                                atoms.extend(res.get_atoms())
        else:
            # search by resname across all chains
            for chain in model:
                for res in chain:
                    if res.resname == token:
                        atoms.extend(res.get_atoms())
                        chains_involved.add(chain.id)
    return atoms, list(chains_involved)


def compute_contacts(model, target_atoms: List, target_chains: List[str], cutoff: float) -> Dict[Tuple[str, str, int], float]:
    all_atoms = [atom for chain in model for res in chain for atom in res]
    search = NeighborSearch(all_atoms)

    contact_residues = {}
    for chain in model:
        for res in chain:
            if not is_aa(res):
                continue
            if res.id[1] > 252:
                continue
            atoms = list(res.get_atoms())
            min_dist = float('inf')
            for atom in atoms:
                neighbors = search.search(atom.coord, cutoff, level='A')
                for nb in neighbors:
                    if nb in atoms:
                        continue
                    if nb not in target_atoms:
                        continue
                    dist = atom - nb
                    if dist < min_dist:
                        min_dist = dist
            if min_dist != float('inf'):
                contact_residues[(chain.id, res.resname, res.id[1])] = min_dist
    return contact_residues


def analyze_pdb(pdb_file: str, selection: str, cutoff: float) -> pd.DataFrame:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)
    model = structure[0]

    target_atoms, target_chains = parse_selection(model, selection)
    contacts = compute_contacts(model, target_atoms, target_chains, cutoff)

    rows = []
    for (chain_id, resname, resnum), dist in contacts.items():
        domain = assign_domain(resnum)
        cis = chain_id in target_chains
        rows.append({
            'State': os.path.basename(pdb_file),
            'Chain': chain_id,
            'Residue': resname,
            'Resnum': resnum,
            'Domain': domain,
            'Distance': dist,
            'Interaction': 'cis' if cis else 'trans'
        })
    df = pd.DataFrame(rows)
    return df


def compare_states(dfs: List[pd.DataFrame]) -> pd.DataFrame:
    """Return dataframe summarizing shared and state-specific contacts."""
    all_contacts = defaultdict(set)
    for df in dfs:
        if df.empty:
            continue
        state = df['State'].iloc[0]
        for _, row in df.iterrows():
            key = (row['Chain'], row['Residue'], row['Resnum'], row['Interaction'], row['Domain'])
            all_contacts[key].add(state)

    rows = []
    states = [df['State'].iloc[0] for df in dfs if not df.empty]
    for key, present_states in all_contacts.items():
        chain, resname, resnum, interaction, domain = key
        rows.append({
            'Chain': chain,
            'Residue': resname,
            'Resnum': resnum,
            'Domain': domain,
            'Interaction': interaction,
            'States': ','.join(sorted(present_states)),
            'Shared': len(present_states) == len(states)
        })
    return pd.DataFrame(rows)


def run_analysis(pdb_dir: str = PDB_DIR, selection: str = SELECTION, cutoff: float = CUTOFF,
                 domain: str | None = FILTER_DOMAIN, interaction: str | None = FILTER_INTERACTION) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Analyze contacts for all PDB files in *pdb_dir*.

    Returns a tuple ``(combined, compare_df)`` where ``combined`` is the table of
    contacts for every file and ``compare_df`` indicates which contacts are
    shared across states.
    """
    pdb_files = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
    if not pdb_files:
        raise FileNotFoundError('No PDB files found in directory')

    dfs = []
    for pdb in pdb_files:
        df = analyze_pdb(pdb, selection, cutoff)
        dfs.append(df)
    combined = pd.concat(dfs, ignore_index=True)

    if domain:
        combined = combined[combined['Domain'] == domain]
    if interaction:
        combined = combined[combined['Interaction'] == interaction]

    compare_df = compare_states(dfs)

    print(combined[['State', 'Chain', 'Residue', 'Resnum', 'Domain', 'Interaction', 'Distance']].to_string(index=False))
    print('\nShared vs State-Specific Contacts:')
    print(compare_df.to_string(index=False))

    return combined, compare_df


def main():
    # Wrapper to maintain CLI compatibility. Modify the variables above for notebook use.
    run_analysis()


if __name__ == '__main__':
    main()
