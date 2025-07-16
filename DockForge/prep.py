from rdkit import Chem

def add_hydrogens(input_pdb, output_pdb):
    # Load PDB
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False, sanitize=False)
    if mol is None:
        print("Failed to load molecule.")
        return

    # Add H
    mol_H = Chem.AddHs(mol)

    # Optimize geometry
    try:
        AllChem.EmbedMolecule(mol_H)
        AllChem.UFFOptimizeMolecule(mol_H)
    except Exception as e:
        print(f"Embedding/optimization failed: {e}")

    # Save molecule as PDB
    Chem.MolToPDBFile(mol_H, output_pdb)
    print(f"Saved molecule with hydrogens to {output_pdb}")