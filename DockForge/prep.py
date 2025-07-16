from openbabel import pybel

def add_hydrogens(input_pdb, output_pdb):
    # Load PDB
    mol = next(pybel.readfile("pdb", input_pdb))

    # Add hydrogens
    mol.addh()

    # Write out the new PDB with hydrogens added
    mol.write("pdb", output_pdb, overwrite=True)
    print(f"Saved molecule with hydrogens to {output_pdb}")
