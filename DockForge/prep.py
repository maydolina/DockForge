import subprocess
from openbabel import pybel

def add_hydrogens(input_pdb, output_pdb):
    # load PDB
    mol = next(pybel.readfile("pdb", input_pdb))

    # add H
    mol.addh()

    # write the new PDB with H
    mol.write("pdb", output_pdb, overwrite=True)
    print(f"Hydrogens added and saved to {output_pdb}")

def add_charges(input_pdb, output_pdbqt):
    try:
        subprocess.run([
            "obabel",
            "-ipdb", input_pdb,
            "-opdbqt", "-O", output_pdbqt,
            "-h",  # add H
            "--partialcharge", "gasteiger"
        ], check=True)
        print(f"Charges assigned and saved to {output_pdbqt}")
    except subprocess.CalledProcessError as e:
        print("Error running Open Babel:", e)

# add H, remove waters
def prepare_protein(input_pdb, output_pdbqt):
    subprocess.run([
        "obabel",
        "-ipdb", input_pdb,
        "-opdb", "-O", output_pdbqt,
        "-h",  # add H
        "--delete", "HOH"  # remove water
    ], check=True)
    print(f"Protein prepared and saved to {output_pdbqt}")


