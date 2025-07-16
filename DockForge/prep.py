import subprocess
import os
from openbabel import pybel
from pdbfixer import PDBFixer
from openmm.app import PDBFile

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
def simple_prepare_protein(input_pdb, output_pdbqt):

    subprocess.run([
        "obabel",
        "-ipdb", input_pdb,
        "-opdb", "-O", output_pdbqt,
        "-h",  # add H
        "--delete", "HOH"  # remove water
    ], check=True)

    print(f"Protein prepared and saved to {output_pdbqt}")


# add H at ph 7, add charges, remove water&ligands, add missing atoms&residudes
def advanced_prepare_protein(input_pdb, output_pdbqt, ph=7.0):

    output_pdb = "temp_prepared.pdb"  # temporary intermediate

    fixer = PDBFixer(filename=input_pdb)

    print("Finding missing residues...")
    fixer.findMissingResidues()

    print("Finding missing atoms...")
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    print("Removing water & ligands...")
    fixer.removeHeterogens(keepWater=False)

    print(f"Adding hydrogens at pH {ph}...")
    fixer.addMissingHydrogens(pH=ph)

    print(f"Writing intermediate PDB to {output_pdb}...")
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"Converting {output_pdb} to PDBQT {output_pdbqt} with charges...")
    subprocess.run([
        "obabel",
        "-ipdb", output_pdb,
        "-opdbqt", "-O", output_pdbqt,
        "-h",
        "--partialcharge", "gasteiger"
    ], check=True)

    # remove temporary PDB file
    if os.path.exists(output_pdb):
        os.remove(output_pdb)
        print(f"Temporary file {output_pdb} deleted.")

    print(f"Protein prepared and saved to {output_pdbqt}")


def prepare_ligand(input_pdb, output_pdbqt):

    output_H = "temp_output.pdb" #temporary intermediate

    try:
        print(f"Preparing ligand {input_pdb}...")

        # add H
        add_hydrogens(input_pdb, output_H)

        # assign charges and save as PDBQT
        add_charges(output_H, output_pdbqt)

    finally:
        if os.path.exists(output_H):
            os.remove(output_H)
            print(f"Deleted temporary file: {output_H}")

    print(f"Ligand prepared and saved as {output_pdbqt}")

