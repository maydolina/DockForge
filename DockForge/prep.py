import subprocess
import os
from pathlib import Path
from openbabel import pybel
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def add_hydrogens(input_pdb, output_pdb):
    # load PDB
    mol = next(pybel.readfile("pdb", input_pdb))

    # add H
    mol.addh()

    # save the new PDB with H
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
        # saves file as PDBQT

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
        "--delete", "HOH"  # remove waters
    ], check=True)
    # saves file as PDBQT

    print(f"Protein prepared and saved to {output_pdbqt}")



# add H at ph 7, add charges, remove water&ligands, add missing atoms&residudes
def advanced_prepare_protein(input_pdb, output_pdbqt, ph=7.0):

    temp_pdb = "temp_prepared.pdb"  # temporary intermediate

    fixer = PDBFixer(filename=input_pdb)

    # find missing residues
    fixer.findMissingResidues()

    # find & add missing atoms
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # remove any ligands & water (non-protein)
    fixer.removeHeterogens(keepWater=False)

    # add H at {ph}
    fixer.addMissingHydrogens(pH=ph)

    # save intermediate as PDB to temp_pdb
    with open(temp_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    # convert intermediate PDB to PDBQT, to output_pdbqt
    subprocess.run([
        "obabel",
        "-ipdb", temp_pdb,
        "-opdbqt", "-O", output_pdbqt,
        "-h",
        "--partialcharge", "gasteiger"
    ], check=True)

    # remove temp_pdb
    if os.path.exists(temp_pdb):
        os.remove(temp_pdb)
        print(f"Temporary file {temp_pdb} deleted.")

    print(f"Protein prepared and saved to {output_pdbqt}")



def prepare_ligand(input_pdb, output_pdbqt):

    temp_pdb = "temp_output.pdb" # temporary intermediate

    try:
        print(f"Preparing ligand {input_pdb}...")

        # add H
        add_hydrogens(input_pdb, temp_pdb)

        # assign charges and save as PDBQT
        add_charges(temp_pdb, output_pdbqt)

    # remove temporary PDB intermediate
    finally:
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)
            print(f"Deleted temporary file: {temp_pdb}")

    print(f"Ligand prepared and saved as {output_pdbqt}")


import subprocess
from pathlib import Path



def convert_to_pdb(input_file, output_file=None):

    input_path = Path(input_file)
    input_type = input_path.suffix.lower().lstrip(".")

    if not output_file:
        output_file = input_path.with_suffix(".pdb")

    try:
        subprocess.run([
            "obabel",
            f"-i{input_type}", str(input_path),
            "-opdb", "-O", str(output_file)
        ], check=True)
        print(f"Converted {input_file} to {output_file}")
        return str(output_file)
    except subprocess.CalledProcessError as e:
        print(f"Open Babel failed: {e}")
        return None



def batch_convert_to_pdb(input_dir, output_dir):

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    supported_types = {".sdf", ".mol2", ".smi", ".pdb", ".xyz", ".inchi"}
    files = [f for f in input_dir.iterdir() if f.suffix.lower() in supported_types]

    if not files:
        print("No supported input files found")
        return

    print(f"Converting {len(files)} files to PDB......\n")

    for file in files:
        input_type = file.suffix.lower().lstrip(".")
        output_file = output_dir / (file.stem + ".pdb")

        try:
            subprocess.run([
                "obabel",
                f"-i{input_type}", str(file),
                "-opdb", "-O", str(output_file)
            ], check=True)
            print(f"Converted: {file.name} to {output_file.name}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to convert {file.name}: {e}")

    print("\n Batch conversion complete.")



def batch_prepare_ligands(input_dir, output_dir):

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    valid_extensions = {".pdb", ".mol2", ".sdf", ".smi"}
    ligand_files = [f for f in input_dir.iterdir() if f.suffix.lower() in valid_extensions]

    if not ligand_files:
        print("No valid ligand files found")
        return

    print(f"Found {len(ligand_files)} ligands in {input_dir}. Starting preparation......\n")

    for file in ligand_files:
        output_file = output_dir / (file.stem + ".pdbqt")

        try:
            print(f"Preparing: {file.name}")
            prepare_ligand(str(file), str(output_file))
            print(f"Saved: {output_file.name}\n")
        except Exception as e:
            print(f"Failed to process {file.name}: {e}\n")

    print("Batch ligand preparation complete")