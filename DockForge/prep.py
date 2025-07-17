import subprocess
import os
from pathlib import Path
from openbabel import pybel
from pdbfixer import PDBFixer
from openmm.app import PDBFile


def check_openbabel():
    try:
        subprocess.run(["obabel", "-V"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        print("\033[91mERROR: Open Babel (obabel) not found. Please install it and ensure it's in PATH.")
        exit(1)


def add_hydrogens(input_pdb, output_pdb=None):
    # load PDB
    mol = next(pybel.readfile("pdb", input_pdb))

    # default value for output_pdb
    input_path = Path(input_pdb)
    if not output_pdb:
        output_pdb = input_path.with_name(f"{input_path.stem}_with_h.pdb")

    # add H
    mol.addh()

    # save the new PDB with H
    mol.write("pdb", str(output_pdb), overwrite=True)



def add_charges(input_pdb, output_pdbqt=None):

    # default value for output_pdbqt
    input_path = Path(input_pdb)
    if not output_pdbqt:
        output_pdbqt = input_path.with_name(f"{input_path.stem}_with_charges.pdbqt")

    output_pdbqt = Path(output_pdbqt)
    if output_pdbqt.exists():
        print(f"\033[93mWARNING: Overwriting existing file {output_pdbqt.name}")

    try:
        subprocess.run([
            "obabel",
            "-ipdb", input_pdb,
            "-opdbqt", "-O", output_pdbqt,
            "-h",  # add H
            "--partialcharge", "gasteiger"
        ], check=True)
        # saves file as PDBQT

    except subprocess.CalledProcessError as e:
        print("\033[91mERROR running Open Babel:", e)




# add H, remove waters
def simple_prepare_protein(input_pdb, output_pdbqt=None):

    mol = next(pybel.readfile("pdb", input_pdb))
    print(f"\033[97m{len(mol.atoms)} atoms and {mol.OBMol.NumBonds()} bonds detected")

    # default value for output_pdbqt
    input_path = Path(input_pdb)
    if not output_pdbqt:
        output_pdbqt = input_path.with_name(f"{input_path.stem}_prepared.pdbqt")

    if output_pdbqt.exists():
        print(f"\033[93mWARNING: Overwriting existing file {output_pdbqt.name}")

    subprocess.run([
        "obabel",
        "-ipdb", input_pdb,
        "-opdb", "-O", output_pdbqt,
        "-h",  # add H
        "--delete", "HOH"  # remove waters
    ], check=True)
    # saves file as PDBQT

    print("\033[97mHydrogens added")
    print("\033[97mWaters removed")

    print(f"\033[92mProtein prepared and saved as: {output_pdbqt}")



# add H at ph 7, add charges, remove water&ligands, add missing atoms&residudes
def advanced_prepare_protein(input_pdb, output_pdbqt=None, ph=7.0):

    mol = next(pybel.readfile("pdb", input_pdb))
    print(f"\033[97m{len(mol.atoms)} atoms and {mol.OBMol.NumBonds()} bonds detected")

    # default value for output_pdbqt
    input_path = Path(input_pdb)
    if not output_pdbqt:
        output_pdbqt = input_path.with_name(f"{input_path.stem}_prepared.pdbqt")

    if output_pdbqt.exists():
        print(f"\033[93mWARNING: Overwriting existing file {output_pdbqt.name}")

    temp_pdb = f"temp_{input_path.stem}.pdb"  # temporary intermediate

    fixer = PDBFixer(filename=input_pdb)

    # find missing residudes
    fixer.findMissingResidues()
    print("\033[97mMissing residues added")

    # find & add missing atoms
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    print("\033[97mMissing atoms added")

    # remove any ligands & water (non-protein)
    fixer.removeHeterogens(keepWater=False)
    print("\033[97mAny ligand & water (non-protein) removed")

    # add H at {ph}
    fixer.addMissingHydrogens(pH=ph)
    print(f"\033[97mHydrogens added at pH {ph}")

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

    print(f"\033[92mProtein prepared and saved as: {output_pdbqt}")



def prepare_ligand(input_pdb, output_pdbqt=None):

    # default value for output_pdbqt
    input_path = Path(input_pdb)
    if not output_pdbqt:
        output_pdbqt = input_path.with_name(f"{input_path.stem}_prepared.pdbqt")

    temp_pdb = f"temp_{input_path.stem}.pdb" # temporary intermediate

    try:
        print(f"\033[97mPreparing ligand {input_pdb}......")

        # add H
        add_hydrogens(input_pdb, temp_pdb)
        print("\033[97mHydrogens added")

        # assign charges and save as PDBQT
        add_charges(temp_pdb, output_pdbqt)
        print("\033[97mCharges added")

    # remove temporary PDB intermediate
    finally:
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)

    print(f"\033[92mLigand prepared and saved as: {output_pdbqt}")



def convert_to_pdb(input_file, output_pdb=None):

    input_path = Path(input_file)
    input_type = input_path.suffix.lower().lstrip(".")

    # default value for output_pdb
    if not output_pdb:
        output_pdb = input_path.with_suffix(".pdb")

    if output_pdb.exists():
        print(f"\033[93mWARNING: Overwriting existing file {output_pdb.name}")

    try:
        subprocess.run([
            "obabel",
            f"-i{input_type}", str(input_path),
            "-opdb", "-O", str(output_pdb)
        ], check=True)

        print(f"\033[92mConverted {input_file} to {output_pdb}")
        return str(output_pdb)

    except subprocess.CalledProcessError as e:
        print(f"\033[91mERROR running Open Babel: {e}")
        return None



def batch_convert_to_pdb(input_dir, output_dir):

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    supported_types = {".sdf", ".mol2", ".smi", ".pdb", ".xyz", ".inchi"}
    files = [f for f in input_dir.iterdir() if f.suffix.lower() in supported_types]

    if not files:
        print("\033[91mERROR: No supported input files found")
        return

    print(f"\033[97mConverting {len(files)} files to PDB......\n")

    for file in files:
        input_type = file.suffix.lower().lstrip(".")
        output_file = output_dir / f"{file.stem}.pdb"

        if output_file.exists():
            print(f"\033[93mWARNING: Overwriting existing file {output_file.name}")

        try:
            subprocess.run([
                "obabel",
                f"-i{input_type}", str(file),
                "-opdb", "-O", str(output_file)
            ], check=True)

        except subprocess.CalledProcessError as e:
            print(f"\033[91mERROR running Open Babel: {e}")

    print("\n\033[94mBatch conversion complete.")



def batch_prepare_ligands(input_dir, output_dir, keep_original_name=False):

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    valid_extensions = {".pdb", ".mol2", ".sdf", ".smi"}
    ligand_files = [f for f in input_dir.iterdir() if f.suffix.lower() in valid_extensions]

    if not ligand_files:
        print("\033[91mERROR: No valid ligand files found")
        return

    print(f"\033[97mFound {len(ligand_files)} ligands in {input_dir}. Starting preparation......\n")

    for file in ligand_files:

        if not keep_original_name:
            output_pdbqt = output_dir / f"{file.stem}_prepared.pdbqt"
        elif keep_original_name:
            output_pdbqt = output_dir / f"{file.stem}.pdbqt"

        prepare_ligand(str(file), str(output_pdbqt))


    print("\033[94mBatch ligand preparation complete")



def print_molecule_info(file_path):

    mol = next(pybel.readfile("pdb", file_path))

    num_atoms = len(mol.atoms)
    num_bonds = mol.OBMol.NumBonds()
    formula = mol.formula
    mol_weight = mol.molwt

    print(f"\033[97m---------- Molecule: {file_path} ----------")
    print(f"\033[97mAtoms: {num_atoms}")
    print(f"\033[97mBonds: {num_bonds}")
    print(f"\033[97mFormula: {formula}")
    print(f"\033[97mMolecular Weight: {mol_weight:.2f} g/mol")