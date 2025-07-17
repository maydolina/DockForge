# DockForge
[View DockForge on PyPI](https://pypi.org/project/DockForge/)  
  
A Python package for preparing proteins and ligands for molecular docking.  
Includes: hydrogen addition, charge assignment, format conversion, cleanup.  

## Usage
The package can be installed via:  
```pip install DockForge```  

## Dependencies
Requires **Open Babel** to be installed separately. Install via:  
```brew install open-babel``` or ```conda install -c conda-forge openbabel```  
See [Open Babel installation instructions](https://openbabel.org/docs/Installation/install.html).
  
Also requires **PDBfixer** to be installed. Install via:  
```conda install -c conda-forge pdbfixer```  
See on [conda website](https://anaconda.org/conda-forge/pdbfixer).  
