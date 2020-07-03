# WaterDock2.0

![flake8](https://github.com/bigginlab/WaterDock-2.0/workflows/flake8/badge.svg?branch=master) ![mypy](https://github.com/bigginlab/WaterDock-2.0/workflows/mypy/badge.svg)
[![License](https://img.shields.io/github/license/RMeli/pyrmsd?color=%2333BBFF)](https://opensource.org/licenses/MIT)

# This is an updated to Python3 version of WaterDock2.0 https://github.com/bigginlab/WaterDock-2.0
# If you still need the  python2 version, you should visit:  https://github.com/akshay-sridhar/WaterDock2.0

Original Python2 version of WaterDock2.0 written by Akshay Sridhar (https://doi.org/10.1371/journal.pone.0172743)
based on the original WaterDock method developed by Greg Ross (https://doi.org/10.1371/journal.pone.0032036) 
Comments and queuries should be address to philip.biggin@bioch.ox.ac.uk 

This is a Command line interface of WaterDock2.0 for python



# To run

python waterdock2.py proteinfile.pdbqt ligandfile.pdb

This will output a single file called "predictedwaters.pdb", which contains the predictions of the oxygen positions.


Example input files are included as example-protein.pdbqt and example-ligand.pdb


# Dependencies

-- MDAnalysis (version >= 0.13)

-- numpy (any version compatible with your MDAnalysis build)

-- scipy (version does not matter)

-- vina (note on Mac Catalina, you may have to obtain the 64bit version directly from http://vina.scripps.edu/download.html) 


# File Formats

Protein File -- *pdbqt* - the Autodock format of a PDB with charge (Q) and atom type (T)

Ligand File -- *pdb/mol2*  sybyl mol2 format

# Examples/Tests

An example protein (example-protein.pdbqt) and an example ligand are provide (example-ligand.pdb)

python waterdock2.py example-protein.pdbqt example-ligand.pdb 

will generate predictedwaters.pdb which should contain 3 predicted water molecule locations.

You can check the result of the prediction in pymol or compare it directly to the file "expected-predictedwaters.pdb"



