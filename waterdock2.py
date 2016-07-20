'''
A command line interface for WaterDock 2

HOW TO USE:
python waterdock2.py proteinfile.pdbqt ligandfile.pdb

OUTPUT FILE -- predictedwaters.pdb

DEPENDENCIES:

**MDAnalysis -- (version >= 0.13)
**numpy -- (any version compatible with the MDAnalysis build)
**scipy

'''
import sys
import os

def waterfile():
	f1 = open('water.pdbqt', 'w')

	f1.write('REMARK  The pdbqt file for using water as a ligand\n')
	f1.write('ROOT\n')
	f1.write('ATOM      1  OW  HOH   231       0.950  11.375  16.494  1.00  0.00    -0.411 OA\n')
	f1.write('ATOM      2  HW1 HOH   231       1.766  11.375  17.071  1.00  0.00     0.205 HD\n')
	f1.write('ATOM      3  HW2 HOH   231       0.134  11.375  17.071  1.00  0.00     0.205 HD\n')
	f1.write('ENDROOT\n')
	f1.write('TORSDOF 0')

	f1.close()


proteinfile, ligandfile = sys.argv[1], sys.argv[2]

if int(os.path.isfile(proteinfile)) == 0:
	sys.exit('Protein File does not exist')

if int(os.path.isfile(ligandfile)) == 0:
	sys.exit('Ligand File does not exist')

ext1 = proteinfile[-5:]
ext2 = ligandfile[-4:]

if ext1 != 'pdbqt':
	sys.exit('Protein File must be in pdbqt format')

if ext2 != '.pdb' and ext2 != 'mol2':
	sys.exit('Ligand File must be in pdb/mol2 format')

if os.path.isfile('predictedwaters.pdb'):
	os.rename('predictedwaters.pdb','predictedwaters1.pdb')

waterfile()

import addwater
addwater.main(ligandfile)

import dockcheck
dockcheck.main(proteinfile, ligandfile)

os.remove('waterdetails.txt')
os.remove('placedwaters.pdb')
os.remove('water.pdbqt')