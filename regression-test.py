'''
Simple regression-type test to check output from WaterDock-2.0 is as expected after installation

Runs WaterDock-2.0 and executes python waterdock2.py proteinfile.pdbqt ligandfile.pdb

It will generate predictedwaters.pdb

Test compares this to expected-predictedwaters.pdb


'''

#!/usr/bin/python

import sys
import os

os.system("waterdock2.py example-protein.pdbqt example-ligand.pdb")

f1=open("expected-predictedwaters.pdb","r")
f2=open("predictedwaters.pdb","r")
for line1 in f1:
	for line2 in f2:
		if line1==line2:
			print("Output is as expected - test passed\n")
		else:
			print(line1 + line2)
		break
f1.close()
f2.close()



