# test_package.py

import os
import filecmp

os.system("python ./waterdock/waterdock2.py ./tests/data/example-protein.pdbqt ./tests/data/example-ligand.pdb")


def compare_input_expected():
    assert filecmp.cmp('./tests/data/expected-predictedwaters.pdb', 'predictedwaters.pdb')
