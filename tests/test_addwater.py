# test_addwater.py
#  PCB July 2020

import numpy as np
import MDAnalysis as mda
import pytest
import os
import waterdock.addwater as aw

# eg: test_file = os.path.join(THIS_DIR, os.pardir, 'tests/data/test-carbonyl.pdb')
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
test = os.path.join(THIS_DIR, './data/')

test_carbonyl_file = os.path.join(test, 'test-carbonyl.pdb')
test_carboxyl_file = os.path.join(test, 'test-carboxyl.pdb')
test_nitrile_file = os.path.join(test, 'test-nitrile.pdb')
test_secamine_file = os.path.join(test, 'test-secamine.pdb')
test_priamine_file = os.path.join(test, 'test-priamine.pdb')
test_threetertamine_file = os.path.join(test, 'test-threetertamine.pdb')
test_threesecamine_file = os.path.join(test, 'test-threesecamine.pdb')
test_threepriamine_file = os.path.join(test, 'test-threepriamine.pdb')
test_ammonia_file = os.path.join(test, 'test-ammonia.pdb')
test_ether_file = os.path.join(test, 'test-ether.pdb')
test_imine_file = os.path.join(test, 'test-imine.pdb')
test_hydroxyl_file = os.path.join(test, 'test-hydroxyl.pdb')
test_halogen_file = os.path.join(test, 'test-halogen.pdb')


@pytest.fixture
def ligand_carbonyl():
    u = mda.Universe(test_carbonyl_file)
    return u.atoms


@pytest.fixture
def ligand_carboxyl():
    u = mda.Universe(test_carboxyl_file)
    return u.atoms


@pytest.fixture
def ligand_secamine():
    u = mda.Universe(test_secamine_file)
    return u.atoms


@pytest.fixture
def ligand_priamine():
    u = mda.Universe(test_priamine_file)
    return u.atoms


@pytest.fixture
def ligand_threetertamine():
    u = mda.Universe(test_threetertamine_file)
    return u.atoms


@pytest.fixture
def ligand_threesecamine():
    u = mda.Universe(test_threesecamine_file)
    return u.atoms


@pytest.fixture
def ligand_threepriamine():
    u = mda.Universe(test_threepriamine_file)
    return u.atoms


@pytest.fixture
def ligand_nitrile():
    u = mda.Universe(test_nitrile_file)
    return u.atoms


@pytest.fixture
def ligand_ammonia():
    u = mda.Universe(test_ammonia_file)
    return u.atoms


@pytest.fixture
def ligand_ether():
    u = mda.Universe(test_ether_file)
    return u.atoms


@pytest.fixture
def ligand_imine():
    u = mda.Universe(test_imine_file)
    return u.atoms


@pytest.fixture
def ligand_hydroxyl():
    u = mda.Universe(test_hydroxyl_file)
    return u.atoms


@pytest.fixture
def ligand_halogen():
    u = mda.Universe(test_halogen_file)
    return u.atoms


def test_unitvector():
    assert aw.unitvector(3) == 1


# test against expected-waters.pdb as the gold standard - below is a bit clumsy bit works.
def test_writewaterfile():
    watercoods = np.zeros((1, 3))
    watercoods[0, 0] = 20.673
    watercoods[0, 1] = 13.203
    watercoods[0, 2] = 42.688
    filename = "test-expected-waters.pdb"
    aw.writewaterfile(filename, watercoods)

    with open("tests/data/expected-predictedwaters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-expected-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


def test_carbonyl(ligand_carbonyl):
    bond_dist = 1.7
    atype = aw.carbonylorcarboxyl(ligand_carbonyl, 3, bond_dist)  # 3 is 4th atom (an oxygen in the carbonyl test)
    assert atype == 'carbonyl'


def test_carboxyl(ligand_carboxyl):
    bond_dist = 1.7
    atype = aw.carbonylorcarboxyl(ligand_carboxyl, 2, bond_dist)  # 2 is the 3rd atom (an oxygen in the carboxyl test)
    assert atype == 'carboxyl'


def test_matefinder(ligand_carboxyl):
    mates = aw.matefinder(ligand_carboxyl,  3)  # 3 is the 2nd oxygen - the hydroxyl oxygen
    nummates = np.size(mates)
    assert nummates == 2


def test_Otypefinder(ligand_carboxyl):
    bond_dist = 1.7
    otype = aw.Otypefinder(ligand_carboxyl, 3, bond_dist)  # 3 is the 2nd oxygen - the hydroxyl oxygen
    assert otype == 'hydroxyl'


def test_Ntypefinder(ligand_nitrile):
    bond_dist = 1.7
    ntype = aw.Ntypefinder(ligand_nitrile, 1, bond_dist)  # 1 is the n in the nitrile
    assert ntype == 'nitrile'


def test_carbonylwater(ligand_carbonyl):
    bond_dist = 1.7
    watercoords = aw.carbonylwaters(ligand_carbonyl, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-carbonyl-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-carbonyl-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-carbonyl-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


def test_carboxylwater(ligand_carbonyl):
    bond_dist = 1.7
    watercoords = aw.carboxylwaters(ligand_carbonyl, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-carboxyl-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-carboxyl-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-carboxyl-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


def test_secaminewater(ligand_secamine):
    bond_dist = 1.7
    watercoords = aw.secaminewater(ligand_secamine, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-secamine-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-secamine-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-secamine-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


def test_priaminewater(ligand_priamine):
    bond_dist = 1.7
    watercoords = aw.priaminewater(ligand_priamine, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-priamine-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-priamine-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-priamine-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


def test_threetertaminewater(ligand_threetertamine):
    bond_dist = 1.7
    watercoords = aw.threetertaminewater(ligand_threetertamine, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-threetertamine-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-threetertamine-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-threetertamine-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


def test_threesecaminewater(ligand_threesecamine):
    bond_dist = 1.7
    watercoords = aw.threesecaminewater(ligand_threesecamine, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-threesecamine-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-threesecamine-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-threesecamine-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


def test_threepriaminewater(ligand_threepriamine):
    bond_dist = 1.7
    watercoords = aw.threepriaminewater(ligand_threepriamine, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-threepriamine-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-threepriamine-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-threepriamine-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


def test_ammonia(ligand_ammonia):
    bond_dist = 1.7
    watercoords = aw.ammoniawater(ligand_ammonia, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-ammonia-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-ammonia-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-ammonia-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test

# large atom test should go here - omitted at the minute


def test_ether(ligand_ether):
    bond_dist = 1.7
    watercoords = aw.etherwater(ligand_ether, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-ether-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-ether-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-ether-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


def test_imine(ligand_imine):
    bond_dist = 1.7
    watercoords = aw.iminewater(ligand_imine, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-imine-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-imine-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-imine-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


def test_hydroxyl(ligand_hydroxyl):
    bond_dist = 1.7
    watercoords = aw.hydroxylwater(ligand_hydroxyl, 1, bond_dist)
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    filename = "test-hydroxyl-waters.pdb"
    aw.writewaterfile(filename, watercoords)
    with open("tests/data/test-expected-hydroxyl-waters.pdb") as f:
        first_line_expected = f.readline().rstrip()

    with open("test-hydroxyl-waters.pdb") as f:
        first_line_test = f.readline().rstrip()

    assert first_line_expected == first_line_test


#  test halogen in a different way - just in case.
def test_halogenwater(ligand_halogen):
    bond_dist = 1.7
    watercoords = aw.halogenwater(ligand_halogen, 0, bond_dist)  # 0 is Cl
    x = watercoords[0, 0]
    y = watercoords[0, 1]
    z = watercoords[0, 2]
    assert x == pytest.approx(-8.800)  # floats are hard to assert accurately
    assert y == pytest.approx(2.862)   # numbers here obtain by writing out a file first
    assert z == pytest.approx(0.000)

@pytest.fixture(scope='session')
def clear_files_teardown():
    yield None
    os.system("rm test*.pdb")
