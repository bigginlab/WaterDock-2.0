**This file is a more detailed version of README.md.**
N.B. The command line commands in this text file work for Mac but may not work for Windows and Linux. You might need to look up the equivalent commands for these operating systems.

HOW TO SET UP WATERDOCK-2.0 (DETAILED)

1.   Download the zip file or clone using GitHub.


2.   Have Python 3.7 (or the latest version of Python 3 installed onto your computer).


3.   Ensure that you have the MDAnalysis, numpy and scipy packages installed using either pip or conda.
      - Look up the correct commands for installing each package using Google\


4.   Install Autodock Vina using this website "http://vina.scripps.edu/download.html" to any folder that you like.
     - It installs as a zip file, so open it to export the folder into the same directory as the zip file.


5.   Add Autodock Vina program (found within the bin folder of the folder you have downloaded) to your PATH.
      - You can view your PATH using the command `echo $PATH` in Terminal.
      - To add to your path, use this command `sudo nano /etc/paths`.
      - Type in the password for your computer.
      - Press the down key until you reach the bottom of the text, which should include other directories like `/usr/bin`. 
      - Type in the directory that should end in `.../vina/bin`.
      - Press Control-X.
      - Press Y and then Enter.
      - Restart the Terminal.


6.   Navigate to your WaterDock-2.0 folder using the `cd` command and its directory and run the example.
      - Use the command `python3 waterdock2.py example-protein.pdbqt example-ligand.pdb`.
      - This produces a file called **predicted waters.pdb** in the folder which can be opened with Pymol or just with a Text editor. This should contain 3 waters.


If you get that result, you have WaterDock-2.0 working!


HOW TO RUN WATERDOCK-2.0 ON YOUR OWN PROTEIN-LIGAND COMPLEXES

1.   Download the protein-ligand complex pdf file from the Protein Databank  (https://www.rcsb.org/)
     - 1g9v is being used as an example for the command line commands.


2.   Separate the protein from the protein complex to produce a pdb file that only contains the protein using this command line command (when in its directory).
      - `more 1g9v.pdb |grep ATOM > 1g9v_prot.pdb`
      - Check for multiple conformations by looking out for A and B versions of the same atoms and deleting one version, whilst renaming the other without the letter in front of the atom name.


3.   Separate the ligand from the protein complex to produces its own pdb file that only contains this ligand again using a command line command.
      - `more 1g9v.pdb |grep HETATM |grep *name of ligand* > 1g9v_lig.pdb`
      - *name of ligand* can be found by looking through the original pdb to find under the column that denotes amino acids or names of molecules and finding the name used for the ligand e.g. RQ3.
      - Check for multiple conformations by looking out for A and B versions of the same atoms and deleting one version, whilst renaming the other without the letter in front of the atom name.
      - Also check if you want to use all the ligands if multiple of the same ligand bind to the protein complex or just one, delete the extras if required from the texfile produced.


4.   Ensure that you have ADFR suite installed, if not, install from this website - (https://ccsb.scripps.edu/adfr/downloads/).


5.   Reduce the protein and ligand files using this command line command.
      - `$ADFRsuite_INSTALL_DIR/bin/reduce ~/Download/1g9v_prot.pdb > 1g9v_protH.pdb`
      - `$ADFRsuite_INSTALL_DIR/bin/reduce ~/Download/1g9v_lig.pdb > 1g9v_ligH.pdb`


6.   Convert the reduced protein file to a pdbqt file using the prepare_receptor.
      - `$ADFRsuite_INSTALL_DIR/bin/prepare_receptor -r 1g9v_protH.pdb -o 1g9v_protH.pdbqt`
      - Open-babel normally could be used for this step but the pdbqt file it produces does not have the correct format with the WaterDock-2.0 program, so it is recommended not to use it.


7.   Copy the reduced ligand pdb file and reduced protein pdbqt file to your WaterDock-2.0 program and run the program as stated in the README.md.
      - `python3 waterdock2.py 1g9v_protH.pdbqt 1g9v_ligH.pdb`


8.   Results will be in **predicted_waters.pdb** in the WaterDock-2.0 folder.
