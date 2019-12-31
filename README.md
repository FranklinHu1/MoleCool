# MoleCool
A Python Implemented Molecule Visualizer
**************************************************************************************
**Welcome to MoleCool!**
**************************************************************************************
**************************************************************************************
**What is MoleCool?**
**************************************************************************************

Molecule is a 2D lewis structure drawing software with 3D visualization capabilities. It uses VSEPR geometries to calculate user-drawn molecular structures in real-time, and renders them as ball and stick models using vpython. Users can then interact with the 3D renders by rotating the molecule, and even save screenshots of the molecule. Additionally, there is a gallery of pre-calculated molecules drawn from the qm7 MATLAB dataset (the dataset comes from http://quantum-machine.org/datasets/) which users can visualize using the gallery visualizer feature. 

**************************************************************************************
**What do each of the modules / folders in the code base do?**
**************************************************************************************

*BackgroundImages (directory):*

A directory that contains all the necessary image files for the project. 	Citations for the images can be found in the code itself as comments

*LewisStructs (directory):*

The main directory for saving lewis structures. Although it is not required that users save their lewis structures here, it is highly recommended for organizational purposes.

*MoleculeGallery (directory):*

Directory containing all 100 of the qm7 molecules which can be viewed. 

*\__init__.py:*

Main file for launching the application

*LewisStructureGen.py:*
	
Contains most of the code for the UI and the canvas for drawing lewis structures on

*Geoms.py:*
	
File for calculating molecule positions based off vector mathematics.

*matlabProcessing.py:*

File for extracting data from qm7, does not need to be run / interacted with in any way

*qm7.mat:*

qm7 dataset, requires MATLAB to open

*TabulationsofCutoffBondLengths2.xlsx:*

Excel file tabulating different inter-atom bond-lengths. Data taken from http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html

*targetOutputBaseNew.py:*

File for translating data from qm7 generated files into actual 3D representations of molecules.

**************************************************************************************
**How do I run the project?**
**************************************************************************************

You've downloaded the code base, so you've already done the hard part! Just open the 
__init__.py file in an editor and run it. 

**************************************************************************************
**What are the dependencies required?**
**************************************************************************************

You will need the following packages to run:

Tkinter
numpy
xlrd
pickle
mendeleev
PIL (Specifically, Image and ImageTk)
vpython

Additionally, the matlabProcessing.py file will also require numpy. Other modules (e.g. os, math, re) are all built-ins that python should already have.

**************************************************************************************
**What keyboard shortcuts exist?**
**************************************************************************************

For drawing lewis structures, the following are keyboard shortcuts:

b for bond placement
i for ring placement
a for atom placement
l for lonepair placement
r for reset

All other commands of the program can be interacted with via the UI and the dropdown 
menus. Also, for a comprehensive guide on how to use the program, be sure to check out the help documentation under the 'help' tab on the dropdown menu



 

   



