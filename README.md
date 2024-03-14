**Fortran Code for Radial Distribution Function**

Written by EagleLeslie for VASP post-processing using binning method.

**Input files:** 
  - VASP XDATCAR

**Output files:**
  - atoms
  - rad1

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

**File information**

The XDATCAR file can be in any form. This program takes in NPT, NPH, and NVT simulations.

The output file "atoms" contains legend labels for visualization purposes. It is printed in the order in which atoms are read in from the XDATCAR file, and considers all combinations of atom pairs. The file "rad1" contains 
information for plotting the radial distribution function. The first column is the r positions followed by the g(r) values for each atom-pair (following the same order as the "atoms" output file).

**Comile**
To compile the code, enter the command "make". 
Compiler flags in makefile use gfortran. MacOSX library required and uses the flag: GLIB = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib

Once compiled, enter the command "main" to run. Example command line script to run:

      make clean
      make
      ./main

Enter "make clean" to remove any *.o object files.

**Additional files**

The rdf_smooth.py file is a python script for plotting the radial distribution function of all pairs onto one figure.

**Example radial distribution function figure for visualization purposes**

_MgSiO3 at 6000 K and 60 GPa (run for 15 ps)_

![bin_rdf](https://github.com/EagleLeslie/RDF/assets/120432106/68fb778c-354f-4324-8a46-2ec13ef0a743)




