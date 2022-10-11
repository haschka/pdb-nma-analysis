# pdb-nma-analysis

(c) Thomas Haschka 2009-2022

A tool to perform quick normal mode analysis on a PDB file
according to the method first proposed by
[Tirion et al.](https://doi.org/10.1006/jmbi.1993.1135).

Please cite: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6888913.svg)](https://doi.org/10.5281/zenodo.6888913)

This tool now has a youtube tutorial: [https://youtu.be/jUBWMLQYkKk](https://youtu.be/jUBWMLQYkKk)

### Compiling the Program

The program was almost completely ( by exercise ) 
written using the SSE2 + SSSE3 intrinsics, and as 
such requires an intel/amd processor supporting these 
instructions to work. Anything built after 2008 basically
should do it.

The program further needs the lapack libraries. On 
a contemporary Debian Based system you should be able to obtain
all requirements using:
```
sudo apt-get install build-essential liblapack-dev
```
You can then compile the program:
```
gcc -O2 -march=native nma-double.c -o nma -lm -llapack
```
obtaining the nma binary.

### Performing Normal Mode Analysis

The nma tool is pretty self explanatory. Its arguments are:
```
file: pdb file 
   k: spring constant in atom - atom interaction potential 
 cut: Cutoff of atom - atom interaction in Angstr. 
nmax: highest mode to be visualized. nmax has to satisfy 
       7 < nmax < 100
```

in example if you download the PDB file 
[1uqq](https://files.rcsb.org/download/1QUU.pdb)
you may perform nma analysis in typing:
```
./nma 1quu.pdb pdb 0.5 15 9
```
which should yield the files:
```
7-mode.xyz
8-mode.xyz
9-mode.xyz
eigenvalues
```
The eigenvalues give you an idea of how much each
mode should get occupied as the protein gets excited.
Each mode can be visualized for instance with 
[VMD](https://www.ks.uiuc.edu/Research/vmd/)
using the Van der Waals representation
in this program you might get nice movies of each mode. 

### History
The program was in its original form written in 2009 as part of Thomas Haschka's 
Master Theses
