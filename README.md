Variational Monte-Carlo algorithm for the SU(N) Heisenberg model
===

Introduction
---

The code provided here is what I wrote during my PhD thesis at the École
Polythecnique Fédérale de Lausanne (EPFL) in Frédéric Mila's group. It is aimed
to compute variational energies of the SU(*N*) Heisenberg model in any fully
antisymmetric irreducible representation (one or more particle per site) on
various lattices. To make transparent what the code can do, let us very briefly
explain what this model is about. It describes the Mott insulating phase of
fermionic ultracold atom gas trapped in optical lattices. This field belongs
to quantum physics. The fermionic atoms have *N* internal degrees of freedom
that will be referred to as colors. On each of the *n* sites of the lattice,
there are *m* particles *(N>m)* that can be exchanged with particles on
neighbouring sites.

The code computes the variational energy of a fermionic Gutzwiller projected
wave function thanks to a Monte-Carlo metropolis sampling. It does not
currently support the bosonic variational wave functions based on Jastrow
factors but the structure of the code is such that this type of wave function
can be added with little effort (it has already been tested in this code for
SU(2) and first neighbour Jastrow factors, so the generalisations to *N>2* and
long range Jastrow are possible).

The following scientific articles exemplify the abilities of this code:

+ [Stabilization of the chiral phase of the SU*(6m)* Heisenberg model on the honeycomb lattice with *m* particles per site for *m* larger than 1](https://arxiv.org/abs/1607.05227)
+ [Variational Monte-Carlo investigation of SU*(N)* Heisenberg chains](https://arxiv.org/abs/1502.01895)

and other (unpublished) results can be found in my thesis:

+ [doi:10.5075/epfl-thesis-7309](https://infoscience.epfl.ch/record/222440)

Installation
---

Get all the sources:

	git clone --recursive https://github.com/jduf/VMC-SUN-Heisenberg.git

This should download the main sources and a submodule containing libraries
also available via:

	git clone https://github.com/jduf/lib.git

Enter the main directory:

	cd VMC-SUN-Heisenberg

and simply type:

	make

If you choose to install the submodules elsewhere, just set the variable
"JDLIB" in the makefile to the path where they have been downloaded.

Note that you can choose which executable to create by un-/-commenting one of
the first line of the makefile or setting the "EXEC" variable to one of the
following:

+ mc
+ min
+ study
+ check
+ mcbi

Usage
---

To run an executable, the structure of the command is the following

	./executable -args -t:name1 val1 -t:name2 val2 -t:name3 val3 ...

where

+ args are arguments without options for instance **min** handles:
	- d: displays the results of the simulation in the web browser (can be set
	  in **lib/lib/config.mk**)
	- p: prints the results of the simulation in a pdf file
	- norun: does not do run a Monte-Carlo simulation (useful to display or
	  print the results of a simulation that already exists)
+ t is the c++ type of variable, all possible choices are:
	- s for string
	- d for double
	- u for unsigned int
	- i for int
+ name is the name of the argument, some examples are:
	- wf for the type of wave function
	- *n* for the number of sites
	- *N* for the number of colors
	- *m* for the number of particles per sites
+ val is the value of the argument

Here is an example of a simple command:

	./mc -d -u:tmax 10 -s:wf square-mu -u:N 4 -u:m 2 -u:n 20 -i:bc 1 -d:mu 0.1 -u:obs 1,3

which runs a simulation for 10 seconds on the wave function named square-mu for
SU(4) with *m=2* particles per site and *n=20* sites on a periodic lattice and
an on site chemical potential of *\mu=0.1* and measures the bond energy and
long range-correlations (observables 1 and 3). The results of the simulation
are then saved in a directory tree with the extension **.jdbin** and displayed
in a web browser. The **.jdbin** files can be reloaded for further use, for
instance the following command print the results in a pdf file.

    ./mc -p -s:sim path/to/results.jdbin -norun

Abilities
---

The code currently supports:

+ any value of *N>1*
+ any number of particle per site *m>N*
+ any number of site *n*
+ any simple boundary condition (periodic, antiperiodic, none). Note that if a
  link crosses the boundaries an even number of time, it is considered as if no
  boundaries were crossed (this is relevant for some particular geometries like
  a tilted cluster on the square lattice for which the links at the corners
  might cross twice the boundaries)
+ any number of particle of each color (having more particles of one color is
  allowed)
+ first neighbour Heisenberg model (easily generalisable to further
  neighbours)
+ inhomogeneous coupling between sites

The message returned when an error occurs should be self explanatory, for
instance if *N=4 m=1* and *n=13* the error message will be

	N, M, m and n are incompatible

because if you want to have *N=4* colors and *m=1* particle per site, the you
need a number of sites which is a multiple of 4. Moreover, the number of site
that you set has to be coherent with the lattice geometry you chose (chain,
ladder, square, triangle,...). For instance, on the square lattice, only the
clusters with *n=pp+qq* (*p* and *q* integers) are allowed therefore *n=13=9+4*
would be allowed.

The values that the option *obs* take defines what observables are measured
(the energy per site is always measured):

0. all possible observable
1. bond energy
2. long range correlations and associated structure factor
3. color occupation
4. variance of the energy (development phase)

Beside performing a VMC simulation

+ can visualise the lattice with unit cell, site numbers, basis vectors,
  hopping amplitude, fluxes

The executable **min** intends to perform a minimisation of variational
parameters and **study** allows an analysis of these results in a web browser.

Note on the output files
---

In order to save disk space for large simulations, the results are saved in a
binary file with the extension **.jdbin**. These files contains a footer that can
be read with the utilities in **lib/jdtools/**. For instance :

    jdhtml path/to/results.jdbin

will display the results in the web browser

Requirement
---

+ c++ 4.9.2 (may compile with older versions)

All the following softwares are not essential be should be present to use this
code at its best potential. Indeed, the visualisation of the results, already
implemented in the code, relies on:

+ gnuplot
+ latex
+ ReStructuedText
+ rst2html
+ rst2latex
+ pdflatex
+ latex
+ dvipdf
+ pdfcrop
+ gs
+ firefox

All these software can be set and personalised in the **lib/lib/config.mk** and
**lib/lib/Linux.cpp** files.  For instance if *gs* is not available to convert
images, only these two files need to be modified.

The code is commented and the documentation can be generated via Doxyfile using

	make ref
