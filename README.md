Variational Monte-Carlo algorithm for the SU(N) Heisenberg model
===

Introduction
---

The code provided here is what I wrote during my PhD thesis at the École
Polythecnique Fédérale de Lausanne (EPFL) in Frédéric Mila's group. It is aimed
to compute variational energies of the SU($N$) Heisenberg model in any fully
antisymmetric irreducible representation (one or more particle per site) on
various lattices. To make transparent what the code can do, let us very briefly
explain what this model is about. It describes the Mott insulating phase of
fermionic ultracold atom  gas trapped in optical lattices. This field belongs
to quantum physics . The fermionic atoms have $N$ internal degrees of freedom
that will be referred to as colours. On each of the $n$ sites of the lattice,
there are $m$  particles  $(N>m)$ that can be exchanged with particles on
neighbouring sites.

The code computes the variational energy of a fermionic Gutzwiller projected
wave function thanks to a Monte-Carlo metropolis sampling. It does not
currently support the bosonic variational wave functions based on Jastrow
factors but the structure of the code is such that this type of wave function
can be added with little effort (it has already been tested in this code for
SU(2) and first neighbour Jastrow factors, so the generalisations to $N>2$ and
long range Jastrow are possible).


The following scientific articles exemplify the abilities of this code:


+ [Stabilization of the chiral phase of the SU $(6m)$ Heisenberg model on the honeycomb lattice with $m$ particles per site for $m$  larger than 1](https://arxiv.org/abs/1607.05227) 
+ [Variational Monte-Carlo investigation of SU$(N)$ Heisenberg chains](https://arxiv.org/abs/1502.01895)


and other (unpublished) results can be found in my thesis:

+ [doi:10.5075/epfl-thesis-7309](https://infoscience.epfl.ch/record/222440)

Installation
---

Simply  enter the directory and type 

	make
	
Note that it will create different executables

+ mc
+ min
+ study
+ (some other executable may be created by created by un-/-commenting one of
  the first line of the makefile)

Usage
---

To run an executable the structure is the following

	./executable -args -t:name1 val1  -t:name2 val2  -t:name3 val3 ... 
	
where 

+ args are arguments without options for instance:
	- norun: do not do run a Monte-Carlo simulation
	- d: display the result of the simulation in the browser
	
+  t is the c++ type of variable, all possible choices are:
	- s for string
	- d for double
	- u for unsigned int
	- i for int
+ name is the name of the argument, some examples are:
	- wf for the type of wave function
	- n for the number of sites
	- N for the number of colors
	- m for the number of particles per sites
+ val is the value of the argument

Here is an example of a simple command:

	./mc -d -u:tmax 10 -s:wf square-mu -u:N 4 -u:m 2 -u:n 20 -i:bc 1 -d:mu 0.1 -u:obs 1,3

which runs a simulation for 10 seconds on the wave function named square-mu for
SU(4) with $m=2$ particles per site and $n=20$ sites on a periodic lattice and
an on site chemical potential of $\mu=0.1$ and mesures the bond energy and long
range-correlations (observables 1 and 3).

The code currently supports:

+ any value of $N>1$
+ any number of particle per site $m>N$
+ any number of site $n$
+ any simple boundary condition (periodic, antiperiodic, none)
+ any number of particle of each color (having more particles of one colour is allowed)
+ first neighbour Heisenberg model (easily generalisable to further neighbours)
+ inhomogeneous coupling between sites 

The values that the option obs take defines what observables are measured (the
energy per site is always measured):

0. all possible observable
1. bond energy
2. long range correlations and associated structure factor
3. colour occupation
4. variance of the energy

Beside performing a VMC simulation

+ can visualise the lattice with unit cell, site numbers, basis vectors,
  hopping amplitude, fluxes 


The executable min intends to perform a minimisation of variational parameters
and study allows an analysis of these results in a web browser.

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

All these software can be set and personalised in the lib/config.mk and file.
For instance if is not available to convert images, only these two files
require some changes.

The code is commented and the documentation can be generated via Doxyfile using

	make ref