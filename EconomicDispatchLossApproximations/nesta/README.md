NESTA - v0.7.0
=================
### NICTA Energy Systems Test Case Archive (NESTA)

This repository contains a collection of power system test cases with a focus on creating challenging tests for optimization algorithms on AC transmission systems.

The initial transmission network data is curated from a variety of sources including, academic publications, [IEEE test case archive](http://www.ee.washington.edu/research/pstca/), [Matpower](http://www.pserc.cornell.edu/matpower/), and the [Edinburgh test case archive](http://www.maths.ed.ac.uk/optenergy/NetworkData/introduction.html).
Each file contains detailed information as to the original source(s) of the data.
These networks are modified to produce more challenging test cases for optimization applications, such as Optimal Power Flow (OPF).
The [NESTA report](http://arxiv.org/abs/1411.0359) explains the motivations and procedures for developing these modified test cases in detail.
The revised versions of these networks are provided in the Matpower case format.
This archive is made available strictly for research purposes.

This repository is currently in an alpha release and is subject to significant changes.


## Test Case Overview
The test cases are organized into directories as follows:

* **opf** - AC Optimal Power Flow cases, as originally specified.  The power flow solution provided in the case file is the locally optimal solution produced by solving the case using [IPOPT](https://projects.coin-or.org/Ipopt)
  * **opf/api** - heavily loaded test cases (i.e. binding thermal limit constraints)
  * **opf/nco** - cases that are useful for testing nonconvex optimization methods
  * **opf/rad** - radial topology cases
  * **opf/sad** - small phase angle difference cases (i.e. binding phase angle difference constraints)
  * **opf/utl** - cases that are helpful for testing software and verified network data from publications


## Citation Guidelines
This archive is always improving and growing.  As a result, it is critically important to indicate the version number when referencing the archive.  Dataset citations should reference the [technical report](http://arxiv.org/abs/1411.0359) document.

If only a subset of the archive is used, referencing the original source document indicated in the file headers is encouraged.

----------------------
### Developed by:

NICTA, Optimisation Research Group, Energy Systems Team

June, 2017

Direct questions and comments to the NESTA Administrator:

Carleton Coffrin - cjc@lanl.gov


