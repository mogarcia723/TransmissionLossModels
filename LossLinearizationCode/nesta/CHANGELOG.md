NESTA Change Log 
=================

### v0.7.0
- Renaming util cases to utl and adding case name post-fix
- Adding all cases from ["AC Power Flow Data in MATPOWER and QCQP Format: iTesla, RTE Snapshots, and PEGASE"](https://arxiv.org/abs/1603.01533) 
- Adding selected cases from ["Local Solutions of the Optimal Power Flow Problem"](http://ieeexplore.ieee.org/abstract/document/6581918/)
- Adding the WECC 240 case from ["Reduced Network Modeling of WECC as a Market Design Protype"](https://pserc.wisc.edu/research/public_reports.aspx)
- Increasing the phase angle difference constraints slightly in the Small Angle Difference (SAD) cases, to improve numerical stability


### v0.6.1
- Updated NESTA build scripts to use PowerModels.jl
- Removed case6_ww__api case due to numerical issues


### v0.6.0
- Adding a new category of AC-OPF test cases, the NonConvex Optimization (NCO) cases
- Adding selected cases from ["Bound Tightening for the Alternating Current Optimal Power Flow Problem"](http://ieeexplore.ieee.org/xpl/login.jsp?arnumber=7328765)
- Corrected instance category descriptions in file comments
- Corrected Rate B and Rate C values in the 1979 IEEE Reliability Test System (RTS) to match the original publication
- Added nesta\_case3\_cc to utility instances


### v0.5.0
- Adding a new category of AC-OPF test cases, the Radial Topology (RAD) cases
- Adding selected radial cases from ["Inexactness of SDP Relaxation and Valid Inequalities for Optimal Power Flow"](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=7056568)
- Adding documentation of the phase angle bound used in the thermal limit upper bound model


### v0.4.0
- Adding PEGASE test cases from matpower 5.1
- Changed thermal limit models to round up to the nearest integer, rather than truncate


### v0.3.0
- Adding a new class of AC-OPF test cases, the Small Angle Difference (SAD) cases
- Removed tap ratios of 1.0 from lines in the IEEE RTS96 network


### v0.2.0
- AC power flow set-points updated to best-known AC-OPF solution


### v0.1.5
- Minor corrections on several test case comments
- Added nesta\_case3\_ch to utility instances


### v0.1.4
- Minor updates the line capacities
- Added draft translation of EIRGrid test cases
- Added "util" directory containing useful special cases for robust testing


### v0.1.3
- Adding a new class of AC-OPF test cases, the Active Power Increase (API) cases.
- Minor updates line capacities
- Author-based test case naming convention implemented


### v0.1.0
- Initial check-in
