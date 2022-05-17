# nougat

This script utilizes Dr Henin's vecexpr package to drastically improve the polar height analysis script and allow for the measurement of membrane curvature across a trajectory.

To use: edit the paths at the top of polarH3.tcl, source polarH3.tcl in VMD. If interested in primary membrane fields (z1, z2, z0, z+) use run_field_mult proc. If interested in fields z1, z2, and z+ as defined by different beads along lipid tail, use run_bead_mult proc. 

Plotting/polarplotH_vec.py does all curvature analysis, plotting, and constructs a PDB file to allow further exploration of results in VMD.
