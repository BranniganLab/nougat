# nougat: A Toolkit For Analysis of Membrane Disruption by Proteins and Other Inclusions

Ion channels are known to interact with membrane lipids in numerous ways, including through the differential recruitment of particular lipid species, local realignment of lipid tails, and even inducing comparatively long-range deformations of the membrane surface. Coarse Grain Molecular Dynamics (CG-MD) simulations offer an attractive ‘computational microscope’ through which to observe these phenomena. While some packages do exist for the analysis of lipid binding, we are not aware of any that offer a holistic view of protein-membrane interactions. We introduce nougat, a toolkit for quantitative analysis of several measures of interest local to a membrane protein: membrane thickness, membrane curvature, lipid density, and lipid tilt. 

## Prerequisites:

### Software Dependencies:
```
VMD
Python3
```

### Python Dependencies:
```
numpy
matplotlib
warnings
```

### Github Dependencies:
This script utilizes Dr Jerome Henin's **vecexpr** package to drastically improve the efficiency of binning and allow for the measurement of multiple membrane quantities across a trajectory. It also makes use of his **qwrap** protocol to speed up initial setup. We have provided pre-packaged versions of these utilities, but recommend that you pay the author a visit.
```
github.com/JHenin/qwrap
github.com/JHenin/vecexpr
```

## How To Use Nougat:


To use: edit the paths at the top of polarH3.tcl, source polarH3.tcl in VMD. If interested in primary membrane fields (z1, z2, z0, z+) use run_field_mult proc. If interested in fields z1, z2, and z+ as defined by different beads along lipid tail, use run_bead_mult proc. 

Plotting/polarplotH_vec.py does all curvature analysis, plotting, and constructs a PDB file to allow further exploration of results in VMD.
