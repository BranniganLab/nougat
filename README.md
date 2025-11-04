# nougat: A Toolkit For Analysis of Membrane Disruption by Proteins and Other Inclusions

Proteins and other inclusions are known to interact with membrane lipids in numerous ways, including through the differential recruitment of particular lipid species, local realignment of lipid tails, and even inducing comparatively long-range deformations of the membrane surface. Coarse Grain Molecular Dynamics (CG-MD) simulations offer an attractive ‘computational microscope’ through which to observe these phenomena. We introduce nougat, a toolkit for quantitative analysis of several measures of interest local to a membrane protein.


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
pathlib
```

### Optional Dependencies:
This script can utilize Dr Jerome Henin's **vecexpr** package to improve the efficiency of binning. It can also make use of his **qwrap** protocol to speed up initial setup. Please visit the repos below and follow the installation instructions if you wish to enable these features. Alternatively there are pre-compiled binaries included in the _tcl/utilities_ directory that you can try, but they may not work due to varying system architectures.
```
github.com/jhenin/qwrap
github.com/jhenin/vecexpr
```

## What will it do?
nougat.tcl will analyze your trajectory in VMD and output .dat files that contain the average height per bin per frame. nougat.py will read in those .dat files to extract the following features: 
- leaflet/bilayer thickness
- mean curvature
- gaussian curvature 

nougat.py can then be used to plot, combine, and/or compare those features. 

## How To Use Nougat:

### Make Edits and Check Your Work
- Open examples/nougat_config.txt in the text editor of your choice and make edits:
```
Add the correct paths
Edit all of the indicated fields to match your system
```
- Open VMD and load your structure and trajectory files. 
- (optional) Delete frame 0
- Open the tkConsole
- Source tcl/nougat.tcl (you may need to specify the path)
- Run the cell_prep proc:
```
cell_prep [path to nougat_config.txt] 1     # the 1 runs a leaflet sorting algorithm on all frames
```
- After completion, visually inspect your trajectory in VMD for the following:
```
Is the inclusion (if present) properly centered in your system?
Is the membrane relatively stable in the xy plane?
Color by user with the value 1 and select "Update Selection Every Frame"
Scroll through your trajectory: are you only seeing lipids in the outer leaflet?
Switch to user 2 and repeat; are you only seeing lipids in the inner leaflet?
```
- if you answered "no," to any questions: troubleshoot; otherwise you're ready to use nougat!

### Run nougat.tcl

```
start_nougat [system name] [path to nougat_config.txt] [dr|Nbins] [Ntheta] [start] [end] [step] [polar]
```
An explanation of the options is as follows:
```
- system name       the name you want your files/directories to be saved with
- dr | Nbins        if using cartesian coords, the number of bins in the x/y directions | if polar coords, the radial bin width in Angstroms
- Ntheta            will be ignored if cartesian; the number of bins in the theta direction 
- start             the frame from which you would like to start analysis
- end               the frame where you would like to stop analysis; -1 will run for the whole simulation
- step              allows you to skip frames; otherwise use 1
- polar             0 for cartesian system, 1 for polar system
```

#### Terminology note:
- z1: outer leaflet surface
- z2: inner leaflet surface
- z+: the midplane between z1 and z2
- z0: the _actual_ interface between z1 and z2 lipid tails; classically assumed to be equal to z+, but not always the case!

###  Run nougat.py
An example script _examples/example.py_ has been provided to demonstrate how nougat.py can be used. You can also view the class diagram in the _python_ folder to learn more.

### Explore your results!
