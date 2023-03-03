# nougat: A Toolkit For Analysis of Membrane Disruption by Proteins and Other Inclusions

Proteins and other inclusions are known to interact with membrane lipids in numerous ways, including through the differential recruitment of particular lipid species, local realignment of lipid tails, and even inducing comparatively long-range deformations of the membrane surface. Coarse Grain Molecular Dynamics (CG-MD) simulations offer an attractive ‘computational microscope’ through which to observe these phenomena. While some packages do exist for the analysis of lipid binding, we are not aware of any that offer a holistic view of protein-membrane interactions. We introduce nougat, a toolkit for quantitative analysis of several measures of interest local to a membrane protein: membrane thickness, membrane curvature, lipid density, and lipid tilt.

#### Terminology note:
- z1: outer leaflet surface
- z2: inner leaflet surface
- z+: the midplane between z1 and z2
- z0: the _actual_ interface between z1 and z2 lipid tails; classically assumed to be equal to z+, but not always the case!

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
argparse
glob
os
```

### Github Dependencies:
This script utilizes Dr Jerome Henin's **vecexpr** package to drastically improve the efficiency of binning and allow for the measurement of multiple membrane quantities across a trajectory. It also makes use of his **qwrap** protocol to speed up initial setup. Please visit the repos below and follow the installation instructions.
```
github.com/jhenin/qwrap
github.com/jhenin/vecexpr
```


## What will it do?
nougat.tcl will analyze your trajectory in VMD and output .dat files that contain the average quantity per bin for each measure of interest per frame. nougat.py will read in those .dat files, extract curvature and thickness information from the heights, generate averages across the trajectory, and output its findings in .dat, .npy, and .pdf files. It will also generate a pdb file allowing the user to explore the average membrane surfaces in 3D in VMD. 

## How To Use Nougat:
### Before You Start:
Gromacs coarse-grain trajectories must be made whole using trjconv's `-pbc whole` command if you want your order parameters to make sense

### Make Edits and Check Your Work
- Open nougat_config.txt in the text editor of your choice and make edits:
```
Add the correct paths to the directories where your qwrap.so and vecexpr.so files are stored 
Edit all of the indicated fields to match your system
```
- Open VMD and load your structure and trajectory files. 
- (optional) Delete frame 0
- Open the tkConsole
- Source nougat.tcl (you may need to specify the path)
- Run the cell_prep proc:
```
cell_prep [path to nougat_config.txt] 1
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

###  Run nougat.py
- python3 [path-to]/nougat/plotting/nougat.py [path-to]/nougat/plotting/nougat_plot_config.txt [system name] [-p for polar]

### Explore your results!
All results are available as heatmaps saved as 700 DPI pdf files. To change the scale, adjust the max and min values in nougat.py (lines 68-79).

You can also load the [system name].avgheight.pdb file into VMD to explore the average surfaces in 3d. 
