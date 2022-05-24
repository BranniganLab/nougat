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
```

### Github Dependencies:
This script utilizes Dr Jerome Henin's **vecexpr** package to drastically improve the efficiency of binning and allow for the measurement of multiple membrane quantities across a trajectory. It also makes use of his **qwrap** protocol to speed up initial setup. Please visit the repos below and follow the installation instructions.
```
github.com/jhenin/qwrap
github.com/jhenin/vecexpr
```

## How To Use Nougat:
### What will it do?
nougat.tcl will analyze your trajectory in VMD and output .dat files that contain the average quantity per bin for each measure of interest per frame (in the case of height and tilt) or for the whole trajectory (in the case of density). Plotting/nougat.py will read in those .dat files, extract curvature information from the heights, generate averages across the trajectory, and output its findings in .dat and .png files. It will also generate a pdb file allowing the user to explore the average membrane surfaces in 3D in VMD. Beadcolor.tcl will import curvature and density data into the user fields of the pdb file in VMD.

### Make Edits and Check Your Work
- Open nougat.tcl in the text editor of your choice and make edits to two sections:
```
Add the correct paths to the directories where your qwrap.so and vecexpr.so files are stored (lines 6-7)
Edit all of the indicated fields in cell_prep proc (lines 12-90)
```
- Open VMD and load your structure and trajectory files. 
- Ideally, ensure that your first frame contains a flat bilayer IE the starting condition from your simulation
- Open the tkConsole
- Source nougat.tcl (you may need to specify the path)
- Run the cell_prep proc:
```
cell_prep [the name you wish to use for your files e.g. "PO" for a 100% POPC system] -1
```
- After completion, visually inspect your trajectory in VMD for the following:
```
Is the inclusion properly centered in your system?
Is the membrane relatively stable in the xy plane?
Color by user with the value 1 and select "Update Selection Every Frame"
Scroll through your trajectory: are you only seeing lipids in the outer leaflet?
Switch to user 2 and repeat; are you only seeing lipids in the inner leaflet?
```
- if you answered "no," to any questions: troubleshoot; otherwise you're ready to use nougat!

### Run nougat.tcl
If you are interested in just analyzing the principle membrane fields z1, z2, z0, and z+ using the $headname you specified:
```
nougatByField [system name] [dr] [Ntheta] [start] [end] [step]
```
If you want to analyze along every bead in the lipid chain:
```
nougatByBead [system name] [dr] [Ntheta] [start] [end] [step]
```
An explanation of the options is as follows:
```
- system name       the name you want your files/directories to be saved with
- dr                the radial bin width in Angstroms; we recommend starting with 6
- Ntheta            the number of bins in the theta direction; we recommend starting with 30
- start             the frame from which you would like to start analysis
- end               the frame where you would like to stop analysis; -1 will run for the whole simulation
- step              allows you to skip frames; otherwise use 1
```

### Edit & Run nougat.py
- Open nougat.py in the text editor of your choice
- If you used nougatByBead, set **readbeads** to 1
- If you _did not_ comment out Protein_Position, set **protein_onoff** to 1
- Edit name_list to contain the name(s) of your system(s) as specified by "system name" in nougat.tcl
- If you used nougatByBead, add your lipid(s) to the bead_dict dictionary as needed; ensure the key is equal to the name in name_list!
- run nougat.py from the command line

### Explore your results!
All results are available as heatmaps saved as 700 DPI png files. They will be saved in the folder labelled with your system name. You can also load the [system name].avgheight.pdb file into VMD to explore the average surfaces in 3d. 
