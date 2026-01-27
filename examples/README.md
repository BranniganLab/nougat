This directory contains the outputs from nougat.tcl after it was run on a 
test trajectory of a protein embedded in a membrane. The example.py file
is intended to demonstrate some of the ways in which nougat can be used.

To run the example file, please ensure that nougat is pip installed and
then enter the following command into your terminal:

> python3 ./example.py . -p

This will produce the following files: 

1. a 2D heatmap of the outer leaflet average height (example_image.pdf)
2. a .pdb file containing the average height surfaces (membrane_heights.pdb)
3. text files containing triangle coordinates for each of the surfaces ([surface]_avg_surface.txt)

To draw the triangulated height surfaces in VMD, open VMD and load the
membrane_heights.pdb file. Then, open the tkConsole and enter:

> source [path to]/nougat/tcl/utilities/draw_triangles.tcl  
> drawTriangles [path to surface .txt file] top

For more details, see the docstring in draw_triangles.tcl.
