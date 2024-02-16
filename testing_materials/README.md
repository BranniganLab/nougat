# nougat: Testing

To run tests: 
```
bash run_full_test_suite.bash
```

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
pathlib
pytest
shutil
```

### Github Dependencies:
This script utilizes Dr Jerome Henin's **vecexpr** package to drastically improve the efficiency of binning and allow for the measurement of multiple membrane quantities across a trajectory. It also makes use of his **qwrap** protocol to speed up initial setup. Please visit the repos below and follow the installation instructions.
```
github.com/jhenin/qwrap
github.com/jhenin/vecexpr
```


## What will it do?
The test suite will run unit tests on nougat.tcl and then run nougat on a few tests systems. It will then check to make sure that the results from the test systems match saved reference data.
