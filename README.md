#  Representational drift in the mouse visual cortex
This repository includes analyses and processed data presented in Deitch et al., 2021.

## Usage and documentation
Scripts are provided in the *scripts* directory. Data sets are provided in the *data* directory.</br>
Note that the provided data sets are a processed version of the publicaly available neuronal data published in de Vries et al., 2020 and Siegle et al., 2021.

To perform the analysis on the provided data sets, use the *visual_drift_analysis.m* script.</br>
Before running the script, change the pathway in *scripts_path* and *data_path* variables to the pathway of the corresponding *scripts* and *data* directories on your computer.

## Compatibility
Was tested on Windows 10 using Matlab 2016b, 2017b and 2021a.</br>
Running the scripts on newer versions of Matlab may require downloading additional packages.

## References
1. Deitch et al., Representational drift in the mouse visual cortex, Current Biology (2021), https://doi.org/10.1016/j.cub.2021.07.062
2. de Vries et al., A large-scale standardized physiological survey reveals functional organization of the mouse visual cortex, Nature neuroscience (2020), https://doi.org/10.1038/s41593-019-0550-9
3. Siegle et al., Survey of spiking in the mouse visual system reveals functional hierarchy, Nature (2021), https://doi.org/10.1038/s41586-020-03171-x
