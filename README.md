# Multi-gesture drag-and-drop decoding in a 2D iBCI control task
<img src="/readme_fig.png" width=100% height=100%>

This repo contains code for reproduction of the results and figures presented in the paper ["Multi-gesture drag-and-drop decoding in a 2D iBCI control task"](https://iopscience.iop.org/article/10.1088/1741-2552/adb180).

You'll need to download the dataset to run the code (Dryad):
https://datadryad.org/dataset/doi:10.5061/dryad.98sf7m0v1

Navigate over to `/scripts` to run `DragAndDropMaster.m`.

Upon running, figures (.png and .svg) and statistical outputs (.txt) will populate the `/outputs` folder (Figure 2A is included as an example).

The output directory can be changed to a local directory by updating the variable `saveFiguresFolder` in `DragAndDropMaster.m`. You will also want to update the `dataFolder` variable in `DragAndDropMaster.m` to specify the location of the saved data (from Dryad).

# Requirements
This code was written with **MATLAB R2021b** (v9.11) and relies on the following MATLAB toolboxes: 
- System Indentification Toolbox (v.9.15)
- Statistics and Machine Learning Toolbox (v12.2)
- BioInformatics Toolbox (v4.15.2)
