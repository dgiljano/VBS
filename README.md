# vbs_analysis
The repository contains the principal files used in the analysis.

<b> 4l channel, ZX data driven MC </b>

For this part the following CMSSW version was used: CMSSW_8_0_26_patch1

The folder must be included in the src directory.
The following instructions are necessary.

The discriminant is defined in: src/run_makeZX.cpp

cd ext/

sh compile_ext.sh

cd ..

source set_library.sh

make run_makeZX

./run_makeZX

The necessary variables (ZZMass dbkg_kin weight) will be stored in ZX.root

<b> 4l channel </b>

For this part the following CMSSW version was used: CMSSW_7_4_7

plotter.c -> The file reads 4l MC + data (from relevant repositories), togetherwith the data drive zx component, contained in the zx.root file. It performs relevant selection, mela cuts, kin_variable generation. The output is the 1D histogram portarying all contributions.

new_plotter.c -> Same functions as the previous file, with the exclusion of 1D histogram generation. In addition the file enables the user to generate the mass+kin_variable 2D histograms and to save them on a dedicated root file.

bkg_Workspace.c -> Generation of templates for signal and background (devided by component: vbs, qqZZ etc..). The templates are created at workspace level and saved accordingly.

zx.c -> Reads the ZX.root file, which is the final step of the data-driven ZX monte carlo generation procedure. The output file is zx.root (taken in input by plotter.c and new_plotter.c as explained above).

<b> 2l2q channel </b>

For this part the following CMSSW version was used: CMSSW_8_0_26_patch1

plotMC.c -> This file is based on Roberto's original 1D template generator for the 2l2q channel. It performs the 2D mass +kin_variable histogram generation for the different components. This version is not automated, so for each contribution (for example, ttbar merged and so on) it is necessary to change the output file names). 

bkgWorkspace_2l2q.c -> Generation of templates for signal and background. The templates are created at workspace level and saved accordingly. The file creates the templates for both the merged and the resolved components. However, since COMBINE requires the two contributions to be saved on two different files, it is necessarry to comment the lines (at the end of the file) where the unwanted templates are imported on the workspace. As can be seen, in this file version this is done for the merged component.

alpha.C -> This file has the same structure of the plotMC.c file. It embeds the alpha method used for the the DY component treatment, by operations with 1D histograms. The default output is the whole set of possible 1D plots (all possible contributions and variables). 

<b> COMBINE </b>



<b> note </b>

These are the single files, and the most complete version I was able to find. However I've included a .tar compressed file of the complete directories (quite messy, but containing a long list of different file versions, data cards, plots and so on). It may be useful in case there is some file not included in the above list, but which is used somewhere in the code.
