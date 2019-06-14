# vbs_analysis
The repository contains the principal files used in the analysis.

<b> 4l channel, ZX data driven MC </b>

For this part the following CMSSW version was used: CMSSW_8_0_26_patch1

The folder data_driven_MC must be included in the src directory.
The following instructions are necessary:

The discriminant is defined in: src/run_makeZX.cpp

cd ext/

sh compile_ext.sh

cd ..

source set_library.sh

make run_makeZX

./run_makeZX

The necessary variables (ZZMass dbkg_kin weight) will be stored in ZX(+suffix).root

<b> 4l channel </b>

For this part the following CMSSW version was used: CMSSW_7_4_7

plotter.c -> The file reads 4l MC + data (from relevant repositories), togetherwith the data driven zx component, contained in the zx.root file. It performs relevant selection, mela cuts, kin_variable generation. The output is the 1D histogram portarying all contributions.

bkg_Workspace.c -> Generation of templates for signal and background (devided by component: vbs, qqZZ etc..). The templates are created at workspace level and saved accordingly.

The code is run with:

source runbkg.sh 

zx.c -> Reads the ZX.root file, which is the final step of the data-driven ZX monte carlo generation procedure. The output file is zx.root (taken in input by plotter.c and new_plotter.c as explained above).

<b> Combine </b>

All cards can be found in the COMBINATION_FOLDER, with sub-directories for single channels and combination of channels. I suggest a double check on the systematics and whether all channels and processes are correctly included/positioned in the right column etc. The following are the Combine commands which I used:


To combine multiple cards: 

combineCards.py Name1=old_card1.txt Name2=old_card2.txt .... > new_card.txt

To run the likelihood analysis (expected significance, without systematics): 

combine -M ProfileLikelihood --significance card_name.txt -t -1 --expectSignal=1 -S 0 --toysFreq

To run the likelihood analysis (expected significance, with systematics): 

combine -M ProfileLikelihood --significance card_name.txt -t -1 --expectSignal=1 -S 1 --toysFreq

To run the likelihood analysis (observed significance): 

combine -M ProfileLikelihood --significance card_name.txt 


<b> Note </b>

These are the single files, and the most complete version I was able to find. However I've included .tar compressed file of the complete directories here: https://www.dropbox.com/sh/mahvupqcauchx53/AAAhKMKCSfl-Bhsw1uwJo5HQa?dl=0 (quite messy, but containing a long list of different file versions, data cards, plots and so on). It may be useful in case there is some file not included in the above list, but which is used somewhere in the code.

The Meng.txt contains the instructions that Meng gave me at our CERN meeting last Spring. 
