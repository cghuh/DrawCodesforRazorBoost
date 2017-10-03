# BoostedRazorAnalysis                                                                                                       
Run II Search for Supersymmetry using Razor variables in the boosted regime

## Recipe

```Shell
git clone https://github.com/cghuh/DrawCodesforRazorBoost
```

## Code description

BackgroundEst/MRR2.cc:
	* Make unrolled MR-R2 plot
	* Caluclate the background estimation with bin-by-bin ratio using unrolled MR-R2  

BackgroundEst/MRR2_universialRatio.cc:
	* Caluclate the background estimation with universial ratio using unrolled MR-R2

BackgroundEst/BGest.cc:
	* Caluclate the background estimation using each MR, R2 

StackedPlot/StackedPlot.C:
	* Draw all of 1-D histograms in the root files

	* Rebin.cfg		: rebin the variables what they contain the specific character
	* RootFilePath.cfg	: Specify the root file path
	* Settings.cfg		: Root files name for drawing (You have to fix the color option in StackedPlot.C dirctly

StackedPlot/StackedPlot_2D.C:
	* Draw all of 2-D histograms in the root files
	* Use same configure files with StackedPlot.C

TriggerEff:
	* Calculate the trigger efficiency with 

WTagScaleFactor:

	* Calculate the W/Top tagging(fake-rate) scale factor in data/simulation and FullSim/FastSim
	* 1D/WTagScaleFactor.cpp	: Calculate the W/Top fake rate and that's scale factor of data/simulation in pT
	* 1D/FullFastSimScaleFactor.cpp	: Calculate the W/Top fake rate and that's scale factor of FullSim/FastSim in pT
	* 2D/WTagScaleFactor.cpp	: Calculate the W/Top fake rate and that's scale factor of data/simulation in pT, eta
	* 2D/FullFastSimScaleFactor.cpp	: Calculate the W/Top fake rate and that's scale factor of FullSim/FastSim in pT, eta

cutflow:

	* cutflow*.py		: Make the tex file of cutflow table automatically. Set the background samples and selection.
	* get*BGcomposition.py	: Make the tex file of background composition in each region automatically.

shapecompare :
	* Compare the MR and R^2 shape in signal and background control in each background sample. For example, compare Q control region and signal region with Multijet sample. (TT-TvsS, WJet-WvsS, ZJet+DY-ZvsS)

## Contact Information

Changgi Huh (changgi.huh@cern.ch)
