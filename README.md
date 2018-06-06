# BoostedRazorAnalysis                                                                                                       
Run II Search for Supersymmetry using Razor variables in the boosted regime

## Recipe

```Shell
git clone https://github.com/cghuh/DrawCodesforRazorBoost
```

## Code description

### BackgroundEst
#### MRR2.cc:
	* Make unrolled MR-R2 plot
	* Caluclate the background estimation with bin-by-bin ratio using unrolled MR-R2  

#### MRR2_universialRatio.cc:
	* Caluclate the background estimation with universial ratio using unrolled MR-R2

#### BGest.cc:
	* Caluclate the background estimation using each MR, R2 

### StackedPlot
#### StackedPlot.C:
	* Draw all of 1-D histograms in the root files

	* Rebin.cfg		: rebin the variables what they contain the specific character
	* RootFilePath.cfg	: Specify the root file path
	* Settings.cfg		: Root files name for drawing (You have to fix the color option in StackedPlot.C dirctly

#### StackedPlot_2D.C:
	* Draw all of 2-D histograms in the root files
	* Use same configure files with StackedPlot.C

### TriggerEff:
	* Calculate the trigger efficiency

### WTagScaleFactor:

	* Calculate the W/Top tagging(fake-rate) scale factor in data/simulation and FullSim/FastSim
	
	* 1D/WTagScaleFactor.cpp	: Calculate the W/Top fake rate and that's scale factor of data/simulation in pT
	* 1D/FullFastSimScaleFactor.cpp	: Calculate the W/Top fake rate and that's scale factor of FullSim/FastSim in pT
	* 2D/WTagScaleFactor.cpp	: Calculate the W/Top fake rate and that's scale factor of data/simulation in pT, eta
	* 2D/FullFastSimScaleFactor.cpp	: Calculate the W/Top fake rate and that's scale factor of FullSim/FastSim in pT, eta

### ZBackgroundEstimation:

	* Calculate the prompt photon purity using template fitting for z->vv background estimation
	
	* ZpurityFit.cpp : fit the charged isolation in each MR, R^2, # of mass tagged W and ECAL barrel/endcap. Prompt template is gamma + jet sample (medium id), fake template is QCD sample (medium id, but invert sigma_ietaieta)

### cutflow:

	* cutflow*.py		: Make the tex file of cutflow table automatically. Set the background samples and selection.
	* get*BGcomposition.py	: Make the tex file of background composition in each region automatically.
	* list*			: List of samples of background and data for W and Top analysis.
	1. If you get the code first time, you have to fix the data location(inputdir) of signal in python file.
	2. When you use use this, type as "python cutflow_W.py list". You should python *W* list_W for W ana, python *Top* list_Top for Top ana.
	3. After running the code, you should fix the cut name as our AN for cutflow table and put the name of signal and background region for the background estmation table also as AN.

### shapecompare :
	* Compare the MR and R^2 shape in signal and background control in each background sample. For example, compare Q control region and signal region with Multijet sample. (TT-TvsS, WJet-WvsS, ZJet+DY-ZvsS)

## Contact Information

Changgi Huh (changgi.huh@cern.ch)
