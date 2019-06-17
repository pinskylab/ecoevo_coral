# ecoevo_coral
Eco-evo models for coral adaptation to climate change

File Descriptions:

EcoEvoModel_MultipleScenarios_EndResults.r
-	This file runs the eco-evo model adapted from Norberg et al in parallel across many biological scenarios. It stores the model results for set time points (20, 50, 100, and 500 years).

EcoEvoModel_SingleScenario_FullTimeSeries.r
-	This file runs a single biological scenario (V and D parameter set) in the eco-evo model adapted from Norberg et al in parallel across 100 iterations and all prioritization strategies. The entire time series of species cover, trait values, temperature conditions, and managed area layout are returned.

NorbergFunctions.r
-	Functions needed to run the model adapted from Norberg et al. This file is run within each of the two previous files in order for the functions to be available in R.

Norberg_StartingParameters_iterationlist_reservesize.csv
-	Parameter values for different biological scenarios. These are read into R by row and passed to the eco-evo model.

Norberg_StartingParameters_test.csv
-	this file is used only as a template for storing parameter values. The values are filled in from the previous scenario file (Norberg_StartingParameters_iterationlist_reservesize.csv).
