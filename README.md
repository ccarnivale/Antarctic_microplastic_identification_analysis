# Antarctic_microplastic_identification_analysis

This is a culmination of final and initial working scripts identifying, visualizing, and analyzing microplastics from the Southern Ocean within the Western Antarctic Penninsula Region, specifically waters from Paradise Bay.

There are several packages utilized to make this project work. None more important than the Hyperspec package to read in and aid in the microplastic identification automation and the OpenSpecy spectral libraries. I didn't use the OpenSpecy package as the aim was to automate the identification process and reduce the time to plastic identification with the shear number of raw spectra and samples assessed in this project.

The R markdown file contains the original implementation of the identification algorithm (Raman_pipelinem.Rmd).
The Rscript that goes with the .Rmd file is the Raman_functions.R. This contains initial attempts to create fuctions to increase reproducibility and aid in readability of the identification pipeline.

Additional Rscripts wiil A) contain the identification algorithm B) contain helper functions to be sourced into the algorithm file 
A) Full_lib_plastic_ID_algOpt.R
B) Plastic_Identification_functions.R [Currently in prep and basically empty]

#*NOTE* There was an initial version of script A included here as well but was an unoptimized version of the algorithm.
