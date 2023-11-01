#Load all packages necessary for hyperspec analysis--------

#Running script on Temple High-performance Computing Server (HPC)
#https://www.hpc.temple.edu/compute/
#Linux based server and run from command line

#Need to install packages first though
#Need to direct with repos to download the install packages since we are running 
#on HPC servers

#install.packages("png", repos = "http://cran.us.r-project.org")
#install.packages("bmp", repos = "http://cran.us.r-project.org")
#install.packages("pixmap", repos = "http://cran.us.r-project.org")

#Packages I acutally use
#install.packages("baseline", repos = "http://cran.us.r-project.org")
#install.packages("hyperSpec", repos = "http://cran.us.r-project.org")
#install.packages("prospectr", repos = "http://cran.us.r-project.org")
#install.packages("ggplot2", repos = "http://cran.us.r-project.org")
#install.packages("magrittr", repos = "http://cran.us.r-project.org")
#install.packages("ggpubr", repos = "http://cran.us.r-project.org")
#install.packages("signal", repos = "http://cran.us.r-project.org")

#install.packages("OpenSpecy")
#install.packages("tidyverse")

#library(tidyverse)
#library(png)
#library(bmp)
#library(pixmap)
library(baseline)
library(hyperSpec)
library(prospectr)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(signal)
#library(OpenSpecy)

#Unfortunately Tidyverse can't be loaded onto HPC server so all indexing and 
#data manipulation will have to be done in base R.

#All functions for the project are based on OpenSpecy but optimized and automated
# to run Thousands of spectra instead of indivual spectra

#This version of the script will be used with the original version of the 
#OpenSpecy database I downloaded. The oringinal had ~6,000, while the new version
#has 40000 reference spectra. I will likely filter out non-plastic spectra
#to reduce the amount of correlations needed to be calculated.

# Main purpose: This script takes folders (1 depth/size/station sampling) of 
#indivial spectra (1 per .txt file) and converts into objects for each sampling 
# the counts of the: total plastic count and count for each individual polymer

#Next script will convert each sampling's count of polymers into a single df
# that will perform feature engineering extraction based on the microplastic 
# composition of each sample.

#Basics outline of the code

# 1) Environment setup with installation and calling of the necessary dependent
# packages

# 2) Read in of both the reference data and meta data
# This needs to be cleaned and merged with original dataset and will likely
# filter out non-plastic spectra as they will be removed if it doesn't make the
# Hit-index threshold regardless
# Data will also be tidyed and smoothed along the same wavelength range to keep
# the data consistent for identification

# **ALL DATA are cleaned**, and smoothed along the same range and number of steps
# cor() analysis needs the same number of observations from each distribution
# to properly work.

# 3) Read in files from folders into individual objects with functions sourced
# from another file
#Will start with the raw for loops in the original scripts for trouble shooting
#and switched to the sourced functions I make after the first one works
# This function will read in the raw .txt, perform the data tidying, and smoothing
# to be in congruence with the reference df all in one step

# 4) Microplastic Identification of the raw spectra and counts of each will be 
# converted

#The highest correlation value will be assigned to each individual spectra
# and will filter out all spectra that don't reach the threshold correlation value
# the remaining will be considered plastics and counted resulting in a final
# summary table

# 5) export tables as csvs or combine first with meta data and then export.
# currently undecided here.

#Setting and checking working directory during hpc computing
setwd("/home/tud11809/Antarctic_microplastics")

print(getwd())

#Read in files for spectral library-------
#OpenSpec library
#Must have these files in working directly to call them like I have here
open_spec_reflib <- read.csv("raman_library.csv")
#open_spec_reflib <- read.csv("/Users/christophercarnivale/Desktop/Dissertation_data/Microplastic_antarctica/raman_library.csv")

open_spec_reflib_metadata <- read.csv("raman_metadata.csv")
#open_spec_reflib_metadata <- read.csv("/Users/christophercarnivale/Desktop/Dissertation_data/Microplastic_antarctica/raman_metadata.csv")

#Original Plastic Library
ref_lib <- read.csv("spectra_db.csv")
#ref_lib <- read.csv("/Users/christophercarnivale/Desktop/Dissertation_data/spectra_db.csv")

#str(ref_lib)

plastic_ref_lib <- ref_lib %>% 
  dplyr::filter(SOURCE == "plastic")
#Need to merge libraries and add information to the metadata csv file for my use

#To merge I need to match the 2 formats of the datasets
plastic_ref_lib_for_merge <- plastic_ref_lib[,c("WAVE", "INTENSITY", "poly_lab")]

#Add an column with a fake name to overwrite
plastic_ref_lib_for_merge['sample_name'] <- 1

#Polymer labels to be interated over
poly_labs <-as.character(unique(plastic_ref_lib$poly_lab))

#Filling the column with the new identifyer numbers to match OpenSpecy ref df
#I need to download the updated database to make sure these numbers arent taken
#and for loop still works.
for(i in seq_along(poly_labs)){
  plastic_ref_lib_for_merge[which(plastic_ref_lib_for_merge$poly_lab == 
                                    poly_labs[i]),]['sample_name'] <- rep(seq(623, 623+length(poly_labs), by = 1)[i], 
                                                                          nrow(plastic_ref_lib_for_merge[which(plastic_ref_lib_for_merge$poly_lab 
                                                                                                               == poly_labs[i]),]))
} 

#Checking to see If the output is in the proper format for databased merger
#print(head(plastic_ref_lib_for_merge))
#print(head(open_spec_reflib))

#Order and names were not exactly the same...change to match for merger
colnames(plastic_ref_lib_for_merge) <- c("wavenumber","intensity", "poly_lab","sample_name")

#Checking to make sure it is now compatible
print(head(plastic_ref_lib_for_merge[,c("wavenumber","intensity", "sample_name")]))
print(head(open_spec_reflib[,c("wavenumber","intensity", "sample_name")]))

#Changing the shape of the df to merge
plastic_ref_lib_for_merge <- plastic_ref_lib_for_merge[,c("wavenumber","intensity", 
                                                          "sample_name")]
open_spec_reflib <- open_spec_reflib[,c("wavenumber","intensity", "sample_name")]

#Merge the dataframes to make a complete reference library
merged_plastic_ref_lib <- rbind(plastic_ref_lib_for_merge,open_spec_reflib)

#Commenting out all of the debug checks as we move forward
#print(head(merged_plastic_ref_lib))

#print(unique(merged_plastic_ref_lib$sample_name))

#To loop over files I need to setup the list of filenames to properly loop over
plastic_test_files <- list.files(path = "/home/tud11809/Antarctic_microplastics/7_06_23/Stn_R/Surface/less_20/", pattern = "*.txt")
#SSH filepath                            
#"/home/tud11809/Antarctic_Microplastic/7_06_23/Stn_R/Surface/Less_20/"
#local filepath
#"/Users/christophercarnivale/Desktop/Dissertation_data/Raman Spectral Data copy/Antarctic_Microplastic_Raman_Spectra/7_06_23/Stn_R/Surface/less_20/"
file_numbers <- seq(plastic_test_files)

#Blank matrix table to be inserted into an hyperspec object
#plastic_test_spc_table <- matrix()
#Ended up not using this 

#Sequence range for which all spectra will be processed by
#I could make this a function call option in the future functions but I will
#hard code it in for now...
freq_index <- seq(780, 1750, by = 0.1)

#Empty dataframe from which to add in each spectra
all_sample_particle_WV <- as.data.frame(freq_index)

print(plastic_test_files)

#for loop over .txt files and create a df of spectra
#This loops and adds the smoothed version of the spectra within the wavelength 
#ranges desired, which can be changed at any time by changing the "freq_index"
#object
for(filename in plastic_test_files){
  
  temp_csv <- read.delim(file = paste0("~/Antarctic_microplastics/7_06_23/Stn_R/Surface/less_20/",filename), header = F, sep = "\t")
  colnames(temp_csv) <- c("Wavelength", "Counts")
  temp_csv_matrix <- as.matrix(temp_csv)
  approx_test_fun <- approxfun(x = temp_csv_matrix[,'Wavelength'],
                               y = temp_csv_matrix[,'Counts'], rule = 2)
  new_counts <- approx_test_fun(freq_index)
  new_counts_smoother <- sgolayfilt(new_counts, 3,11, m = 1)
  all_sample_particle_WV[filename] <- new_counts_smoother
}

#All print statements prior to end count is debug catch steps

#print(str(all_sample_particle_WV))
#print(head(all_sample_particle_WV))

#Matrix is in the wrong format and needs to by transposed to fit the hyperspec obj
#Does all functions on a per row basis instead of a percolumn basis, which is 
#counterintuitive

#print(head(t(all_sample_particle_WV))) 
#To view the structure of the transposed matrix before commiting to it

all_sample_particle_WV_t <- t(all_sample_particle_WV)

#Converting the first ROW into column names and then removing the row from the df
colnames(all_sample_particle_WV_t) <- all_sample_particle_WV_t[1,]

all_sample_particle_WV_t <- all_sample_particle_WV_t[-1,]

#Extracting the filenames or name of specific spectra being analyzed
#hyperspec (HS) class is an S4 object with very strict parameters for creating a 
#hyperspec obj and the spectra metadata, in this case filenames (particle #),
#are saved in a separate vector and called to be added to the HS obj.
potplastic_filename <- rownames(all_sample_particle_WV_t)

#Making hyperspec obj with matrix of spectra with filename meta data
potp_hyper <- as.hyperSpec(all_sample_particle_WV_t, data = data.frame(potplastic_filename))

#Baseline correction
potp_baselines <-  spc.fit.poly.below(potp_hyper, poly.order = 4)
#Plotting baseline corrections
#plot(potp_hyper-potp_baselines)
#title(xlab = "Wavelength",ylab = "Arbitrary unit of excitation")

potp_hyper_corrected <- potp_hyper-potp_baselines

#Standard Normal variate normalization
potp_hyper_corrected@data[["spc"]] <- standardNormalVariate(potp_hyper_corrected@data[["spc"]])

#min-max scaling
potp_hyper_corrected@data[["spc"]] <- (potp_hyper_corrected@data[["spc"]] - min(potp_hyper_corrected@data[["spc"]])) / (max(potp_hyper_corrected@data[["spc"]])-min(potp_hyper_corrected@data[["spc"]]))
#}

#plot(potp_hyper_corrected, col = rainbow(39),title.args = list(xlab = "Wavelength",ylab = "Arbitrary unit of Intensity"))+legend(500, 1, legend = potplastic_filename,col = rainbow(39))
#pdf("test_hyperspec_plot.pdf")
#plot(potp_hyper_corrected,col = rainbow(40),title.args = list(xlab = "Wavelength",ylab = "Arbitrary unit of Intensity"))
#dev.off()

#Extracting unique names from reference library
#poly_labs <- as.character(unique(ref_lib$poly_lab))
#remove the last column (why?)
#poly_labs <- poly_labs[-15]

#This previous version took all samples from the study and not just the reference platic. 
#All of the categorized "." were environmental samples thus necessitating removal. 
#I had a filtered df for plastic reference only so use that df for plastic label names.

#Creating vector of 
reference_labs <-as.character(unique(merged_plastic_ref_lib$sample_name))
#Create an correlation matrix to fill that are the dimensions of the plastic reference vs the number of spectra collected within a single sample
cor_table_manyplastic <- as.data.frame(matrix(ncol = length(reference_labs), nrow = nrow(potp_hyper_corrected@data[["spc"]]), dimnames = list(NULL, reference_labs)))
#Calculating correlations and filling the DF with the values

#To optimize I am going to move the library smoothing to another df as a saved obj
#and that obj to be called for the correlation analysis

merged_plastic_ref_lib_smoothed <- as.data.frame(freq_index)

for(i in seq_along(merged_plastic_ref_lib$sample_name)){
temp_filter <- dplyr::filter(merged_plastic_ref_lib, sample_name == merged_plastic_ref_lib$sample_name[i])
approx_test_fun <- approxfun(x = temp_filter[,'wavenumber'],
                             y = temp_filter[,'intensity'], rule = 2)
temp_new_counts <- approx_test_fun(freq_index)
temp_new_counts_smoother <- sgolayfilt(temp_new_counts, 3,11, m = 1)
merged_plastic_ref_lib_smoothed[i] <- temp_new_counts_smoother
}

print(merged_plastic_ref_lib_smoothed[1:5,1:5])

#for(row in 1:nrow(potp_hyper_corrected@data[["spc"]])){
#  for(i in seq_along(merged_plastic_ref_lib$sample_name)){
#    cor_table_manyplastic[row,i] <- round(cor(potp_hyper_corrected@data[["spc"]][row,], temp_new_counts_smoother, method = "pearson"), digits = 3)
#  }
#}
#Creating an empty df to hold the identified plastic info in
#plastic_lab <- as.data.frame(matrix(nrow = nrow(cor_table_manyplastic), ncol = 2))

#colnames(plastic_lab)  <- c("Hit_Index_Value","Polymer_ID")

#rownames(plastic_lab) <- plastic_test_files
#Search the correlation matrix for plastic that are the highest correlated to the spectra with the filename (or really particle name).
#for (d in 1:nrow(cor_table_manyplastic)){
#  plastic_lab[d,1]<-cor_table_manyplastic[d,which.max(cor_table_manyplastic[d,])]
#  plastic_lab[d,2]<-colnames(cor_table_manyplastic[which.max(cor_table_manyplastic[d,])])
  
#}
#Filter spectra that don't meet the spectra limit
#plastic_lab_filt <- plastic_lab %>% 
#  rownames_to_column('filename') %>% 
#  dplyr::filter(Hit_Index_Value > 0.22) %>% 
#  column_to_rownames('filename')

#Summarise the counts for each polymer - But I need to have a way to add to a df with polymers unidentified in the sample but we have references for.
#plastic_lab_count <- group_by(plastic_lab_filt, Polymer_ID) %>% summarise(n())

#plastic_lab_count
#plastic_lab_filt

#ggtexttable(plastic_lab_count)

#ggtexttable(plastic_lab_filt)