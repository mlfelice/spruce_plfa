# Date created: 5/1/2020
# Created by: Mark Felice
# Last Modified: 5/7/2020
# Last Modified by: Mark Felice

# This script is for importing Cameron's 2014 PLFA data and combining with 2016
# data

# Try to maintain consistent variable naming convention in this script and 
# across all of our PLFA analyses. This will make it easier to merge scripts w/
# different data sets and to track script more easily.
#
# raw named peak lists: pks_<YY...>
#   Mulitple years get listed consecutively
#   eg. for 2014/2015 data: pks_1415
# for a filepath, suffix with fp
# for the qc-check output, suffix with qc
# 
# File with d13C values calculated: d13c_<YY...>
#   Mulitple years get listed consecutively
#   eg. for 2014/2015 data: d13c_1415
#
# Dataframe of methanol d13C signatures: meth_<YY...>
#
#wts_1415
#md_1415
#ind13_1415


# Maybe we could cut down on the length of the script if we merged everything
# early on and then ran calcs?


# plfayer package used for analysis is in active development and may need to be
# reinstalled to incorporate changes. Because it is present locally, it 
# requires slightly different installation procedure. 
# To install latest version of plfayer:
devtools::install("C:/Users/Mark/Desktop/plfayer", upgrade_dependencies = FALSE)

#------------------------------------------------------------------------------


# 1. Load packages
#------------------------------------------------------------------------------

library(plfayer)
library(tidyverse)
library(readxl)

# Functions, which should be moved to source or package
iqr <- function(df, cols){
  for(c in cols){
    q_75 <- quantile(df[[c]], 0.75, na.rm = TRUE)
    q_25 <- quantile(df[[c]], 0.25, na.rm = TRUE)
    iqr <- q_75 - q_25
    
    df[which(df[[c]] < q_25 - 1.5 * iqr | df[[c]] > q_75 + 1.5 * iqr), c] <- NA  
  }
  return(df)
}

sdr <- function(df, cols){
  for(c in cols){
    min_out <- mean(df[[c]], na.rm = TRUE) - sd(df[[c]], na.rm = TRUE) * 3
    max_out <- mean(df[[c]], na.rm = TRUE) + sd(df[[c]], na.rm = TRUE) * 3
    
    df[which(df[[c]] < min_out | df[[c]] > max_out), c] <- NA
  }
  return(df)
}

filter_outliers <- function(df, cols, method = 'iqr'){
  switch(method,
         iqr = iqr(df, cols),
         sd = sdr(df, cols))
}


# 2. Import data
#------------------------------------------------------------------------------

##########################
##### 2014/2015 data #####
##########################

# Jess did statistical analyses of Cameron's data, so if we go beyond 
# exploratory data analysisneed to verify our processing pipepline and all 
# choices made regarding lipids to include peak height cutoffs, outlier 
# removal, etc. 

# plfayer is designed to work with peak lists in specific format, so will use
# modified import steps for Cameron's data, manually instead of 
# import_batch_base()

# Specify filepath for Cameron's PLFA data
pks_1415_fp <- paste0('C:/Users/Mark/Dropbox/umn_gutknecht_postdoc/', 
                       'spruce_project/plfa_13C/data/SPRUCE DATA Nov ', 
                       '2016_all the final calculations.xls')

# Sheet has all (and only) named peaks w/height, area, d13C
pks_1415 <- read_excel(pks_1415_fp, sheet = 'all samples compiled_named only')

# Rename columns to be consistent with data imported by import_batch_base()
names(pks_1415) <- c('BatchDataFileName', 'DataFileName', 'RetTimeSecs',
                      'MajorHeightnA', 'DisplayDelta1', 'Name', 
                      'TotalPeakArea1') # For Cameron's data, the correction for standard is already done, just needs methylation correction
                      
# Create Batch column and reorder columns to match plfayer format
pks_1415 <- pks_1415 %>% 
             mutate(Batch = str_extract(BatchDataFileName, 'B[0-9]+(?=_)')) %>%
             select('BatchDataFileName', 'Batch', 'DataFileName', 'RetTimeSecs', 
                    'MajorHeightnA', 'TotalPeakArea1', 'DisplayDelta1', 'Name')

# Can run checks on dataframe or list (below)
pks_1415_qc <- quality_check_base(pks_1415)


# Convert dataframe to list to work in 
ls <- sprintf('B%d', seq(1:14))  #gen list of batches to loop over
#ls <- sapply(seq(1:10), function(x){paste0('B', x)})  # alternate method

# split 2014 dataframe into list of df by batch
pks_1415_ls <- lapply(ls, function(x){filter(pks_1415, Batch == x)})

pks_1415_qc <- quality_check_base(pks_1415_ls)

do.call(rbind, lapply(pks_1415_qc, function(x){x[['duplicate_lipids']]}))  # Only M1M2 have dupl peaks
do.call(rbind, lapply(pks_1415_qc, function(x){x[['missing_lipids']]}))
do.call(rbind, lapply(pks_1415_qc, function(x){x[['lipid_frequency']]}))

# Create dataframe linking the methanol 13C to PLFA batches
# This is used in the isotope calculations
meth_1415 <- data.frame(Batch = ls, Meth13C = -55.84)

# Calculate corrected d13C for all samples. This can also be perfomed directly
# on a df (pks_1415)

d13c_1415 <- lapply(pks_1415_ls, 
                   function(x){x <- x[!x[['Name']] %in% c('13:0', '16:0', '19:0'), ]
                   correct_iso_base2(x, d13c_correction = 0, 
                                     meth_13c_df = meth_1415,
                                     min_height = 1)}
                        )
# Double check warnings thrown here: 

d13c_1415[[1]]

# Merge into dataframe
d13c_1415_df <- do.call('rbind', d13c_1415)


# for 2014/2015 data, important metadata mapped to sample names available in 
# 'Data for pivot tables' sheet in 
# 'SPRUCE DATA Nov 2016_all the final calculations.xls'
# 'Stats Ready Data Values' sheet in same xlsx file contains moisture and temp
# Can probably use Stats Ready Data Values to compare to my calcs

c_types1 <- c('text', 'text', 'text', 'text', 'numeric', 'numeric', 'text', 
              'numeric', 'text', 'text', 'numeric', 'numeric', 'numeric',
              'text', 'numeric')
# I think we should actually use this for initial data processing input, as 
# that would skip the step of importing and appending another df, but we don't
# need the weights for d13C calcs
# may also want to reorder the rename items
wts_1415 <- readxl::read_xls(path = paste0('C:\\Users\\Mark\\Dropbox\\',                                                                          'umn_gutknecht_postdoc\\',
                                           'spruce_project\\', 
                                           'plfa_13c\\data\\',
                                           'SPRUCE DATA Nov 2016_all the ',
                                           'final calculations.xls'), 
                             sheet = 'Data for pivot tables', na = 'NA',
                             col_types = c_types1) %>%
  dplyr::rename(BatchDataFileName = `ID_with Batch ID`, DataFileName = Name, 
                TotalPeakArea1 = `Total Peak Area`, RetTimeSecs = `RT (Sec)`,
                BatchID = `Batch ID`, SampleWt = `Weight (g)`, 
                PeakHeightnA = `Peak Height`, DisplayDelta1 = `Corrected 13C`)
# Cell M2261 is blank, gives warning about expecting numeric


# Import df with temp and moisture data for plotting with isotope data
md_1415 <- readxl::read_xls(path = paste0('C:\\Users\\Mark\\Dropbox\\',                                                                          'umn_gutknecht_postdoc\\',
                                             'spruce_project\\', 
                                             'plfa_13c\\data\\',
                                             'SPRUCE DATA Nov 2016_all the ',
                                             'final calculations.xls'), 
                               sheet = 'Stats Ready Data Values', na = 'NA',
                               range = cell_cols('A:I')) %>%
  dplyr::rename(BatchDataFileName = `Sample ID`, BatchID = `Batch ID`, 
                Temp = `Soil Probe Temp`, DepthIncrement = `Depth Increment`,
                GWC_perc = `% Gravimetric Water Content`) %>%
  dplyr::mutate(DataFileName = str_extract(BatchDataFileName, 
                                           '(?<=[A-Za-z][0-9]{1,2}_).+')) %>%
  dplyr:: select(-BatchDataFileName)  # indic_iso_base matches on DataFileName, this is dup
md_1415[['Month']] <- factor(md_1415[['Month']], 
                           levels = c('June', 'July', 'August', 
                                      'October'))


# Calculate d13C for bioindicator groups
ind13_1415 <- lapply(d13c_1415, 
                        function(x){x <- x[!x[['Name']] %in% c('13:0', '16:0', '19:0'), ]
                        indic_iso_base(x, md_1415)}
)

# Received errors stating that there were multiple rows of same biomarker name, so troubleshooting
# Looking for which batch, so I can see if this is an issue
ind13_1415 <- lapply(d13c_1415, 
                          function(x){x <- x[!x[['Name']] %in% c('13:0', '16:0', '19:0'), ]
                          tryCatch(indic_iso_base(x, md_1415),
                                   warning = function(c){
                                     message(paste0('Batch ', x[['Batch']][1],
                                                    ' contains duplicate', 
                                                    ' lipids'))})
                                   }  # If use print instead of message(), output gets incorporated into dataframe
)


ind13_1415[[3]][ind13_1415[[3]][['Name']] %in% c('22:1', '15:0', '15:0 anteiso', '16:1 w7c'), ]
# Batch 3 has multiple 22:1 peaks in B3_M1M2_1 half dilution_3.raw

ind13_1415[[10]][ind13_1415[[3]][['Name']] %in% c('22:1', '15:0', '15:0 anteiso', '16:1 w7c'), ]
# Batch 10 has duplicate 15:0 peaks in B10_M1M2 half dilution.raw

ind13_1415[[11]][ind13_1415[[3]][['Name']] %in% c('22:1', '15:0', '15:0 anteiso', '16:1 w7c'), ]
# Batch 11 has duplicate 16:1 w7c peaks in B11_M1M2 half dilution.raw
# Not sure where the 15:0 anteiso duplicate is

ind13_1415[[11]][ind13_1415[[3]][['Name']] =='16:1 w7c', ]


# Merge lists into single dataframe for each type of analysis
# This faciliates further analysis and visualization
ind13_1415_df <- do.call('rbind', ind13_1415)  # d13C of individual lipids


# Some samples with very low total lipids didn't have any lipids associated
# with specific microbial groups. This means that they remain in the df, but
# all of their d13C values end up as NA. May want to verify that this is what's 
# going on will all values.


d13c_1415_df <- filter_outliers(d13c_1415_df, cols = 'd13C_corrected', 
                               method = 'iqr')

ind13_1415_df <- filter_outliers(ind13_1415_df, cols = 'avg_d13C_corrected', 
                              method = 'iqr')


# All dates, plots, depths
ind13_1415_df %>% ggplot(aes(x = Temp, y = avg_d13C_corrected, color = Indicator)) +
  geom_point() +
  geom_smooth(aes(color = Indicator), 
              method = 'lm', alpha = 0.3) +
  facet_wrap(~Indicator)

d13c_1415_df %>% ggplot(aes(x = Temp, y = d13C_corrected, color = Plot)) +
  geom_point() +
  geom_smooth(aes(color = Indicator), 
              method = 'lm', alpha = 0.3) #+
  facet_wrap(~Indicator)

# Pre-DPH only
ind13_1415_df %>% 
  filter(Month == 'June' & Year == 2014) %>%
  ggplot(aes(x = Temp, y = avg_d13C_corrected, color = Indicator)) +
  geom_point() +
  geom_smooth(aes(color = Indicator), 
              method = 'lm', alpha = 0.3) +
  facet_wrap(~Indicator) + 
  ggtitle('June 2014')

# Post-DPH only
ind13_1415_df %>% 
  filter(Month == 'Sept' & Year == 2014) %>%
  ggplot(aes(x = Temp, y = avg_d13C_corrected, color = Indicator)) +
  geom_point() +
  geom_smooth(aes(color = Indicator), 
              method = 'lm', alpha = 0.3) +
  facet_wrap(~Indicator) + 
  ggtitle('June 2014')

# Equilibration (June 2015) only
ind13_1415_df %>% 
  filter(Month == 'June' & Year == 2015) %>%
  ggplot(aes(x = Temp, y = avg_d13C_corrected, color = Indicator)) +
  geom_point() +
  geom_smooth(aes(color = Indicator), 
              method = 'lm', alpha = 0.3) +
  facet_wrap(~Indicator) + 
  ggtitle('June 2014')


ref <- lipid_reference


# Still want to add acrotelm/mesotelm/catotelm designations (prob after merge)
# Add depth number and top depth/bottom depth

#####################
##### 2016 data #####
#####################

md_16 <- readxl::read_xlsx(path = paste0('C:\\Users\\Mark\\Dropbox\\',                                                                          'umn_gutknecht_postdoc\\',
                                           'spruce_project\\', 
                                           'plfa_13c\\data\\',
                                           '20190916_spruce_plfa_metadata.xlsx'), 
                             sheet = 'Sheet1', na = 'NA') %>%
  dplyr::rename(Batch = GCBatchNum, SampleWt = SampleWeight)
# Make Month a factor and set levels to display in proper order on plots
md_16[['Month']] <- factor(md_16[['Month']], 
                             levels = c('June', 'July', 'August', 
                                        'October'))

# Select directory holding batch files
pks_16_dir <- paste0('C:\\Users\\Mark\\Desktop\\', 
                     'temp plfa practice - can delete whenever')

# Load peak lists for all batches as a list of dataframes (actually tibbles)
pks_16_ls <- import_batch_multi_base(source_dir, keyword = 'Batch')  

pks_16_qc <- quality_check_base(pks_16_ls)

# Display dataframe of duplicate lipids in batch 7
# Need to address these duplicates
pks_16_qc[[1]]['duplicate_lipids']
pks_16__qc[[7]]['duplicate_lipids']

# Create dataframe with methanol lots used for each batch
# Need to resolve 13C correction for some batches --> no data for it
meth_16 <- data.frame(Batch = c(1, 2, 3, 4, 5, 6, 7), 
                          Meth13C = c(-0.74, -0.89, -0.74, -0.89, -0.74, -0.74, 
                                      -0.74))

# Calculate d13C for individual lipids
# Filter out 13:0, 16:0 and 19:0 lipids, as we know some are non-microbial 
# and there's no (easy) way to correct for that.

# 2016 data had pretty weak signals, so majority of peaks get filtered out if
# min_height = 1
d13c_16 <- lapply(pks_16_ls, 
                   function(x){x <- x[!x[['Name']] %in% c('13:0', '16:0', '19:0'), ]
                   correct_iso_base2(x, d13c_correction = 37.433, 
                                     meth_13c_df = meth_16,
                                     min_height = 1)}
)

# Merge list of df into single df for further analysis
# TO DO: should we/do we need to change indic_iso_base() to accept df?
d13c_16_df <- do.call(rbind, d13c_16)

# Calculate d13C for bioindicator groups
ind13_16 <- indic_iso_base(d13c_16_df, md_16)

#Filter outliers
d13c_all_df <- filter_outliers(d13c_all_df, cols = 'd13C_corrected', 
                               method = 'iqr')

indic_d13c <- filter_outliers(indic_d13c, cols = 'avg_d13C_corrected', 
                              method = 'iqr')

# These 13C values are uncorrected. Methanol d13C = -55.84. Need to double check which biomarkers included
l <- lipid_reference


# Base R option
#plfa_2014 <- transform(plfa_2014, Batch = regmatches(BatchDataFileName, 
#                                        gregexpr('B[0-9]+(?=_)', 
#                                                 BatchDataFileName, 
#                                                 perl = TRUE))[[1]])  # regmatches() returns a list, so need to subset




# base version. Could use str_extract() from stringr
regmatches('B1_M1M2_1 half dilution', gregexpr('B[0-9]+(?=_)', 'B1_M1M2_1 half dilution', perl = TRUE))

column_names <- c('Batch', 'DataFileName', 'RetTimeSecs', 'MajorHeightnA',
                  'TotalPeakArea1', 'DisplayDelta1', 'Name')
column_types <- c('text', 'text', 'numeric', 'numeric', 'numeric', 'guess',
                  'text')