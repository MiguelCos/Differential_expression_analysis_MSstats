## Differential expression analysis with MSstats (LFQ) ----

# Load required packages ----

library(tidyverse)
library(MSstatsTMT)
library(MSstats)
library(here)

## Load the required data 

lf_evid <- read.delim(file = here::here("Data/Label-free/evidence.txt"),
                      sep = "\t", stringsAsFactors = FALSE)

lf_protg <- read.delim(file = here::here("Data/Label-free/proteinGroups.txt"),
                       sep = "\t", stringsAsFactors = FALSE)

lf_annot <- read.delim(file = here::here("Data/Label-free/ABRF2015_MaxQuant_annotation.csv"),
                       sep = ",", stringsAsFactors = FALSE)

## Run `MaxQtoMSstatsFormat()` (label-free analysis)  

if (file.exists(here::here("Data/Label-free/LF_ABRF_data_formated.Rds")) == FALSE){
  
  lf_mstsformat <- MaxQtoMSstatsFormat(evidence = lf_evid,
                                       proteinGroups = lf_protg,
                                       annotation = lf_annot)
  
  
  
  saveRDS(lf_mstsformat,
          file = "Data/Label-free/LF_ABRF_data_formated.Rds")
} else {
  lf_mstsformat <- readRDS(file = "Data/Label-free/LF_ABRF_data_formated.Rds")
}

## Run `dataProcess()` (label-free analysis)  

if (file.exists(here::here("Data/Label-free/LF_ABRF_data_formated.Rds")) == FALSE){
  
  lf_proc <- dataProcess(lf_mstsformat)
  
  saveRDS(lf_proc,
          file = "Data/Label-free/LF_ABRF_data_processed.Rds")
} else {
  lf_proc <- readRDS("Data/Label-free/LF_ABRF_data_processed.Rds")
}


## Quality control plots  

dataProcessPlots(data.peptide = tmt_mstsformat,
                    data.summarization = tmt_procs,
                    type = "QCPlot",
                    address = FALSE,
                    originalPlot = FALSE,
                    summaryPlot = TRUE,
                    which.Protein = "0000",
                    width = 21,
                    height = 5)

## Profile Plots  
dataProcessPlots(data.peptide = tmt_mstsformat,
                    data.summarization = tmt_procs,
                    type = "ProfilePlot",
                    address = FALSE,
                    originalPlot = FALSE,
                    summaryPlot = TRUE,
                    which.Protein = "00000",
                    width = 21,
                    height = 5)

## Creating a contrast matrix  

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)


comparison_all <- rbind(comparison1, 
                        comparison2, comparison3, 
                        comparison4, comparison5, comparison6)

colnames(comparison_all)<- c("0.125", "0.5", "0.667", "1")

row.names(comparison_all)<-c("0.5-0.125","0.667-0.125","1-0.125",
                             "0.667-0.5","1-0.5","1-0.667")

## Executing `groupComparison()`

diffexpr <- groupComparisonTMT(data = lf_proc,
                               contrast.matrix = comparison_all)
