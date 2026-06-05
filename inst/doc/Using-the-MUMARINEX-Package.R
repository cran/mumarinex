## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(mumarinex)

## -----------------------------------------------------------------------------
# Load example dataset
data("Simulated_data")

# Display the first rows
head(Simulated_data)

# Definition of the reference position
ref_idx <- 1:10 # row number of the reference samples

## -----------------------------------------------------------------------------
# Compute MUMARINEX and sub-indices
rMUM <- mumarinex(x = Simulated_data, ref = ref_idx, subindices = TRUE, log = FALSE)

# Extract MUMARINEX
rMUMARINEX<-rMUM$MUMARINEX

# Extract sub-indices
Subind<-rMUM$Subindices

## ----fig.width=7, fig.height=5------------------------------------------------
stations<-matrix(unlist(strsplit(rownames(Simulated_data),".",fixed=TRUE)),ncol=2,byrow=TRUE)[,1] # get station labels from data rownames

stations<-factor(stations,levels=unique(stations)) # setting station names as factor to specify in which order it must display it in the boxplot

boxplot(rMUMARINEX~stations,ylim=c(0,1)) # ylim is set in the interval 0-1 as it is the range of MUMARINEX

## ----fig.width=10, fig.height=5-----------------------------------------------

decomplot(x = Simulated_data, g = stations, ref = ref_idx, log = FALSE, main = "Artificial data")

## -----------------------------------------------------------------------------
diagnostic_tool(x = Simulated_data, g = stations, ref = ref_idx, log = FALSE)

