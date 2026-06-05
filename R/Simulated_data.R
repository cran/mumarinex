#' Simulated_data
#'
#' @name Simulated_data
#' @docType data
#' @author N. Chauvel
#'
#' @usage data(Simulated_data)
#'
#' @format A data frame with 50 rows (samples) and 12 species:
#' \describe{
#'   \item{Sp_A}{Species tolerant to different impacts, mean abundance = 500}
#'   \item{Sp_B}{Species tolerant to different impacts, mean abundance = 50}
#'   \item{Sp_C}{Species highly sensitive to impacts at stations RD, RDI, and M, leading to their disappearance, mean abundance = 500}
#'   \item{Sp_D}{Species highly sensitive to impacts at stations RD, RDI, and M, leading to their disappearance, mean abundance = 50}
#'   \item{Sp_E}{Species favored by impacts at stations RI, RDI, and M, leading to their appearance, mean abundance = 500}
#'   \item{Sp_F}{Species favored by impacts at stations RI, RDI, and M, leading to their appearance, mean abundance = 50}
#'   \item{Sp_G}{Species sensitive to impacts at stations AD, ADI, and M, leading to decreased abundance, mean abundance = 1000}
#'   \item{Sp_H}{Species sensitive to impacts at stations AD, ADI, and M, leading to decreased abundance, mean abundance = 500}
#'   \item{Sp_I}{Species stimulated by impacts at stations AI, ADI, and M, leading to increased abundance, mean abundance = 1000}
#'   \item{Sp_J}{Species stimulated by impacts at stations AI, ADI, and M, leading to increased abundance, mean abundance = 500}
#'   \item{Sp_K}{Species strongly stimulated by impacts at stations D and M, leading to a shift in dominance, mean abundance = 1000}
#'   \item{Sp_L}{Species strongly stimulated by impacts at stations D and M, leading to a shift in dominance, mean abundance = 500}
#' }
#'
#' @description
#' This dataset was constructed to simulate the theoretical impacts of a disturbance on the abundances of 12 hypothetical species.
#' Reference stations REF1 and REF2 represent control sites, with alternating taxa showing normally distributed abundances (means of 500 and 50, standard deviation of 10% of mean).
#' Stations RI, RD and RDI illustrate changes in community richness, represented by the loss of species C and D (station RD), the appearance of species E and F (station RI), or both simultaneously (station RDI).
#' Stations AD, AI and ADI represent changes in community structure, expressed as a decrease in the abundance of species G and H (station AD), an increase in the abundance of species I and J (station AI), or both (station ADI).
#' Station D focuses on community dominance shifts, characterized by a strong increase in the abundance of species K and L.
#' Station M1 combines all of these effects (loss, gain, structural shifts, and dominance changes).
#' Species A and B are insensitive to impacts and act as tolerant taxa.
#'
NULL
