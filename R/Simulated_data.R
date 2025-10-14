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
#'   \item{Sp_A}{Species tolerant to different impacts, mean abundance = 50}
#'   \item{Sp_B}{Species tolerant to different impacts, mean abundance = 10}
#'   \item{Sp_C}{Species highly sensitive to impacts at stations A1, A3, and M1, leading to their disappearance, mean abundance = 50}
#'   \item{Sp_D}{Species highly sensitive to impacts at stations A1, A3, and M1, leading to their disappearance, mean abundance = 10}
#'   \item{Sp_E}{Species favored by impacts at stations A2, A3, and M1, leading to their appearance, mean abundance = 50}
#'   \item{Sp_F}{Species favored by impacts at stations A2, A3, and M1, leading to their appearance, mean abundance = 10}
#'   \item{Sp_G}{Species sensitive to impacts at stations S1, S3, and M1, leading to decreased abundance, mean abundance = 50}
#'   \item{Sp_H}{Species sensitive to impacts at stations S1, S3, and M1, leading to decreased abundance, mean abundance = 10}
#'   \item{Sp_I}{Species stimulated by impacts at stations S2, S3, and M1, leading to increased abundance, mean abundance = 50}
#'   \item{Sp_J}{Species stimulated by impacts at stations S2, S3, and M1, leading to increased abundance, mean abundance = 10}
#'   \item{Sp_K}{Species strongly stimulated by impacts at stations D1 and M1, leading to a shift in dominance, mean abundance = 50}
#'   \item{Sp_L}{Species strongly stimulated by impacts at stations D1 and M1, leading to a shift in dominance, mean abundance = 10}
#' }
#'
#' @description
#' This dataset was constructed to simulate the theoretical impacts of a disturbance on the abundances of 12 hypothetical species.
#' Stations R1 to R3 illustrate changes in community richness, represented by the loss of species C and D (station R1), the appearance of species E and F (station R2), or both simultaneously (station R3).
#' Stations S1 to S3 represent changes in community structure, expressed as a decrease in the abundance of species G and H (station S1), an increase in the abundance of species I and J (station S2), or both (station S3).
#' Station D1 focuses on community dominance shifts, characterized by a strong increase in the abundance of species K and L.
#' Station M1 combines all of these effects (loss, gain, structural shifts, and dominance changes).
#' Reference stations REF1 and REF2 represent control sites, with alternating taxa showing normally distributed abundances (means of 50 and 10, standard deviation of 10% of mean).
#' Species A and B are insensitive to impacts and act as tolerant taxa.
#'
NULL
