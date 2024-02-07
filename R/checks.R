#' Make table of distribution of prescription patterns (high or low intensity) by exposure definition
#' @param flat Dataset with one row per participant.
#' @return A dataframe with distribution of prescription patterns
checks_exposure_dist <- function(flat) {
	library(janitor)
	library(tidyverse)
	
	temp <- flat %>% 
		select(pattern_rt_window, pattern_rt_window0,
					 pattern_rt_sum_gaps, pattern_rt_sum_gaps0,
					 pattern_rt_n_gaps, pattern_rt_n_gaps0)
	
	map2_df(temp, names(temp), ~tibble(.y, tabyl(.x, show_na=FALSE), intensity=c("low intensity", "high intensity")))
}



#' Histograms of continuous measures of prescription pattern
#' @param flat 
#' @return A list of ggplot objects
checks_histograms <- function(flat) {
	library(ggplot2)
	library(tidyverse)
	plot_hist <- function(x, label) {ggplot(data=flat) +
			geom_histogram(aes_string(x=x), 
										 binwidth = 1,
										 color="#0D5257",
										 fill="#0D5257") +
			theme_bw() +
			labs(x=label)}
	
	list(
	plot_hist("riskthreshold_window_length", 
						"Days taken to reach the risk threshold of 450mg PED"),
	plot_hist("riskthreshold_sum_gaps", 
						"Sum of gap days before reaching the risk threshold of 450mg PED"),
	plot_hist("riskthreshold_n_gaps", 
						"Number of gaps before reaching the risk threshold of 450mg PED")
	)
}


#' Make survfits
#' @param cohort_steroids 
#' @param slices_outcome_analysis 
#' @return A list of survfits
checks_survfits <- function(cohort_steroids, slices_outcome_analysis) {
	library(survival)
	library(survminer)
	library(tidyverse)
	# fit model for FPC as outcome
	data_bp <- cohort_steroids[[with(slices_outcome_analysis, which(outcome == "fpc" & analysis == "main"))]]$data
	fit_bp <- survfit(outcome_surv ~ pattern_rt_window, data=data_bp, id=patid, timefix = FALSE)

	# fit model for fracture as outcome
	data_fx <- cohort_steroids[[with(slices_outcome_analysis, which(outcome == "fract_composite" & analysis == "sens_all_follow_up"))]]$data
	fit_fx <- survfit(outcome_surv ~ pattern_rt_window, data=data_fx, id=patid, timefix = FALSE)

	list(data_bp, fit_bp, data_fx, fit_fx)
	
}

#' Make survplots
#' @param survfits 
#' @return A list of survplots
checks_survplots <- function(survfits) {
	library(survminer)
	survplot_bp <- ggsurvplot(survfits[[2]], data=survfits[[1]],
														legend.labs = c("low intensity OCS use", "high intensity OCS use"),
														ylim=c(0.5,1))
	survplot_fx <- ggsurvplot(survfits[[4]], data=survfits[[3]],
														legend.labs = c("low intensity OCS use", "high intensity OCS use"),
														ylim=c(0.5,1))
	list(survplot_bp, survplot_fx)
	
}


#' Test the proportional hazards assumption for Cox models
#' @param cohort_steroids 
#' @param slices_outcome_analysis 
#' @return A list of cox.zph objects
checks_cox_zph <- function(cohort_steroids, slices_outcome_analysis) {
	# fit model for FPC as outcome
	cox_zph0 <- cohort_steroids[[with(slices_outcome_analysis, which(outcome == "fpc" & analysis == "main"))]]$data %>% 
		coxph(outcome_surv ~ pattern_rt_window, data=., id=patid, timefix=FALSE) %>% 
		cox.zph()
	
	# fit model for fracture as outcome
	cox_zph1 <- cohort_steroids[[with(slices_outcome_analysis, which(outcome == "fract_composite" & analysis == "sens_all_follow_up"))]]$data %>% 
		coxph(outcome_surv ~ pattern_rt_window, data=., id=patid, timefix=FALSE) %>% 
		cox.zph()

	list(cox_zph0, cox_zph1)
}
