#' Run Cox regression
#' @description Run Cox regression, tidy the output and add names of model, outcome and analysis.
#' @param cohort_steroids A list containing a dataset for a specified outcome and analysis.
#' @return A dataframe containing results.

analysis_regression <- function(cohort_steroids, slices_exposure_model) {

		 			coxph(formula(paste("outcome_surv ~", slices_exposure_model$exposure, slices_exposure_model$model)), 
		 						data=cohort_steroids$data, 
		 						cluster=pracid, 
		 						timefix=FALSE) %>% 
		 				tidy(exponentiate=TRUE, conf.int=TRUE) %>%
		 				mutate(outcome=cohort_steroids$outcome, 
		 							 exposure=slices_exposure_model$exposure,
		 							 model=slices_exposure_model$model_name,
		 							 analysis=cohort_steroids$analysis)
		 	
}



#' Run Cox regression for Competing risks models
#' @description Run Cox regression, tidy the output and add names of model, outcome and analysis.
#' @param cohort_steroids A list containing a dataset for a specified outcome and analysis.
#' @return A dataframe containing results.

analysis_regression_cr <- function(cohort_steroids_cr, slices_exposure_model) {
	
	coxph(formula(paste("outcome_surv ~", slices_exposure_model$exposure, slices_exposure_model$model)), 
				data=cohort_steroids_cr$data, 
				id=patid,
				cluster=pracid, 
				timefix=FALSE) %>% 
		tidy(exponentiate=TRUE, conf.int=TRUE) %>%
		mutate(outcome=cohort_steroids_cr$outcome, 
					 exposure=slices_exposure_model$exposure,
					 model=slices_exposure_model$model_name,
					 analysis=cohort_steroids_cr$analysis)
	
}



#' Calculate rates.
#' @description Calculate rates, person-years and number of fractures for different variables.
#' @param datasets AA lists of lists containing datasets for different outcomes.
#' @param rate_vars A character vector of variables for which rates should be calculated.
#'
#' @return A dataframe containing results.
#' @export
#'
#' @examples
analysis_rates <- function(datasets, rate_vars) {
	results_rates <- pyears(outcome_surv ~ eval(as.symbol(rate_vars)), data = datasets$data, scale = 365.25) %>%
		tidy() %>% 
		mutate(rate=(event/pyears)*1000) %>% 
		cbind(tibble(term=paste0(rate_vars, levels(as.factor(datasets$data[[rate_vars]]))), 
								 group=rate_vars, 
								 outcome=datasets$outcome,
								 analysis=datasets$analysis))
	#Get unexposed person years in the same row
	results_rates <- results_rates %>% 
		group_by(outcome, group) %>%
		mutate(pyears_unexposed=pyears[1],
					 ratio=1/(pyears[2]/pyears[1]),
					 ratio_text=paste0("1:", round(ratio))) %>%
		ungroup()
	
	results_rates
}



#' Calculate total and average follow-up times
#' @param cohort 
#' @return
analysis_futime <- function(cohort) {
	tibble(
		total = survival::pyears(cohort$outcome_surv ~ 1)$pyears,
		n_pats = length(unique(cohort$patid)),
		avg = total/n_pats,
		
		total_cont = survival::pyears(filter(cohort, pattern_rt_window=="0-89")$outcome_surv ~ 1)$pyears,
		n_pats_cont = length(unique(filter(cohort, pattern_rt_window=="0-89")$patid)),
		avg_cont = total_cont/n_pats_cont,
		
		total_int = survival::pyears(filter(cohort, pattern_rt_window=="90-180")$outcome_surv ~ 1)$pyears,
		n_pats_int = length(unique(filter(cohort, pattern_rt_window=="90-180")$patid)),
		avg_int = total_int/n_pats_int,
	)
}



#' Calculate participant counts at various stages
#' @param files 
#' @param extract 
#' @param extract_ocs 
#' @return
analysis_cohort_flow <- function(files, extract, extract_ocs) {
	db_pop <- files$patient %>% 
		open_dataset() %>% 
		summarise(n=n()) %>%
		collect() %>% 
		mutate(step="Database population")
	
	extract_pop <- extract %>% 
		count() %>% 
		mutate(step="with Eczema/Asthma/COPD code")
	
	extract_ocs_pop <- extract_ocs %>% 
		count() %>% 
		mutate(step="with >1 OCS prescription")
	
	bind_rows(db_pop, extract_pop, extract_ocs_pop)
	
}



#' Get baseline characteristics.
#' @description Calculate descriptive statistics at baseline using the gtsummary package.
#' @param cohort A dataframe containing the cohort of interest with only one row per participant.
#' @return A dataframe with baseline characteristics.
analysis_baseline_characteristics <- function(cohort) {
	library(gtsummary)
	library(tidyverse)
	
	cohort %>% 
		mutate(pattern_rt_window=case_when(pattern_rt_window=="90-180" ~ "low intensity",
																			 pattern_rt_window=="0-89" ~ "high intensity"),
					 carstairs=fct_explicit_na(fct_rev(carstairs))) %>% 
		select(age_at_index, sex, carstairs, eczema, asthma, 
					 copd, rheumatoid_arthritis, non_steroid_fx_drugs, 
					 timegroup, pattern_rt_window) %>% 
		tbl_summary(by="pattern_rt_window",
								
								statistic = list(where(is.numeric) ~ "{median} ({p25} - {p75})",
																 where(is.factor) ~ "{n} ({p}%)"),
								
								value = list(sex ~ "Male", eczema ~ 1, asthma ~ 1, copd ~ 1, 
														 rheumatoid_arthritis ~ 1, non_steroid_fx_drugs ~ 1),
								
								label = list(age_at_index ~ "Age", carstairs ~ "Deprivation", sex ~ "Male",
														 eczema ~ "Eczema", asthma ~ "Asthma", copd ~ "COPD",
														 rheumatoid_arthritis ~ "Rheumatoid arthritis",
														 non_steroid_fx_drugs ~ "Other drugs increasing fracture risk",
														 timegroup ~ "Year of indexdate"),
							
									digits = list(all_categorical() ~ c(0, 1))
								) %>% 
		as_tibble()
	
}



#' Get counts of exposed and unexposed in each rolling window.
#' @param cohort_steroids 
#' @param exposure 
#' @return A list of dataframes
analysis_exposed_counts <- function(cohort_steroids, exposure) {
	cohort_steroids$data %>% 
		group_by(patid, rolling_window) %>% 
		filter(!duplicated(rolling_window)) %>% 
		ungroup() %>% 
		select(rolling_window, all_of(exposure)) %>% 
		table() %>% 
		as.data.frame() %>% 
		mutate(rolling_window=as.numeric(rolling_window),
					 outcome=cohort_steroids$outcome, 
					 analysis=cohort_steroids$analysis,
					 exposure=exposure)
}



