library(targets)
library(tarchetypes)
# This _targets.R file defines the {targets} pipeline.
# Run tar_make() to run the pipeline, tar_make(target) to run up to a defined target, and tar_read(target) to view the results.

# # Parallel computing with the `future` framework (run branches of targets on multiple cores)
# library(future.callr)
# plan(callr)

# Define external paths
# Paths to the following need to be defined:
# the raw data (path_in)
# the linked Carstairs data (path_z_drive_linked)
# prescription dosage infromation (path_common_dosages)
source("paths/paths.R")

# Source all functions from the "R" folder
sapply(list.files("R", full.names = TRUE) ,source, .GlobalEnv)

# Set target-specific options such as packages.
tar_option_set(
	packages = c(
		"rio", #For reading various filetypes
		"arrow", #For reading parquet files
		"cli", #For producing command line interface elements such as loading bars
		"zoo", # For rolling windows
		"survival", # For survival analysis including Cox regression
		"magrittr", # To improve readability with pipes
		"haven", # To import Stata .dta
		"readstata13", # To import Stata .dta
		"broom", # To clean regression output
		"epiDisplay", # For data exploration
		"biostat3", # For data exploration
		"lubridate", # To manage dates
		"summarytools", # For data exploration
		"gt", # To create publication-ready tables
		"ggplot2", # To make plots
		"forestplot", # To make forest plots
		"DiagrammeR", # To make flow diagrams
		"DiagrammeRsvg", # To export flow diagrams
		"rsvg", # To export flow diagrams
		"dtplyr", #data.table backend for dplyr
		"qs", #To store targets in qs format for quicker reading and writing
		"tidyverse", # For data management
		"collapse" #For fast data management
	),
	format="qs" #Save all targets in qs format, which should be quicker to read and write than the default RDS
) 



# List of target objects.
list(
	
	# Specifications -------------------------------------------------------------
	
	# ┠ Outcomes -----
	tar_target(
		outcome, 
		head(n=99, #TEMPORARILY LIMIT TO SPEED UP RUNTIME
				 c("fpc",
				 	"bisphosphonate",
				 	"fract_composite", 
				 	"fract_hip",
				 	"calcium_and_vit_d",
				 	"dxa_scan",
				 	"anxiety_meds",
				 	"epilepsy_meds",
				 	"migraine_meds"
				 ))
	),
	
	# ┠ Exposures -----
	tar_target(
		exposure,
		head(n=99, #TEMPORARILY LIMIT TO SPEED UP RUNTIME
				 c("pattern_rt_window",
				 	"pattern_rt_sum_gaps", 
				 	"pattern_rt_n_gaps",
				 	"pattern_rt_window0",
				 	"pattern_rt_sum_gaps0", 
				 	"pattern_rt_n_gaps0",
				 	"logdays"
				 ))
	),
	
	# ┠ Imputation -----
	#Specify under which conditions the daily glucocorticoid dose should be imputed
	tar_target(
		imputation_condition,
		c("is.na(daily_dose) & !is.na(impute)", #Impute for prescriptions where daily dose is missing
			"(is.na(daily_dose) | daily_dose==0) & !is.na(impute)") #Impute for prescriptions where daily dose is missing or 0
	),
	
	# ┠ Analyses -----
	tar_target(
		analysis,
		head(n=99, #TEMPORARILY LIMIT TO SPEED UP RUNTIME
		c("main", #Main analysis
			"sens_all_follow_up", #Sensitivity analysis: Follow up not limited to 1 year
			"sens_asthma", #Sensitivity analysis: Exclude people who don't have asthma at baseline
			"sens_eczema", #Sensitivity analysis: Exclude people who don't have eczema at baseline
			"sens_copd", #Sensitivity analysis: Exclude people who don't have copd at baseline
			"sens_all_ages", #Sensitivity analysis: don't restrict to only those aged 66+ when they cross the risk threshold for the first time
			"sens_all_ped_imputed" #Sensitivity analysis: Impute all values for PED
			) 
		)
	),
	
	
	# ┠ Models -----
	tar_target(# Specify models for regression
		model, 
		head(n=99, #TEMPORARILY LIMIT TO SPEED UP RUNTIME
		c("crude" = "",
			"adjusted_comorb" = "+ age_at_index + sex + carstairs + eczema + asthma + copd + rheumatoid_arthritis",
			"adjusted_age_sex_carstairs" = "+ age_at_index + sex + carstairs"
			#"adjusted_comorb_cumdose" = "+ age_at_index + sex + carstairs + asthma + copd + rheumatoid_arthritis + cumdose_cont",
			#"adjusted_comorb_non_steroid_fx_drugs" = "+ age_at_index + sex + carstairs + asthma + copd + rheumatoid_arthritis + non_steroid_fx_drugs"
			))
	),
	
	# ┠ Rate vars -----
	tar_target( # Specify variables for which rates should be calculated
		rate_vars, 
		c(exposure, "eczema", "asthma", "copd")
	),
	
	# ┠ Main inputs -----
	tar_target(# Specify the main outcomes and analyses here (i.e. these will be run in combination with all models, exposures)
		slices_outcome_analysis,
		crossing(outcome, analysis) %>% 
			filter(outcome %in% c("fpc", "fract_composite", "bisphosphonate") | 
						 	analysis %in% c("main"))
	),
	tar_target(# Specify the main models and exposures here (i.e. these will be run in combination with all outcomes, analyses)
		slices_exposure_model,
		crossing(exposure, model) %>% 
			filter(exposure %in% c("pattern_rt_window") | 
						 	model %in% c("", "+ age_at_index + sex + carstairs + eczema + asthma + copd + rheumatoid_arthritis")) %>% 
			mutate(model_name=names(model))
	),

	
	
	# Data paths -----------------------------------------------------------------
	# External linked data
	tar_file(path_carstairs_pat, paste0(path_z_drive_linked, "patient_carstairs_20_051.txt")),
	tar_file(path_carstairs_prac, paste0(path_z_drive_linked, "practice_carstairs_20_051.txt")),
	# Other
	tar_file(path_ped_reference, "input/ped_reference.csv"),
	
	
	
	
	
	
	# Codelists ---------------------------------------------------------
	#Codelists must be in the correct folder as .csv files and have the same name as specified here
	
	tar_target( # Codelists used to extract the cohort
		codelists_extract,
		tribble(~name, ~codevar, ~extract_from,
						"asthma",	"medcode",	"clinical",
						"copd",	"medcode",	"clinical",
						"eczema",	"medcode",	"clinical") %>% 
			mutate(path=paste0("codelists/", codevar, "/", name, ".csv"),
						 full=map(path, ~read_csv(.x)),
						 codes=map2(path, codevar, ~read_csv(.x)[[.y]]))
			),
		
	tar_target( # Codelists used to extract eventdata
		codelists, 
		tribble(~name, ~codevar, ~extract_from,
						#medcodes
						"asthma",	"medcode",	"clinical",
						"cancer",	"medcode",	"clinical",
						"copd",	"medcode",	"clinical",
						"dxa_scan",	"medcode",	"clinical_referral",
						"eczema",	"medcode",	"clinical",
						"fractures_all",	"medcode",	"clinical",
						"fractures_hip",	"medcode",	"clinical",
						"fractures_pelvis",	"medcode",	"clinical",
						"fractures_spine",	"medcode",	"clinical",
						"fractures_wrist",	"medcode",	"clinical",
						"rheumatoid_arthritis",	"medcode",	"clinical",
						# prodcodes
						"anxiety_meds",	"prodcode",	"therapy_reduced",
						"bisphosphonates",	"prodcode",	"therapy_reduced",
						"calcium_and_vit_d",	"prodcode",	"therapy_reduced",
						"epilepsy_meds",	"prodcode",	"therapy_reduced",
						"migraine_meds",	"prodcode",	"therapy_reduced",
						"non_steroid_fx_drugs",	"prodcode",	"therapy_reduced",
						"oral_glucocorticoids",	"prodcode",	"therapy") %>% 
			mutate(path=paste0("codelists/", codevar, "/", name, ".csv"),
						 full=map(path, ~read_csv(.x)),
						 codes=map2(path, codevar, ~read_csv(.x)[[.y]]))
	),
	
	
	
	# Extract -----------------------------------------------------------------

	tar_target( # Define search patterns for raw files
		patterns,
		c(clinical_referral = "Clinical[0-9]{1,99}_reduced.parquet|Referral[0-9]{1,99}_reduced.parquet",
			therapy = "Therapy[0-9]{1,99}.parquet",
			clinical = "Clinical[0-9]{1,99}.parquet",
			additional = "Additional[0-9]{1,99}.parquet",
			consultation = "Consultation[0-9]{1,99}.parquet",
			immunisation = "Immunisation[0-9]{1,99}.parquet",
			patient = "Patient[0-9]{1,99}.parquet",
			referral = "Referral[0-9]{1,99}.parquet",
			test = "Test[0-9]{1,99}.parquet",
			therapy_reduced = "Therapy[0-9]{1,99}_reduced.parquet")
	),
	
	tar_target( # Make list where each element contains the file paths of one of the above specified file types
		files,
		map(patterns, ~dir(path_in, pattern = .x, full.names = TRUE)),
		iteration = "list"
	),
	
	tar_target( # Extract a cohort of people that have a code in the files
		extract,
		extract_cohort(
			extract_files=files$clinical_referral,
			extract_codes=unlist(codelists_extract$codes),
			study_start="1998-01-02",
			study_end="2020-01-31",
			practice_info_file=paste0(path_in, "Practice.parquet"),
			patient_info_files=files$patient)
		),
	
	tar_target( # Get eventdata for every codelist
		eventdata,
		open_dataset(files[[codelists$extract_from]]) %>% 
			select(any_of(c("patid", "eventdate", "medcode", "prodcode", "dosageid", "qty", "numdays"))) %>% 
			filter(eventdate>=as.Date("1998-01-02"),
						 eventdate<=as.Date("2020-01-31"),
						 patid %in% extract$patid,
						 eval(sym(codelists$codevar)) %in% unlist(codelists$codes)) %>% 
			collect(),
		pattern = map(codelists),
		iteration = "list"
	),
	
	
	# Cohort specific filtering -----------------------------------------------

	# Get people who have at least 1 prescriptions for an oral corticosteroid (OCS)
	tar_target(
		extract_ocs,
			lazy_dt(eventdata[[which(codelists$name=="oral_glucocorticoids")]]) %>% # Take the eczema Rx eventdata
			group_by(patid) %>% #For each patient ...
			filter(!duplicated(eventdate)) %>% #...keep only one prescription per day...
			slice(1) %>%  #...and keep only the first prescription
			select(patid, eventdate) %>% 
			left_join(extract, by="patid") %>%  #Join the extract...
			mutate(eventdate=pmax(eventdate.x, eventdate.y)) %>% #...and keep the later eventdate
			select(-c(eventdate.x, eventdate.y)) %>% 
			as_tibble()
	),
	
	# Time-based filtering of cohort 
	#Startdate: latest of: date of any of the diseases from the extract, date of OCS prescription, 18th birthday, study start date (January 2, 1998), date their practice met CPRD quality-control standards, or practice registration date plus 1 year.
	#Enddate: earliest of: date of outcome, date of death, date the participant left the practice, date of the last data collection from the practice, or the end of the study (31st January 2020).
	tar_target(
		extract_ocs_filtered,
		lazy_dt(extract_ocs) %>% 
			mutate(adult_date=as.Date(paste0(realyob+18, "-06-01")),
						 indexdate=pmax(adult_date, eventdate, crd+365, uts, as.Date("1998-01-01"), na.rm = TRUE), 
						 enddate=pmin(tod, deathdate, lcd, as.Date("2020-01-31"), na.rm = TRUE)) %>% 
			filter(enddate>indexdate) %>% 
			as_tibble()
	),

	
	
	
	
	# Data management ------------------------------------------------------------
	
	# ┠ Prescriptions -----
	tar_target(
		prescriptions,
		eventdata[[which(codelists$name=="oral_glucocorticoids")]] %>% 
			filter(patid %in% main_cohort$patid) %>% 
			left_join(main_cohort[c("patid", "sex", "realyob", "indexdate")], by="patid") %>% 
			left_join(codelists$full[[which(codelists$name=="oral_glucocorticoids")]], by="prodcode") %>% 
			left_join(read.delim(path_common_dosages), by="dosageid") %>% 
			select(patid, eventdate, qty, drugsubstancename, daily_dose, dosage_text, substancestrength, sex, realyob, indexdate) %>% 
			mutate(eventdate=ymd(eventdate), pracid=patid %% 1000, dob=date(paste0(realyob, "-01-01")),
						 age_at_event=(as.numeric(eventdate)-as.numeric(dob))/365.25,
						 age_grp=cut_width(age_at_event, 5), qty=as.numeric(qty))
	),
	
	tar_target(
		prescriptions_ped,
		create_prescriptions_ped(prescriptions, path_ped_reference, imputation_condition),
		pattern = head(map(imputation_condition), n=99), #TEMPORARILY LIMIT TARGET TO SPEED UP RUNTIME
		iteration = "list"
	),
	
	tar_target(
		prescriptions_risk_threshold,
		create_prescriptions_risk_threshold(prescriptions_ped),
		pattern = map(prescriptions_ped),
		iteration = "list"
	),
	
	tar_target(
		prescriptions_expanded,
		create_prescriptions_expanded(prescriptions_ped),
		pattern = map(prescriptions_ped),
		iteration = "list"
	),
	
	tar_target(
		non_steroid_fx_drugs,
		eventdata[[which(codelists$name=="non_steroid_fx_drugs")]] %>%
			filter(patid %in% main_cohort$patid) %>%
			collect() %>% 
			left_join(main_cohort[c("patid", "indexdate")], by="patid") %>% 
			filter(eventdate >= indexdate) %>% #Keep only the first after indexdate for each patient
			arrange(eventdate) %>% 
			filter(!duplicated(patid)) %>% 
			select(patid, eventdate)
	),
	
	# ┠ Study cohort -----
	tar_target(
		main_cohort, 
		extract_ocs_filtered %>% 
			left_join(read.delim(path_carstairs_pat) %>% mutate(patid=as.numeric(str_sub(as.character(patid), end=-4)))) %>% 
			left_join(read.delim(path_carstairs_prac) %>% rename(carstairs2011_5_prac=carstairs2011_5)) %>% 
			mutate(carstairs=ifelse(!is.na(carstairs2011_5), carstairs2011_5, carstairs2011_5_prac)) %>% 
			rename(sex=gender) %>% 
			select(patid, indexdate, enddate, sex, realyob, pracid, carstairs)
	),
	

	
	tar_target(
		cohort_tmerged, 
		create_cohort_tmerged(main_cohort, 
												 prescriptions_risk_threshold,
												 prescriptions_expanded,
												 non_steroid_fx_drugs,
												 eventdata,
												 codelists),
		pattern = map(prescriptions_risk_threshold, prescriptions_expanded),
		iteration = "list"
	),

	tar_target(
		cohort_steroids,
		create_cohort_steroids(cohort_tmerged, 
													 slices_outcome_analysis, 
													 eventdata_cancer=eventdata[[which(codelists$name=="cancer")]]),
		pattern = map(slices_outcome_analysis),
		iteration = "list"
	), 
	
	# ┠ Additional  -----
	tar_target(
		flat,
		cohort_steroids[[1]]$data %>% filter(!duplicated(patid))
	),
	tar_target(
		baseline_characteristics,
		analysis_baseline_characteristics(flat)
	),

	
	
	
	
	
	# Analysis -------------------------------------------------------------------
	tar_target(
		results_regression, 
		analysis_regression(cohort_steroids, slices_exposure_model),
		pattern = cross(cohort_steroids, slices_exposure_model)
		),
	tar_target(
		results_rates, 
		analysis_rates(datasets=cohort_steroids, rate_vars=rate_vars),
		pattern = cross(cohort_steroids, rate_vars)
		),
	tar_target(
		baseline_characteristics_cohort_steroids, 
		analysis_baseline_characteristics(cohort=cohort_steroids[[1]]$data)
	),
	tar_target(
		exposed_counts,
		analysis_exposed_counts(cohort_steroids, exposure),
		pattern = cross(cohort_steroids, exposure)
	),
	tar_target(
		cohort_flow,
		analysis_cohort_flow(files, extract, extract_ocs)
	),
	tar_target(
		participant_flow,
		map(cohort_steroids, magrittr::extract, c("counts", "outcome", "analysis")) %>%
			map_df(~ .x$counts %>% mutate(outcome=.x$outcome, analysis=.x$analysis))
	),
	tar_target(
		futime,
		analysis_futime(cohort_steroids$data) %>% 
			mutate(outcome=cohort_steroids$outcome, analysis=cohort_steroids$analysis),
		pattern = cohort_steroids
	),

	
	

	# Multistate models ------------------------------------------------
	
	tar_target(
		cohort_steroids_ms,
		create_cohort_steroids_ms(cohort_tmerged=cohort_tmerged[[1]], 
															analysis, 
															eventdata_cancer=eventdata[[which(codelists$name=="cancer")]],
															outcome="event"),
		pattern = slice(map(analysis), 1:2),
		iteration = "list"
	),
	tar_target(
		results_regression_ms, 
		coxph(formula(paste("outcome_surv ~", "pattern_rt_window")), 
					data=cohort_steroids_ms$data, 
					id=patid, 
					timefix=FALSE) %>% 
			tidy(exponentiate=TRUE, conf.int=TRUE) %>%
			mutate(outcome=cohort_steroids_ms$outcome, 
						 exposure="pattern_rt_window",
						 model="crude",
						 analysis=cohort_steroids_ms$analysis),
		pattern = map(cohort_steroids_ms)
	),
	

	# Checks ------------------------------------------------------------------

	#Make lists of the most common codes for each eventdata
	tar_target(
		common_codes,
		map2(eventdata, codelists$full,
				 ~.x %>% 
				 	rename(code=matches("medcode|prodcode")) %>% 
				 	left_join(
				 		.y %>% rename(code=matches("medcode|prodcode")),
				 		by="code") %>% 
				 	group_by(across(matches("readterm|desc|productname"))) %>% 
				 	tally() %>% 
				 	arrange(desc(n)) %>% 
				 	filter(n>10)) %>% 
			set_names(codelists$name)
	),
	
	tar_target(exposure_dist, checks_exposure_dist(flat)),
	tar_target(histograms, checks_histograms(flat)),
	tar_target(survfits, checks_survfits(cohort_steroids, slices_outcome_analysis)),
	tar_target(survplots, checks_survplots(survfits)),
	tar_target(cox_zph, checks_cox_zph(cohort_steroids, slices_outcome_analysis)),
	
	
	# # Redaction --------------------------------------------------------------------
	# 
	# tar_target(
	# 	results_rates_redacted,
	# 	results_rates %>% 
	# 		mutate(across(c(pyears, n, event, pyears_unexposed), ~if_else(.x>=10, .x, NA_real_)))
	# ),
	# tar_target(
	# 	exposed_counts_redacted,
	# 	exposed_counts %>% 
	# 		mutate(across(c(Freq), ~if_else(.x>=10, .x, NA_integer_)))
	# ),
	# tar_target(
	# 	participant_flow_redacted,
	# 	participant_flow %>% 
	# 		mutate(n=as.double(n), pyrs=as.double(pyrs)) %>% 
	# 		mutate(across(c(n, pyrs), ~if_else(.x>=10, .x, NA_real_)))
	# ),
	

	# Write outputs -----------------------------------------------------------
	
	#Tables as .csv
	tar_file(exposed_counts_file, write_and_return_path(exposed_counts)),
	tar_file(results_rates_file, write_and_return_path(results_rates)),
	tar_file(participant_flow_file, write_and_return_path(participant_flow)),
	tar_file(results_regression_file, write_and_return_path(results_regression)),
	tar_file(results_regression_ms_file, write_and_return_path(results_regression_ms)),
	tar_file(cohort_flow_file, write_and_return_path(cohort_flow)),
	tar_file(baseline_characteristics_file, write_and_return_path(baseline_characteristics)),
	tar_file(exposure_dist_file, write_and_return_path(exposure_dist)),
	tar_file(slices_outcome_analysis_file, write_and_return_path(slices_outcome_analysis)),
	tar_file(futime_file, write_and_return_path(futime)),
	
	#Vectors as text
	tar_file(exposure_file, write_lines_and_return_path(exposure)),
	tar_file(outcome_file, write_lines_and_return_path(outcome)),
	tar_file(analysis_file, write_lines_and_return_path(analysis)),
	tar_file(model_file, write_lines_and_return_path(model)),
	tar_file(model_names_file, write_lines_and_return_path(names(model))),
	
	#Other data (e.g. lists, ggplot objects) as R data
	tar_file(histograms_file1_svg, write_svg_and_return_path(histograms[[1]])),
	tar_file(histograms_file2_svg, write_svg_and_return_path(histograms[[2]])),
	tar_file(histograms_file3_svg, write_svg_and_return_path(histograms[[3]])),
	
	tar_file(survplots_file1_svg, write_ggsurvplot_svg_and_return_path(survplots[[1]])),
	tar_file(survplots_file2_svg, write_ggsurvplot_svg_and_return_path(survplots[[2]])),
	
	tar_file(histograms_file1_png, write_png_and_return_path(histograms[[1]])),
	tar_file(histograms_file2_png, write_png_and_return_path(histograms[[2]])),
	tar_file(histograms_file3_png, write_png_and_return_path(histograms[[3]])),
	
	tar_file(survplots_file1_png, write_ggsurvplot_png_and_return_path(survplots[[1]])),
	tar_file(survplots_file2_png, write_ggsurvplot_png_and_return_path(survplots[[2]])),
	
	tar_file(cox_zph_file, save_and_return_path(cox_zph))
	
	
)
