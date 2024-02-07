#' Create time merged cohort.
#' @description Make a start-stop dataset from the cohort and the various files of events. See: Therneau T, Crowson C, Atkinson E. Using Time Dependent Covariates and Time Dependent Coefficients in the Cox Model.

#' @param path_main_cohort Path to file.
#' @param path_prescriptions Path to file.
#' @param path_cprdfract Path to file.
#' @param path_hesfract Path to file.
#' @param path_alc Path to file.
#' @param path_asthma Path to file.
#' @param path_ethn Path to file.
#'
#' @return A dataframe containing start-stop observations.

create_cohort_tmerged <- function(main_cohort, 
																 prescriptions_risk_threshold,
																 prescriptions_expanded,
																 non_steroid_fx_drugs,
																 eventdata,
																 codelists) {

library(arrow)
library(rio)
library(haven)
library(survival)
library(lubridate)
library(tidyverse)
	
# Read data ----------------------------------------------------

#Import other data (and only keep IDs present in the main cohort)
fract_all <- eventdata[[which(codelists$name=="fractures_all")]]
fract_hip <- eventdata[[which(codelists$name=="fractures_hip")]]
fract_spine <- eventdata[[which(codelists$name=="fractures_spine")]]
fract_wrist <- eventdata[[which(codelists$name=="fractures_wrist")]]
fract_pelvis <- eventdata[[which(codelists$name=="fractures_pelvis")]]
asthma <- eventdata[[which(codelists$name=="asthma")]]
copd <- eventdata[[which(codelists$name=="copd")]]
eczema <- eventdata[[which(codelists$name=="eczema")]]
rheumatoid_arthritis <- eventdata[[which(codelists$name=="rheumatoid_arthritis")]]
bisphosphonate <- eventdata[[which(codelists$name=="bisphosphonates")]]
calcium_and_vit_d <- eventdata[[which(codelists$name=="calcium_and_vit_d")]]
dxa_scan <- eventdata[[which(codelists$name=="dxa_scan")]]
anxiety_meds <- eventdata[[which(codelists$name=="anxiety_meds")]]
epilepsy_meds <- eventdata[[which(codelists$name=="epilepsy_meds")]]
migraine_meds <- eventdata[[which(codelists$name=="migraine_meds")]]






# Make composite fracture outcome (hip, spine, wrist, and pelvis)
fract_composite <- rbind(fract_hip[c("patid", "eventdate")], 
												 fract_spine[c("patid", "eventdate")],
												 fract_wrist[c("patid", "eventdate")],
												 fract_pelvis[c("patid", "eventdate")])

# Make composite fracture preventive care outcome (bisphosphonates and DXA Scans)
fpc <- rbind(bisphosphonate[c("patid", "eventdate")], dxa_scan[c("patid", "eventdate")])

# Make composite fracture preventive care outcome (bisphosphonates and calcium and vitamin D)
bp_cal_vit_d <- rbind(bisphosphonate, calcium_and_vit_d)

#Make a table containing the timegroup cutpoints
timegroups <- tibble(
	patid=rep(unique(main_cohort$patid), each=2),
	eventdate=rep(date(c("2006-01-01", "2013-01-01")), length(unique(main_cohort$patid)))
)	 

#Make a table containing the agegroup cutpoints
agecuts <- c(40, 50, 66, 80)
agegroups <- main_cohort[rep(seq_len(nrow(main_cohort)), each = length(agecuts)), c("patid", "realyob")]
agegroups$eventdate <- agegroups$realyob + agecuts #The eventdate is the date at the agecut
agegroups$eventdate <- ymd(agegroups$eventdate, truncated = 2L)



#Make competing risks dataset (with both bisphsosphoantes and fractures) -------------------------

bp_single <- bisphosphonate %>% 
	group_by(patid) %>% 
	filter(!duplicated(eventdate)) %>% 
	select(patid, eventdate) %>% 
	mutate(event="bp")

fract_single <- fract_composite %>% 
	group_by(patid) %>% 
	filter(!duplicated(eventdate)) %>% 
	select(patid, eventdate) %>% 
	mutate(event="fract")

cr_fract_bp <- rbind(bp_single, fract_single) %>% 
	mutate(event=factor(event, levels = c("censor", "fract", "bp"))) %>% 
	arrange(patid, eventdate, event)

#Check distribution of event types
table(cr_fract_bp$event)

#Check if any events on same day
cr_fract_bp %>% 
	group_by(patid) %>% 
	filter(duplicated(eventdate))

#Move the event in the second factor level back one day if on the same day as other event type
cr_fract_bp <- cr_fract_bp %>% 
	group_by(patid) %>% 
	mutate(eventdate=if_else(duplicated(eventdate), eventdate+1, eventdate)) %>% 
	ungroup()



# Create a start-stop dataset using the tmerge function ------------------------

tmerged <- tmerge(main_cohort, main_cohort, id=patid, tstart = indexdate-(100*365), tstop = enddate) #Therneau: "The first call sets the time range for each subject to be from 0 (default) to last follow-up. If a later call tried to add an event outside that range, at time = -2 say, that addition would be ignored." In our case we want to capture everything before the end of the study, since we will set our follow up window later, so here we set tstart to 100 years before the indexdate.
tmerged <- tmerge(tmerged, prescriptions_risk_threshold, id=patid,  
									 #rollgap=tdc(start, rollgap),
									 #rollmeanped=tdc(start,rollmeanped),
									 rollsumped=tdc(start, rollsumped),
									 cumdose_cont=cumtdc(start, ped*duration, 0),
									 cumdays_cont=cumtdc(start, duration, 0),
									 riskthreshold=tdc(start, riskthreshold, 0))
tmerged <- tmerge(tmerged, prescriptions_expanded %>% mutate(ped=ifelse(is.na(ped), 0, ped)), id=patid, 
									active=tdc(start, active, "not active"),
									daily_ped=tdc(start, ped, 0))
tmerged <- tmerge(tmerged, eczema, id=patid, eczema=tdc(eventdate))
tmerged <- tmerge(tmerged, asthma, id=patid, asthma=tdc(eventdate))
tmerged <- tmerge(tmerged, copd, id=patid, copd=tdc(eventdate))
tmerged <- tmerge(tmerged, rheumatoid_arthritis, id=patid, rheumatoid_arthritis=tdc(eventdate))
tmerged <- tmerge(tmerged, non_steroid_fx_drugs, id=patid, non_steroid_fx_drugs=tdc(eventdate))
tmerged <- tmerge(tmerged, timegroups, id=patid, timegroup=cumtdc(eventdate))
tmerged <- tmerge(tmerged, agegroups, id=patid, agegroup=cumtdc(eventdate))
tmerged <- tmerge(tmerged, bisphosphonate, id=patid, bisphosphonate=event(eventdate))
tmerged <- tmerge(tmerged, calcium_and_vit_d, id=patid, calcium_and_vit_d=event(eventdate))
tmerged <- tmerge(tmerged, dxa_scan, id=patid, dxa_scan=event(eventdate))
tmerged <- tmerge(tmerged, bp_cal_vit_d, id=patid, bp_cal_vit_d=event(eventdate))
tmerged <- tmerge(tmerged, fpc, id=patid, fpc=event(eventdate))
tmerged <- tmerge(tmerged, fract_all, id=patid, fract_any=event(eventdate))
tmerged <- tmerge(tmerged, fract_hip, id=patid, fract_hip=event(eventdate))
tmerged <- tmerge(tmerged, fract_spine, id=patid, fract_spine=event(eventdate))
tmerged <- tmerge(tmerged, fract_wrist, id=patid, fract_wrist=event(eventdate))
tmerged <- tmerge(tmerged, fract_pelvis, id=patid, fract_pelvis=event(eventdate))
tmerged <- tmerge(tmerged, fract_composite, id=patid, fract_composite=event(eventdate))
tmerged <- tmerge(tmerged, cr_fract_bp, id=patid, cr_fract_bp=event(eventdate, event))
tmerged <- tmerge(tmerged, anxiety_meds, id=patid, anxiety_meds=event(eventdate))
tmerged <- tmerge(tmerged, epilepsy_meds, id=patid, epilepsy_meds=event(eventdate))
tmerged <- tmerge(tmerged, migraine_meds, id=patid, migraine_meds=event(eventdate))





tmerged <- left_join(tmerged, prescriptions_risk_threshold[!duplicated(prescriptions_risk_threshold$patid), c("patid", "riskthreshold_window_length", "riskthreshold_sum_gaps",  "riskthreshold_n_gaps")])


# Categorise variables and make factors ----------------------------------------------------

tmerged %>% 
	mutate(
		logdays=log10(riskthreshold_window_length + 1),
		cumdose=case_when(cumdose_cont == 0 ~ "0 g",
											cumdose_cont < 450 ~ "1 to 449 mg",
											cumdose_cont < 900 ~ "450 to 899 mg",
											cumdose_cont < 1500 ~ "900 to 1500 mg",
											cumdose_cont >= 1500 ~ ">1500 mg"),
		cumdose=factor(cumdose, levels = c("0 g", 
																			 "1 to 449 mg", 
																			 "450 to 899 mg", 
																			 "900 to 1500 mg", 
																			 ">1500 mg")),
		cumdays=case_when(cumdays_cont < 90 ~ "few",
											cumdays_cont >= 90 ~ "many"),
		cumdays=factor(cumdays, levels = c("few", "many")),
		dob=ymd(realyob, truncated = 2L), #Assume 1st January as Birthday
		age=as.numeric(tstart-dob)/365.25,
		sex=factor(sex),
		asthma=factor(asthma),
		copd=factor(copd),
		eczema=factor(eczema),
		rheumatoid_arthritis=factor(rheumatoid_arthritis),
		non_steroid_fx_drugs=factor(non_steroid_fx_drugs),
		carstairs=factor(carstairs),
		timegroup=factor(timegroup, levels = c(0, 1, 2),labels = c("1997-2005", "2006-20012", "2013-2020")),
		agegroup=factor(agegroup, levels = c(0, 1, 2, 3, 4), labels = c("18-39", "40-49", "50-65", "66-79", "80+")),
		
		#Exposures
		
		# Days taken to reach the risk threshold of 450mg PED
		pattern_rt_window=factor(ifelse(riskthreshold_window_length>=90, "90-180", ifelse(riskthreshold_window_length<90, "0-89", NA)), 
														 levels=c("90-180", "0-89")),
		pattern_rt_window0=factor(ifelse(riskthreshold_window_length>=1, "1-180", ifelse(riskthreshold_window_length==0, "0", NA)), 
															 levels=c("1-180", "0")),
		
		# Sum of gap days before reaching the risk threshold of 450mg PED
		pattern_rt_sum_gaps=factor(ifelse(riskthreshold_sum_gaps>=90, "90-180", ifelse(riskthreshold_sum_gaps<90, "0-89", NA)), 
															 levels=c("90-180", "0-89")),
		pattern_rt_sum_gaps0=factor(ifelse(riskthreshold_sum_gaps>=1, "1-180", ifelse(riskthreshold_sum_gaps==0, "0", NA)), 
															 levels=c("1-180", "0")),
		
		
		# Number of gaps before reaching the risk threshold of 450mg PED
		pattern_rt_n_gaps=factor(ifelse(riskthreshold_n_gaps>=2, "2+", ifelse(riskthreshold_n_gaps<2, "0-1", NA)), 
														 levels=c("2+", "0-1")),
		pattern_rt_n_gaps0=factor(ifelse(riskthreshold_n_gaps>=1, "1+", ifelse(riskthreshold_n_gaps==0, "0", NA)), 
														 levels=c("1+", "0")),
		
		
		
		)


}
