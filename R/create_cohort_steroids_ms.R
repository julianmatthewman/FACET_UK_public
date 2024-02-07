#' Make cohort of people with eczema with high cumulative steroid use.
#' @description For people with eczema with high cumulative corticosteroid use make survival datasets for different outcomes
#' @param cohort_tmerged A start-stop dataset with time updated variables for all people with eczema.
#' @param outcome A character vector containing names of outcome variables.
#' @param analysis A character vector containing names of analyses.
#' @return A list containing the survival datasets, the outcome name, participant counts for each step of restricting the dataset and tables for the distribution of exposures
create_cohort_steroids_ms <-  function(cohort_tmerged, analysis, eventdata_cancer, outcome) {
	library(survival)
	library(lubridate)
	library(tidyverse)
	library(magrittr)
	library(collapse)

	temp0 <- cohort_tmerged
	
	
	
	# Remove those that never cross the risk threshold --------------------------------------------------------
	
	# temp1 <- temp0 %>%
	# 	group_by(patid) %>% # Group by patient ID ...
	# 	filter(any(riskthreshold == 1 & tstart>=indexdate)) %>% # ... keep only those that cross the risk threshold after joining the eczema cohort ...
	# 	mutate(startfudate=min(tstart[riskthreshold==1 & tstart>=indexdate]), enddate=max(tstop)) %>% # ...set startfudate (for start follow up) to be the date they crossed the risk threshold.
	# 	ungroup()
	
	
	include <- temp0 %>% 
		fsubset(riskthreshold==1 & tstart>=indexdate) %>%
		fsubset(!duplicated(patid), patid)
	
	temp1 <- temp0 %>% 
		fsubset(patid %in% include$patid) %>% 
		fgroup_by(patid) %>% 
		fmutate(startfudate=min(tstart[riskthreshold==1 & tstart>=indexdate]),
						enddate=max(tstop)) %>% 
		fungroup() %>% 
		fmutate(startfudate=as_date(startfudate),
						riskthreshold_crossed=ifelse(tstart>=startfudate, 1, 0))
	
	
	
	
	# SENSITIVITY ANALYSIS: restrict to those aged 66+ when they cross the risk threshold for the first time
	if (analysis!="sens_all_ages") {
		# temp1 <- temp1 %>% 
		# 	group_by(patid) %>% 
		# 	filter(min(age[riskthreshold==1 & tstart>=indexdate])>=66) %>% 
		# 	ungroup()
		
		include_66 <- temp1 %>% 
			fgroup_by(patid) %>% 
			fsummarise(include=min(age[riskthreshold==1 & tstart>=indexdate])>=66) %>% 
			subset(include==TRUE, patid)
		
		temp2 <- temp1 %>% 
			fsubset(patid %in% include_66$patid)
	} else {
		temp2 <- temp1
	}
	
	
	# Exclude people with fractures or bisphosphonate prescriptions before the index date --------------------------------------------------------------
	
	# temp2 <- temp1 %>% 
	# 	group_by(patid) %>% 
	# 	mutate(riskthreshold_crossed=ifelse(tstart>=startfudate, 1, 0)) %>% 
	# 	filter(!any((bisphosphonate == 1 & riskthreshold_crossed==0) | (fract_composite == 1 & riskthreshold_crossed==0))) %>% 
	# 	ungroup()
	
	exclude_prev_event <- temp2 %>% 
		fgroup_by(patid) %>% 
		fsummarise(exclude=any((bisphosphonate == 1 & riskthreshold_crossed==0) | (fract_composite == 1 & riskthreshold_crossed==0))) %>% 
		fsubset(exclude==TRUE, patid)
	
	temp3 <- temp2 %>% 
		fsubset(!(patid %in% exclude_prev_event$patid))
	
	
	
	# Exclude people with cancer in the 6 month preceding index date --------------------------------------------------------------
	
	exclude_cancer <- temp3 %>% 
		fgroup_by(patid) %>% 
		fsummarise(indexdate=min(indexdate)) %>% 
		left_join(eventdata_cancer) %>% 
		filter(eventdate>=(indexdate-180), eventdate<=indexdate)
	
	temp4 <- temp3 %>% 
		fsubset(!(patid %in% exclude_cancer$patid)) %>% 
		fmutate(startfudate0=startfudate)
	
	
	
	# Make rolling windows ----------------------------------------------------------
	
	if (analysis!="sens_shorter_windows") {
		cutpoints <- c(-90,0,90*1:400) #Split every 90 days
	} else { # SENSITIVITY ANALYSIS: restrict to those aged 66+ when they cross the risk threshold for the first time
		cutpoints <- c(-45,0,45*1:800) #Split every 45 days
	}
	
	#Then split into windows of 3 months (the survival object will be attached to the split data)
	temp5 <- survSplit(Surv(time = as.numeric(tstart), 
													time2 = as.numeric(tstop),
													origin = as.numeric(startfudate),
													event = fract_composite) ~ ., 
										 data=temp4, 
										 cut=cutpoints,
										 episode="rolling_window",
										 event = "fract_composite")
	
	
	
	# Create multistate event variable ----------------------------------------
	temp6 <- temp5 %>% 
		mutate(ever_bisphosphonate=ifelse(bisphosphonate==0, NA, bisphosphonate)) %>% 
		group_by(patid) %>% 
		fill(ever_bisphosphonate) %>% 
		ungroup() %>% 
		mutate(ever_bisphosphonate=ifelse(is.na(ever_bisphosphonate), 0, ever_bisphosphonate),
					 event=ifelse(fract_composite==1, 2, ever_bisphosphonate),
					 event=factor(event, levels=0:2, labels = c("censor", "fpc", "fract")))
	
	
	# Stop follow up at fracture but not at bp ----------------------------------
	temp7 <- temp6 %>%
		group_by(patid) %>%
		slice(if(any(event == "fract" & riskthreshold_crossed==1)) 1:which.max(event == "fract"  & riskthreshold_crossed==1) else row_number()) %>% 
		mutate(event_2=as.character(event),
					 event_2=ifelse(event_2=="fract" & lag(event_2)=="fpc", "fract_after_fpc", event_2),
					 event_2=factor(event_2, levels=c("censor", "fpc", "fract_after_fpc", "fract")),
					 enddate=max(tstop)) %>% 
		ungroup()
	
	
	#Make survival object
	temp8 <- temp7
	
	temp8$outcome_surv <- Surv(time = as.numeric(temp7$tstart), 
											 time2 = as.numeric(temp7$tstop),
											 origin = 0,
											 event = temp7[[outcome]])

	
	
	
	# Limit follow up time (skip for sens analysis) ---------------------------
	if (analysis!="sens_all_follow_up") {
		#Limit follow up time to 6 90 day windows (2 windows are before the indexdate)
		temp8 <- temp8 %>% 
			filter(rolling_window<7)
		
		# Set enddate to be the stop of the remaining observations
		temp8 <- temp8 %>% 
			group_by(patid) %>% 
			mutate(enddate=max(outcome_surv[,"stop"])) %>% 
			ungroup()
	}
	
	# Exclude people who ever have a record for asthma ---------------------------
	if (analysis=="sens_no_asthma") {
		temp8 <- temp8 %>% 
			group_by(patid) %>% 
			filter(!any(asthma==1))
	}
	
	
	# Make the exposure variable (prescription pattern) -----------------------------------
	pattern <- temp8 %>% 
		mutate(time=outcome_surv[,"stop"]-outcome_surv[,"start"]) %>%  #Calculate time for each observation
		select(patid, rolling_window, riskthreshold_crossed, time, active, outcome_surv)
	
	pattern <- pattern %>% 
		group_by(patid, rolling_window) %>% 
		summarise(active_prop=sum(time[active=="active"])/sum(time)) %>% #Calculate the proportion of active time as the sum of events where there is an active prescription divided by the total time of the window
		mutate(pattern=factor(ifelse(active_prop>=.66, "continuous", "intermittent"), levels=c("intermittent", "continuous")), #Set pattern to "continuous" if active time is more than the defined threshold
					 pattern_lagged=lag(pattern), #SENSITIVITY ANALYSIS: Set pattern to "continuous" if active time in the previous window is more than the defined threshold
					 pattern_itt=nth(pattern, 3)) #SENSITIVITY ANALYSIS: Set pattern to "continuous" if active time in the second window is more than the defined threshold
	
	
	data <- temp8 %>% left_join(pattern, by=c("patid", "rolling_window"))
	
	#Get number of 90 day windows with continuous and intermittent use
	pattern_table <- map(pattern[c("pattern", "pattern_lagged", "pattern_itt")], table)
	
	
	# Start follow up on startfudate --------------------------------------------
	
	data <- data %>%
		group_by(patid) %>%
		slice(which.max(outcome_surv[,"start"]==0):n())
	
	
	# Drop missing levels (e.g. for cumulative dose variable there are no longer people without steroid use) -----------------
	
	data$cumdose <- droplevels(data$cumdose)
	
	
	# Get counts and person years ---------------------------------------------
	
	# Function to get person years
	getpyrs <- function(x, start, end) {
		x[!duplicated(x$patid),] %>% 
			mutate(futime=as.numeric(({{ end }} - {{ start }})/365.25)) %>% 
			pull(futime) %>% 
			sum()
	}
	
	# Get participant counts and person years for each step
	counts <- tibble(
		n = c(length(unique(temp0$patid)),
					length(unique(temp1$patid)),
					length(unique(temp2$patid)),
					length(unique(temp3$patid)),
					length(unique(temp4$patid)),
					length(unique(temp5$patid)),
					length(unique(data$patid))),
		step = c("Eczema cohort",
						 "Ever cross risk-threshold",
						 "Excluded individuals aged <66",
						 "Excluded individuals with previous major osteoporotic fracture or bisphosphonate prescription",
						 "Excluded individuals with cancer in the previous 6 months",
						 "Restricted to first occurence of outcome",
						 "From crossing the risk threshold to a maximum of one year thereafter"))
	
	
	# Create an object containing checks for the survival dataset ----------------------------------------------
	
	sc <- survcheck(data$outcome_surv ~ 1, data=data, id=patid, timefix=FALSE)
	
	
	
	
	# This will be saved to the list: -----------------------------------------
	set_names(list(data, outcome, analysis, counts, pattern_table, sc),
						c("data", "outcome", "analysis", "counts", "pattern_table", "survcheck"))
	
}

# # View structure of resulting list ----------------------------------------
# str(cohort_steroids, max.level=2, give.attr = FALSE)