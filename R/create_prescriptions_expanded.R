#' Expand prescriptions dataset to include periods of non-use
#' @param prescriptions_ped A dataframe containing patient IDs, event dates, and other information on corticosteroid prescriptions
#' @return The prescriptions dataset expanded by periods of non-use
#'
create_prescriptions_expanded <- function(prescriptions_ped) {
	library(dtplyr)
	# Expand prescriptions data to include periods of non-use ---------------------------------------------------------------------
	
	#Expand the prescriptions dataset
	presc_exp <- prescriptions_ped[rep(1:nrow(prescriptions_ped), each = 2), ] #Duplicate each row
	presc_exp[1:nrow(presc_exp) %% 2 == 0, ] <- NA #Set every second row to NA
	presc_exp <- fill(presc_exp, indexdate, patid, .direction = "down") #Copy values down to empty cells
	presc_exp$active <- rep(c("active", "not active"), nrow(prescriptions_ped)) #Make every first row "active", and every second row "not active"
	
	
	#Make starts and ends for the non active periods
	lazy_dt(presc_exp) %>% 
		group_by(patid) %>%
		mutate(end=case_when(!is.na(end) ~ as_date(end), #If its already there
												 is.na(lead(start)) ~ as_date(100000), #If its a participants last observation
												 is.na(end) ~ as_date(lead(start))), #Else, make it the start date of the next observation
					 start=case_when(!is.na(start) ~ as_date(start),
					 								is.na(start) ~ as_date(lag(end)))) %>%
		filter(start<end & end>=indexdate) %>%  #Keep only observations where the start is before the end and that end within or after the follow up period
		as_tibble()
	
	
}