#' Define measures for glucocorticoid use 
#' @details Define various measures for glucocorticoid prescriptions including if a prescription crosses a risk threshold of a predefined cumulative prednisolone equivalent dose (PED) using a 6 month rolling window
#' @param prescriptions_ped A dataframe containing patient IDs, event dates, and other information on corticosteroid prescriptions
#' @return The prescriptions dataset enhanced with additional variables on cumulative doses and the risk threshold.
create_prescriptions_risk_threshold <- function(prescriptions_ped) {
	library(zoo)
	library(tidyverse)
	
	prescriptions_ped %>% 
		group_by(patid) %>%
		mutate(
			gap=start-lag(end),
			gap=ifelse(gap<0, 0, gap), #Some prescriptions are prescribed before an active one ends; set the gap for these to 0
			rollmeangap=rollapply(gap,
														width=row_number() - findInterval(start - 180, start), #specify the number of prior observations that fall within the time window for each observation.
														FUN=mean, #apply the "mean" function to the rolling window
														align='right'), #(i.e.: the current day is at the right side of the window)
			peakdose=cummax(ped),
			ped_whole_duration=ped*duration,
			rollmeanped=rollapply(ped,
														width=row_number() - findInterval(start - 180, start),
														FUN=mean,
														align='right'),
			rollsumped=rollapply(ped_whole_duration,
													 width=row_number() - findInterval(start - 180, start),
													 FUN=sum,
													 align='right'),
			riskthreshold=ifelse(rollsumped>=450 & start >= indexdate, 1, 0),
			riskthreshold_window_stop=if(any(riskthreshold==1)) start[which.max(riskthreshold==1)] else NA,
			riskthreshold_window_start=if(any(riskthreshold==1)) start[which.max(start>=(riskthreshold_window_stop-180) & start<=riskthreshold_window_stop)] else NA,
			riskthreshold_window_length=as.numeric(riskthreshold_window_stop-riskthreshold_window_start),
			riskthreshold_sum_gaps=if(any(riskthreshold==1)) sum(gap[start>riskthreshold_window_start & start<=riskthreshold_window_stop]) else NA,
			riskthreshold_n_gaps=if(any(riskthreshold==1)) length(gap[gap>0 & start>riskthreshold_window_start & start<=riskthreshold_window_stop]) else NA,
			
		) %>% 
		ungroup()
}
