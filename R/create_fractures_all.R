#' Join CPRD and HES fractures.
#'
#' @param path_cprdfract 
#' @param path_hesfract 
#' @param main_cohort 
#'
#' @return A dataset including both CPRD and HES fractures.
#' @export
#'
#' @examples
create_fractures_all <- function(path_cprdfract, path_hesfract, main_cohort) {
	#Import fractures data (and only keep IDs present in the main cohort)
	cprdfract <- read_dta(path_cprdfract) %>%
		filter(patid %in% main_cohort$patid)
	hesfract <- read_dta(path_hesfract) %>%
		filter(patid %in% main_cohort$patid) %>%
		rename(eventdate=epistart, #Rename variables to match CPRD fractures
					 prox_hum=Proximalhumerus,
					 Hip=Hipfractures,
					 spine=Spine,
					 wrist=Wrist)
	#Join
	allfract <- full_join(cprdfract, hesfract) %>% mutate(include_either=ifelse(include==1 | All_fractures==1, 1, 0))
}