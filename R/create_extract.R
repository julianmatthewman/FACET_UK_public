extract_cohort <- function(extract_files,
													 extract_codes,
													 study_start="1998-01-02",
													 study_end="2020-01-31",
													 practice_info_file,
													 patient_info_files) {
	
	#Get all rows with a code in the above defined
	extract <- open_dataset(extract_files) %>% 
		select(patid, eventdate, medcode) %>% 
		filter(medcode %in% extract_codes,
					 eventdate>=as.Date(study_start),
					 eventdate<=as.Date(study_end)) %>%
		collect()
	
	#Limit to first row per patient
	extract <- extract %>% 
		filter(!duplicated(patid)) %>% 
		mutate(pracid=as.numeric(str_sub(as.character(patid), start=-3)),
					 patid_short=as.numeric(str_sub(as.character(patid), end=-4)))
	
	#Get all patient information with a patid in the extract
	patient_info <- open_dataset(patient_info_files) %>% 
		select(patid, gender, realyob, crd, tod, deathdate) %>% 
		filter(patid %in% extract$patid) %>% 
		collect()
	
	#Get all practice information
	practice_info <- read_parquet(practice_info_file)
	
	#Join to the extract
	extract %>% 
		left_join(patient_info, by="patid") %>% 
		left_join(practice_info, by="pracid")
	
}
