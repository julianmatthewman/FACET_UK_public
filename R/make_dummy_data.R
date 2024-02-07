#' Make dummy data
#' @description Load dummy datasets and path names, and writes them to files

make_dummy_data <- function() {
	library(dplyr)

# Specifications ----------------------------------------------------------

n <- 10000 #Set number of required participants
prop_presc <- 0.5 #Set number of participants with at least one prescription


# dummy main_cohort -------------------------------------------------------
main_cohort <- data.frame(patid=seq(n))
main_cohort$sex <- rbinom(n, 1, 0.5)
main_cohort$dob <- sample(seq(as.Date('1930/01/01'), as.Date('1990/01/01'), by="day"), n, replace = TRUE)
main_cohort$realyob <- year(main_cohort$dob)
main_cohort$indexdate <- sample(seq(as.Date('1998/01/01'), as.Date('2020/01/01'), by="day"), n, replace = TRUE)
main_cohort$enddate <- main_cohort$indexdate + abs(rnorm(n, mean=365*6.5, sd=365*2))
main_cohort$bmi <- rnorm(n, mean = 27, sd = 6)
main_cohort$bmi_cat <- rbinom(n, 3, 0.3)
main_cohort$smokstatus <- rbinom(n, 1, 0.5)
main_cohort$carstairs <- sample(seq(5), n, replace = TRUE)
main_cohort$eth5 <- rbinom(n, 4, 0.1)
main_cohort$eth16 <- rbinom(n, 4, 0.1)
main_cohort$exposed <- 1
main_cohort$pracid <- sample(seq(500), n, replace = TRUE)

main_cohort <<- main_cohort
path_main_cohort <<- "dummy_data/dummy_main_cohort.csv"
write.csv(main_cohort, path_main_cohort, row.names = FALSE)



# dummy prescriptions -----------------------------------------------------
#Specify a proportion of the original patids
prescriptions_ped <- data.frame(patid=sample(sample(main_cohort$patid, n*prop_presc), (n*prop_presc)*30, replace = TRUE))
#hist(as.integer(table(prescriptions_ped$patid))) #How many people have one prescription, two prescriptions, three prescriptions, etc...
prescriptions_ped$ped <- rgamma(10000, shape = 13)
prescriptions_ped$high <- ifelse(prescriptions_ped$ped>=20, 1, 0)
prescriptions_ped$duration <- ceiling(abs(rnorm(10000, mean=22, sd=30)))+1
prescriptions_ped$start <- sample(seq(as.Date('1998/01/01'), as.Date('2020/01/01'), by="day"), n, replace = TRUE)
prescriptions_ped$end <- prescriptions_ped$start + prescriptions_ped$duration

prescriptions_ped <- arrange(prescriptions_ped, patid, start) %>% left_join(main_cohort[c("patid", "indexdate")], by="patid")

prescriptions_ped <<- prescriptions_ped
path_prescriptions_ped <<- "dummy_data/dummy_prescriptions_ped.csv"
write.csv(prescriptions_ped, path_prescriptions_ped, row.names = FALSE)



# dummy severity ---------------------------------------------------------------
severity <- data.frame(patid=sample(main_cohort$patid, n*0.16))
severity$date <- sample(seq(as.Date('1998/01/01'), as.Date('2020/01/01'), by="day"), n*0.16, replace = TRUE)
severity$modsevere <- rbinom(n*0.16, 2, 0.23)+1

severity <- arrange(severity, date)

severity <<- severity
path_severity <<- "dummy_data/dummy_severity.csv"
write.csv(severity, path_severity, row.names = FALSE)



# Dummy time updated variables -------------------------------------------------

# Specify name and prevalence
vars <- c(
	"fract_all" = 0.01,
	"fract_hip" = 0.008,
	"fract_spine" = 0.006,
	"fract_wrist" = 0.004,
	"fract_pelvis" = 0.002,
	"bisphosphonate" = 0.2,
	"calcium_and_vit_d" = 0.03,
	"alc" = 0.001,
	"antabuse" = 0.001,
	"asthma" = 0.01,
	"copd" = 0.001,
	"rheumatoid_arthritis" = 0.005
)

for (i in seq_along(vars)) {
	# Sample from main cohort with given prevalence
	data <- data.frame(patid=sample(main_cohort$patid, n*vars[[i]]))
	
	# Assign random date to events
	data$eventdate <- sample(seq(as.Date('1998/01/01'), as.Date('2020/01/01'), by="day"), n*vars[[i]], replace = TRUE)
	data <- arrange(data, eventdate)
	
	# Assign object to global enivronment
	assign(names(vars[i]), data, envir = .GlobalEnv)
	
	# Make path, assign to global environment and write to file
	path <- paste0("dummy_data/", names(vars[i]), ".csv")
	assign(paste0("path_", names(vars[i])), path, envir = .GlobalEnv)
	write.csv(data, path, row.names = FALSE)
}


# dummy carstairs ---------------------------------------------------------

patient_carstairs <- tibble(
	patid=1:n,
	pracid=sample(1:1000, n, replace = TRUE),
	carstairs2011_5=sample(1:5, n, replace = TRUE)
)

patient_carstairs <<- patient_carstairs
path_patient_carstairs <<- "dummy_data/patient_carstairs.csv"
write.csv(patient_carstairs, path_patient_carstairs, row.names = FALSE)


practice_carstairs <- tibble(
	pracid=sample(1:1000, n, replace = TRUE),
	country=sample(c("England", "Ireland", "Scotland", "Wales"), n, replace = TRUE),
	carstairs2011_5=sample(1:5, n, replace = TRUE)
)

practice_carstairs <<- practice_carstairs
path_practice_carstairs <<- "dummy_data/practice_carstairs.csv"
write.csv(practice_carstairs, path_practice_carstairs, row.names = FALSE)

}
