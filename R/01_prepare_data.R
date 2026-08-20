source("R/00_functions.R")

input_csv <- "data/df_propen_matched.csv"
dta <- read.csv(input_csv)

trim_endpoint <- function(time, status, cap = CAP_MONTHS) {
  list(time = pmin(time, cap), status = ifelse(time > cap, 0L, as.integer(status)))
}
os <- trim_endpoint(dta$Time_OS, dta$Status)
pfs <- trim_endpoint(dta$Time_PFS, dta$Recurrence_progression)
dta$Time_OS <- os$time
dta$Status <- os$status
dta$Time_PFS <- pfs$time
dta$Recurrence_progression <- pfs$status

stopifnot(max(dta$Time_OS) <= CAP_MONTHS, max(dta$Time_PFS) <= CAP_MONTHS)
saveRDS(dta, "data/dta_matched_72.rds")
