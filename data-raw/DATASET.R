## Set working directory to fdpbandsdata/
## Creates the file "R/sysdata.rda", which will contain all generated tables

stband_table_dir  <- "data-raw/quantile-tables/"
uniband_table_dir <- "data-raw/u_gamma-tables/"
save_to           <- "R/sysdata.rda"

for (file in list.files(stband_table_dir)) {
  if (file != "README.txt") {
    load(paste0(stband_table_dir, file)); assign(file, qtable)
  }
}

for (file in list.files(uniband_table_dir)) {
  if (file != "README.txt") {
    load(paste0(uniband_table_dir, file)); assign(file, utable)
  }
}

dont_save <- c(
               "file",
               "qtable",
               "save_to",
               "stband_table_dir",
               "uniband_table_dir",
               "utable",
               "dont_save"
             )
save_this <- ls()
save_this <- save_this[!save_this %in% dont_save]

save(list = save_this, file = save_to, envir = .GlobalEnv)

