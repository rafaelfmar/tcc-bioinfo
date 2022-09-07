Log <- function(message) {
  current_timestamp <- paste("[", format(Sys.time(), '%Y-%m-%d %H:%M:%S'), "]", sep = "")
  log_file_location <- "data/log.txt"
  
  cat("\n", current_timestamp, " ", message, sep = "")
}
