options(biomas.launch=FALSE,warn=1)
source("app.R")
biomas_test_root<-getwd()
dir.create("work/biomas",recursive=TRUE,showWarnings=FALSE)
notice_log<-character()
showNotification<-function(ui,...) {
  notice_log<<-c(notice_log,paste(as.character(ui),collapse=" "))
  shiny::showNotification(ui,...)
}
results<-testthat::test_dir("tests",reporter="summary",stop_on_failure=TRUE)
d<-as.data.frame(results)
cat(sprintf("\nTOTAL: %d tests; %d passed expectations; %d failures; %d errors; %d warnings.\n",
  nrow(d),sum(d$passed),sum(d$failed),sum(d$error),sum(d$warning)))
write.csv(d[,c("test","nb","failed","error","warning","passed")],"work/biomas/test_results.csv",row.names=FALSE)
