testthat::test_that("propose_clock_times returns every requested example", {
  cases <- list(
    list(c("1h00","2h00","18h00"), c(1,2,18)),
    list(c("08:00","10:30","23:00","00:00"), c(8,10.5,23,0)),
    list(c("KSS_9u_dag1","KSS_21u_dag1","KSS_3u_dag2"), c(9,21,3)),
    list(c("11PM","12AM","1:30am","12PM"), c(23,0,1.5,12)),
    list(c("Base7h30","Base9h"), c(7.5,9)),
    list(c("T1","T2","T3"), NULL),
    list(c("VAS_10","VAS_12"), NULL))
  for(case in cases)testthat::expect_identical(propose_clock_times(case[[1]]),case[[2]],info=paste(case[[1]],collapse=", "))
})
testthat::test_that("clock parsing validates whole tokens and never partially guesses", {
  testthat::expect_identical(propose_clock_times(c("KSSS_8h30","VAS_9u30","PVT_1:30PM")),c(8.5,9.5,13.5))
  testthat::expect_identical(propose_clock_times(c("8h","9u"," 11 pm ")),c(8,9,23))
  for(bad in c("13PM","0AM","12:60am","25h","24h00","8h75","08:90","108:00","8h123","08:00:30","8h_10h","VAS_10","T1"))
    testthat::expect_null(propose_clock_times(bad),info=bad)
  testthat::expect_null(propose_clock_times(c("8h","unknown","10h")))
  testthat::expect_null(propose_clock_times(c("8h",NA_character_)))
  testthat::expect_null(propose_clock_times(character()))
})
testthat::test_that("midnight unwrapping advances all following times by 24 hours", {
  x<-c(seq(8,22,2),seq(0,14,2),15)
  testthat::expect_identical(unwrap_clock_times(x),c(seq(8,22,2),seq(24,38,2),39))
  testthat::expect_identical(unwrap_clock_times(c(22,0,2,0,1)),c(22,24,26,48,49))
  testthat::expect_identical(unwrap_clock_times(c(8,10,12)),c(8,10,12))
  testthat::expect_identical(unwrap_clock_times(c(8,32,56)),c(8,32,56))
  testthat::expect_identical(unwrap_clock_times(8),8)
  testthat::expect_error(unwrap_clock_times(c(8,8)),"equal.*recording order")
  testthat::expect_error(unwrap_clock_times(propose_clock_times(c("day1_8h","day2_8h"))),"equal")
  testthat::expect_error(unwrap_clock_times(c(22,0,0)),"equal")
  testthat::expect_error(unwrap_clock_times(c(48,0)),"recording order")
  testthat::expect_error(unwrap_clock_times(c(8,Inf)),"finite")
})
testthat::test_that("clock input counts and elapsed alignment are explicit", {
  testthat::expect_identical(parse_clock_times("8, 10.5, 12",3),c(8,10.5,12))
  testthat::expect_error(parse_clock_times("8, 10",3),"must match")
  testthat::expect_error(parse_clock_times("8, 10, 12",2),"must match")
  testthat::expect_error(parse_clock_times("",3),"No equally spaced")
  testthat::expect_error(parse_clock_times("8,,10",2),"empty value")
  testthat::expect_error(parse_clock_times("8, 10:30",2),"10.5")
  testthat::expect_error(parse_clock_times("8, NA",2),"missing")
  testthat::expect_error(parse_clock_times("-1,8",2),"negative")
  m<-clock_time_mapping("8, 22, 0, 2, 15",5,8)
  testthat::expect_identical(m$unwrapped,c(8,22,24,26,39))
  testthat::expect_identical(m$elapsed,c(0,14,16,18,31))
  testthat::expect_true(m$adjusted)
  testthat::expect_identical(clock_time_mapping("7.5,9",2,8)$elapsed,c(-.5,1))
})

testthat::test_that("numeric clock interpretation is opt-in and supports midnight", {
  testthat::expect_null(propose_clock_times(c("KSS_10","KSS_12")))
  testthat::expect_null(propose_clock_times(c("VAS_10","VAS_12")))
  for(prefix in c("KSS","KSSS","VAS","PVT","PVT_Lapses","Lapses")) {
    testthat::expect_identical(propose_numeric_clock_times(paste0(prefix,c("_22","_24","_2","_4_SD"))),c(22,0,2,4))
  }
  for(bad in list(c("T1","T2"),c("KSS_10","unknown"),"KSS_25","KSS_10.5","KSS_8h75",NA_character_))
    testthat::expect_null(propose_numeric_clock_times(bad))
  testthat::expect_identical(propose_numeric_clock_times(c("KSS_10","KSS_12h")),c(10,12))
  testthat::expect_error(unwrap_clock_times(propose_numeric_clock_times(c("KSS_8","KSS_8_SD"))),"equal")
})
