bed_setup<-function(session)session$setInputs(
  t0_clock=8,duration=48,dated=TRUE,start_iso="2026-10-20T08:00:00+02:00",input_mode="bedtime",auto_recalc=FALSE,
  model_version="2024fe",custom_params=FALSE,baseline_on=22,baseline_off=8,dt_choice="5",
  w1=.35,w2=.4,pvt_ceiling=20,index_threshold=.66,prediction_scale="raw",state_metric="PVT",state_variable="u",
  kss_decimal=".",vas_decimal=".",pvt_decimal=".",vas_min=0,vas_max=100,vas_direction="direct",
  kss_show_mean=TRUE,vas_show_mean=TRUE,pvt_show_mean=TRUE,pvt_lapse_rule=">=500 ms",
  bed_format="local",bed_date_order="dmy")

# Actiware-style rest-interval export: separate dates and times, REST plus
# sleep-scored SLEEP rows that lie inside the REST intervals.
bed_fixture<-function() {
  f<-tempfile(fileext=".csv")
  writeLines(c(
    "Interval Type,Start Date,Start Time,End Date,End Time",
    "REST,17/10/2026,22:45,18/10/2026,06:50",
    "SLEEP,17/10/2026,23:05,18/10/2026,06:30",
    "REST,18/10/2026,23:30,19/10/2026,07:10",
    "REST,19/10/2026,23:15,20/10/2026,07:00",
    "REST,20/10/2026,22:00,21/10/2026,08:00",
    "SLEEP,20/10/2026,22:20,21/10/2026,07:40",
    "REST,21/10/2026,22:30,22/10/2026,07:30",
    "REST,22/10/2026,22:00,23/10/2026,07:00"),f)
  f
}
bed_upload<-function(session,f)session$setInputs(bed_file=data.frame(name=basename(f),size=file.info(f)$size,type="text/csv",datapath=f))

testthat::test_that("local date-times use one explicit offset and an explicit date order", {
  ref<-fr_parse_iso("2026-10-18T23:15:00+02:00")
  testthat::expect_equal(fr_parse_local_datetime("18/10/2026 23:15","+02:00","dmy"),ref)
  testthat::expect_equal(fr_parse_local_datetime("10/18/2026 11:15 PM","+02:00","mdy"),ref)
  testthat::expect_equal(fr_parse_local_datetime("2026-10-18 23:15:00","+02:00","dmy"),ref)
  testthat::expect_equal(fr_parse_local_datetime("18.10.2026 23:15:00","+02:00","dmy"),ref)
  testthat::expect_error(fr_parse_local_datetime("18/10/2026 23:15","+02:00","ymd"),"date order")
  testthat::expect_error(fr_parse_local_datetime("31/02/2026 10:00","+02:00","dmy"),"invalid")
  testthat::expect_error(fr_parse_local_datetime("18/10/2026 13:15 PM","+02:00","dmy"),"12-hour")
  testthat::expect_error(fr_parse_local_datetime("2026-10-18 23:15",NA_character_),"UTC offset")
  testthat::expect_error(fr_parse_local_datetime("yesterday 23:15","+02:00"),"not a recognised")
  o<-fr_origin("dated",start_iso="2026-10-20T08:00:00+02:00")
  testthat::expect_equal(fr_bedtime_elapsed("20/10/2026 22:00","local",o,"dmy"),14)
  testthat::expect_equal(fr_bedtime_elapsed("2026-10-20T22:00:00Z","iso",o),16)
  testthat::expect_equal(fr_bedtime_elapsed(c("14","24.5"),"hours",o),c(14,24.5))
  testthat::expect_error(fr_bedtime_elapsed("20/10/2026 22:00","local",fr_origin(),"dmy"),"dated")
})

testthat::test_that("spreadsheet date and time cells combine without reinterpretation", {
  date<-as.Date(c("2026-10-20","2026-10-21"))
  time<-as.POSIXct(c("1899-12-31 22:00:00","1899-12-31 22:30:00"),tz="UTC")
  testthat::expect_identical(fr_join_date_time(date,time),c("2026-10-20 22:00:00","2026-10-21 22:30:00"))
  excel_date<-as.POSIXct(c("2026-10-20 00:00:00","2026-10-21 00:00:00"),tz="UTC")
  testthat::expect_identical(fr_join_date_time(excel_date,c("22:00","22:30")),c("2026-10-20 22:00","2026-10-21 22:30"))
  testthat::expect_identical(fr_cell_text(as.POSIXct("2026-10-20 22:00:00",tz="UTC")),"2026-10-20 22:00:00")
  testthat::expect_error(fr_join_date_time(date,c("2026-10-20 22:00","2026-10-21 22:30")),"already contain a date")
  testthat::expect_identical(fr_join_date_time(NULL,c("x",NA)),c("x",NA))
})

testthat::test_that("bedtime rows become whole sleep opportunities with an audited role", {
  r<-fr_bedtime_intervals(bed=c(-33,-9,20,-2,46,60),rise=c(-25,-3,30,1,52,70),end=48,t0_clock=8,rows=1:6)
  testthat::expect_identical(r$schedule,data.frame(onset_h=c(0,20,46),offset_h=c(1,30,48)))
  testthat::expect_identical(unname(r$counts),c(2L,3L,2L,1L))
  testthat::expect_identical(r$nights$row,c(1L,2L,4L,3L,5L,6L))
  testthat::expect_identical(r$nights$use,c("Pre-study (baseline)","Pre-study (baseline)","Study schedule (clipped)",
    "Study schedule","Study schedule (clipped)","After study end (ignored)"))
  testthat::expect_equal(r$nights$bed_clock[1:2],c(23,23))
  testthat::expect_equal(r$nights$time_in_bed_h,c(8,6,3,10,6,10))
  testthat::expect_error(fr_bedtime_intervals(c(1,5),c(6,9),48,8),"Overlapping")
  testthat::expect_error(fr_bedtime_intervals(5,4,48,8),"must follow")
  testthat::expect_error(fr_bedtime_intervals(5,30,48,8),"exceeds 24")
  testthat::expect_error(fr_bedtime_intervals(c(1,NA),c(2,3),48,8,rows=c(4,7)),"Row\\(s\\) 7")
  only_pre<-fr_bedtime_intervals(-10,-2,24,8)
  testthat::expect_identical(nrow(only_pre$schedule),0L)
})

testthat::test_that("habitual baseline uses circular means across midnight", {
  h<-fr_habitual_time_in_bed(c(23.5,.5),c(7,8))
  testthat::expect_equal(h$onset,0);testthat::expect_equal(h$offset,7.5);testthat::expect_equal(h$time_in_bed_h,7.5)
  h<-fr_habitual_time_in_bed(c(22.75,23.5,23.25),c(6+50/60,7+10/60,7))
  testthat::expect_equal(h$onset,23+10/60);testthat::expect_equal(h$offset,7)
  testthat::expect_equal(h$mean_night_tib_h,mean(c(8+5/60,7+40/60,7+45/60)))
  testthat::expect_identical(fr_hhmm(c(23+59.9/60,7.5,-1)),c("00:00","07:30","23:00"))
  testthat::expect_error(fr_habitual_time_in_bed(numeric(),numeric()),"No complete")
  testthat::expect_error(fr_habitual_time_in_bed(c(0,12),c(7,8)),"whole clock")
})

testthat::test_that("column proposals cover common actigraphy and diary exports", {
  g<-fr_guess_bedtime_columns(c("Interval Type","Start Date","Start Time","End Date","End Time"))
  testthat::expect_identical(unlist(g),c(bed="Start Time",bed_date="Start Date",rise="End Time",rise_date="End Date",type="Interval Type"))
  g<-fr_guess_bedtime_columns(c("ID","In Bed Date","In Bed Time","Out Bed Date","Out Bed Time"))
  testthat::expect_identical(c(g$bed,g$bed_date,g$rise,g$rise_date),c("In Bed Time","In Bed Date","Out Bed Time","Out Bed Date"))
  g<-fr_guess_bedtime_columns(c("night","lights_out","lights_on","sleep_onset"))
  testthat::expect_identical(c(g$bed,g$rise),c("lights_out","lights_on"))
  testthat::expect_true(is.na(fr_guess_bedtime_columns(c("sleep_onset","final_awakening"))$bed))
})

testthat::test_that("interface defaults to the paper's 10 h baseline time in bed", {
  html<-as.character(ui)
  testthat::expect_match(html,'id="baseline_on"[^>]*value="22"')
  testthat::expect_match(html,'id="baseline_off"[^>]*value="8"')
  testthat::expect_match(html,"Initial conditions: baseline time in bed",fixed=TRUE)
  testthat::expect_match(html,"not total sleep time",fixed=TRUE)
  testthat::expect_match(html,'value="bedtime"',fixed=TRUE)
  shiny::testServer(server,{
    bed_setup(session)
    testthat::expect_match(output$baseline_summary$html,"22:00-08:00: 10.00 h in bed, 14.00 h awake",fixed=TRUE)
    session$setInputs(baseline_on=23.5,baseline_off=7.5)
    testthat::expect_match(output$baseline_summary$html,"23:30-07:30: 8.00 h in bed",fixed=TRUE)
  })
})

testthat::test_that("bedtime-interval workflow simulates REST time in bed and derives the baseline", {
  shiny::testServer(server,{
    bed_setup(session)
    f<-bed_fixture();bed_upload(session,f)
    ui_html<-output$bed_columns_ui$html
    for(x in c('<option value="Start Time" selected>','<option value="Start Date" selected>',
               '<option value="End Time" selected>','<option value="End Date" selected>','<option value="Interval Type" selected>'))
      testthat::expect_match(ui_html,x,fixed=TRUE)
    session$setInputs(bed_col="Start Time",bed_date_col="Start Date",bed_rise_col="End Time",bed_rise_date_col="End Date",
      bed_type_col="Interval Type",bed_type_value="__all__")
    testthat::expect_error(bedtime_result(),"Overlapping")
    testthat::expect_match(output$bed_type_value_ui$html,'<option value="REST" selected>',fixed=TRUE)
    session$setInputs(bed_type_value="REST")
    r<-bedtime_result()
    testthat::expect_equal(active_schedule(),data.frame(onset_h=c(14,38.5),offset_h=c(24,47.5)))
    testthat::expect_identical(unname(r$counts),c(3L,2L,0L,1L))
    testthat::expect_identical(r$nights$row,c(1L,3L,4L,5L,7L,8L))
    testthat::expect_match(output$bed_status$html,"3 pre-study, 2 in the study (0 clipped), 1 after the study end",fixed=TRUE)
    testthat::expect_match(output$bed_status$html,"23:10-07:00 (7.8 h",fixed=TRUE)
    session$setInputs(bed_to_baseline=1)
    testthat::expect_true(any(grepl("Baseline time in bed set to 23:10-07:00 \\(7.8 h\\) from 3 pre-study",notice_log)))
    session$setInputs(baseline_on=23+10/60,baseline_off=7,btn_run=1)
    b<-valid_bundle()
    testthat::expect_equal(b$inputs$schedule,data.frame(onset_h=c(14,38.5),offset_h=c(24,47.5)))
    testthat::expect_equal(b$inputs$baseline,list(onset=23+10/60,offset=7))
    testthat::expect_identical(b$inputs$source$schedule$bedtime$settings$type_value,"REST")
    testthat::expect_identical(nrow(b$inputs$source$schedule$bedtime$nights),6L)
    testthat::expect_match(b$inputs$assumptions$baseline,"23:10-07:00",fixed=TRUE)
    pred<-b$results$predictions
    testthat::expect_false(pred$wake[match(14,pred$time_h)])
    testthat::expect_true(pred$wake[match(24,pred$time_h)])
    # The schedule run equals a direct core run of the same intervals and baseline.
    direct<-fr_run(b$inputs)
    testthat::expect_equal(direct$results$predictions,b$results$predictions)
    writeLines(fr_export_script(b),file.path(biomas_test_root,"work/biomas/reproduce_bedtime.R"))
    session$setInputs(bed_type_value="SLEEP")
    testthat::expect_identical(unname(bedtime_result()$counts),c(1L,1L,0L,0L))
    session$setInputs(bed_date_order="ymd")
    testthat::expect_error(bedtime_result(),"date order")
    session$setInputs(baseline_paper=1)
  })
})

testthat::test_that("pre-study baseline needs pre-study nights", {
  shiny::testServer(server,{
    bed_setup(session)
    f<-tempfile(fileext=".csv")
    writeLines(c("bedtime,risetime","2026-10-20 22:00,2026-10-21 07:00"),f)
    bed_upload(session,f)
    session$setInputs(bed_col="bedtime",bed_date_col="__none__",bed_rise_col="risetime",bed_rise_date_col="__none__",bed_type_col="__none__")
    testthat::expect_equal(active_schedule(),data.frame(onset_h=14,offset_h=23))
    n0<-length(notice_log)
    session$setInputs(bed_to_baseline=1)
    testthat::expect_true(any(grepl("No intervals end before t = 0",notice_log[-seq_len(n0)])))
  })
})
