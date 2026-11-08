setup_inputs<-function(session)session$setInputs(
  t0_clock=8,duration=36,dated=FALSE,start_iso="",input_mode="manual",auto_recalc=FALSE,
  model_version="2024fe",custom_params=FALSE,baseline_on=23,baseline_off=7,dt_choice="1",
  w1=.35,w2=.4,pvt_ceiling=20,index_threshold=.66,prediction_scale="raw",state_metric="PVT",state_variable="u",
  kss_decimal=".",vas_decimal=".",pvt_decimal=".",vas_min=0,vas_max=100,vas_direction="direct",
  kss_show_mean=TRUE,vas_show_mean=TRUE,pvt_show_mean=TRUE,pvt_lapse_rule=">=500 ms")
upload<-function(session,prefix,columns,selection=columns,values=NULL) {
  f<-tempfile(fileext=".csv")
  if(is.null(values))values<-if(prefix=="vas")rep(25,length(columns))else rep(4,length(columns))
  writeLines(c(paste(c("ID",columns),collapse=","),paste(c("s1",values),collapse=",")),f)
  args<-setNames(list(data.frame(name=basename(f),size=file.info(f)$size,type="text/csv",datapath=f),
    selection,"ID","__none__"),paste0(prefix,c("_file","_cols_selected","_id_col","_group_col")))
  do.call(session$setInputs,args)
  invisible(f)
}
testthat::test_that("all three inputs propose clocks in file order and align predictions", {
  shiny::testServer(server,{
    setup_inputs(session)
    for(prefix in c("kss","vas","pvt")){
      columns<-paste0(toupper(prefix),c("_22h","_0h","_2h"))
      upload(session,prefix,columns,rev(columns))
      r<-imported[[toupper(prefix)]]()
      testthat::expect_identical(r$settings$columns,columns)
      testthat::expect_identical(r$settings$clock_hours,c(22,0,2))
      testthat::expect_identical(r$settings$recording_hours,c(22,24,26))
      testthat::expect_identical(r$settings$times,c(14,16,18))
      testthat::expect_match(output[[paste0(prefix,"_time_note")]]$html,"text-success")
      testthat::expect_match(output[[paste0(prefix,"_time_note")]]$html,"Check and edit before running",fixed=TRUE)
      testthat::expect_match(output[[paste0(prefix,"_ui_times")]]$html,'value="22, 0, 2"',fixed=TRUE)
    }
    testthat::expect_true(any(grepl("Midnight adjustment.*22, 24, 26",notice_log)))
    session$setInputs(btn_run=1)
    b<-valid_bundle()
    for(m in c("KSS","VAS","PVT")){
      testthat::expect_identical(b$results$comparisons[[m]]$cells$time_h,c(14,16,18))
      testthat::expect_true(all(c(14,16,18)%in%b$results$predictions$time_h))
    }
    testthat::expect_false(is.null(output$plt_kss_pvt))
    testthat::expect_false(is.null(output$kss_plt_ts))
    testthat::expect_false(is.null(output$vas_plt_ts))
    testthat::expect_false(is.null(output$pvtlap_plt_ts))
    session$setInputs(kss_times_txt="21.5, 0.5, 2.5",t0_clock=9)
    testthat::expect_identical(imported$KSS()$settings$clock_text,"21.5, 0.5, 2.5")
    testthat::expect_identical(imported$KSS()$settings$times,c(12.5,15.5,17.5))
    session$setInputs(kss_btn_redetect=1)
    testthat::expect_identical(imported$KSS()$settings$clock_text,"21.5, 0.5, 2.5")
    session$setInputs(kss_cols_selected=c("KSS_22h","KSS_0h"))
    testthat::expect_error(imported$KSS(),"must match")
    session$setInputs(kss_cols_selected=c("KSS_22h","KSS_0h","KSS_2h"))
    testthat::expect_identical(imported$KSS()$settings$clock_text,"21.5, 0.5, 2.5")
  })
})
testthat::test_that("ambiguous labels stay empty and never invent spacing", {
  shiny::testServer(server,{
    setup_inputs(session)
    for(prefix in c("kss","vas","pvt")){
      upload(session,prefix,paste0(toupper(prefix),c("_10","_12")))
      testthat::expect_error(imported[[toupper(prefix)]](),"No equally spaced")
      testthat::expect_match(output[[paste0(prefix,"_ui_times")]]$html,'value=""',fixed=TRUE)
      testthat::expect_match(output[[paste0(prefix,"_time_note")]]$html,"No explicit clock times could be read",fixed=TRUE)
    }
    testthat::expect_false(sim_inputs_live()$ok)
    session$setInputs(kss_times_txt="10, 12",vas_times_txt="10, 12",pvt_times_txt="10, 12",btn_run=1)
    testthat::expect_false(is.null(valid_bundle()))
    session$setInputs(kss_times_txt="10")
    testthat::expect_error(imported$KSS(),"must match")
    testthat::expect_error(valid_bundle())
    session$setInputs(kss_times_txt="10, 10")
    testthat::expect_error(imported$KSS(),"equal")
    session$setInputs(kss_times_txt="10, 34",btn_run=2)
    testthat::expect_identical(valid_bundle()$inputs$empirical$KSS$settings$times,c(2,26))
  })
})
testthat::test_that("restored entries survive initial loading and origin validation recovers", {
  shiny::testServer(server,{
    setup_inputs(session)
    session$setInputs(kss_times_txt="9, 11")
    upload(session,"kss",c("KSS_8h","KSS_10h"))
    testthat::expect_identical(imported$KSS()$settings$clock_text,"9, 11")
    session$setInputs(btn_run=1)
    testthat::expect_identical(valid_bundle()$inputs$empirical$KSS$settings$times,c(1,3))
    testthat::expect_warning(session$setInputs(dated=TRUE,start_iso=""),NA)
    testthat::expect_false(sim_inputs_live()$ok)
    testthat::expect_error(valid_bundle())
    testthat::expect_warning(session$setInputs(kss_btn_redetect=1),NA)
    testthat::expect_warning(session$setInputs(start_iso="2026-10-20T08:00:00"),NA)
    testthat::expect_false(sim_inputs_live()$ok)
    session$setInputs(start_iso="2026-10-20T09:00:00+02:00",btn_run=2)
    testthat::expect_identical(imported$KSS()$settings$clock_text,"9, 11")
    testthat::expect_identical(valid_bundle()$inputs$empirical$KSS$settings$times,c(0,2))
    upload(session,"kss",c("KSS_12h","KSS_14h"))
    testthat::expect_identical(imported$KSS()$settings$clock_hours,c(12,14))
    upload(session,"kss",c("T1","T2"))
    testthat::expect_error(imported$KSS(),"No equally spaced")
    session$setInputs(kss_times_txt="12, 14",auto_recalc=TRUE)
    session$elapse(1000)
    testthat::expect_false(is.null(valid_bundle()))
    bundle<-valid_bundle()
    saveRDS(bundle,file.path(biomas_test_root,"work/biomas/test_bundle.rds"))
    writeLines(fr_export_script(bundle),file.path(biomas_test_root,"work/biomas/reproduce.R"))
    testthat::expect_match(fr_export_script(bundle),"BioMas reproducible analysis",fixed=TRUE)
    testthat::expect_match(fr_export_script(bundle),"propose_clock_times <-",fixed=TRUE)
  })
})
testthat::test_that("branding and original layout are retained", {
  html<-as.character(ui)
  testthat::expect_match(html,"BioMas 1.3.5",fixed=TRUE)
  testthat::expect_match(html,"Biomathematical Sleep Modeling",fixed=TRUE)
  for(tab in c("schedule","params","results","compare","export"))testthat::expect_match(html,paste0('id="shiny-tab-',tab,'"'),fixed=TRUE)
  testthat::expect_false(grepl("FaRiMa|FaRiMaTo|Provenance",html))
  testthat::expect_match(html,"NOT VALIDATED",fixed=TRUE)
})

testthat::test_that("the supplied file's calculation and calibration functions are preserved", {
  exprs<-parse("fixtures/supplied_app.R",keep.source=FALSE)
  kept<-character()
  for(e in exprs) {
    if(!is.call(e)||!identical(e[[1]],as.name("<-"))||!is.symbol(e[[2]])||
       !is.call(e[[3]])||!identical(e[[3]][[1]],as.name("function")))next
    nm<-as.character(e[[2]])
    if(!(grepl("^fr_",nm)||nm%in%c("mccauley_params","make_calibration_suggester")))next
    if(nm%in%c("fr_time_token","fr_propose_times","fr_bulk_times","fr_export_script"))next
    expected<-eval(e[[3]])
    actual<-get(nm,envir=globalenv())
    testthat::expect_identical(body(actual),body(expected),info=nm)
    testthat::expect_identical(formals(actual),formals(expected),info=nm)
    kept<-c(kept,nm)
  }
  testthat::expect_true(all(c("fr_ccc","fr_inter_ccc_gee","fr_suggest_calibration","fr_compare","fr_run")%in%kept))
  testthat::expect_true(all(c("fr_ccc","%||%","clock_time_mapping")%in%names(fr_code_snapshot())))
})

testthat::test_that("VAS scaling controls and added comparison UI survive clock changes", {
  shiny::testServer(server,{
    setup_inputs(session)
    session$setInputs(vas_target="01")
    upload(session,"vas",c("VAS_22h","VAS_0h","VAS_2h"),values=c(0,50,100))
    r<-imported$VAS()
    testthat::expect_identical(r$settings$times,c(14,16,18))
    testthat::expect_identical(r$settings$target,"01")
    testthat::expect_equal(r$observations$value,c(0,.5,1))
    testthat::expect_match(output$vas_range_note,"0 to 100")
    session$setInputs(vas_min=0,vas_max=200,vas_target="kss")
    testthat::expect_equal(imported$VAS()$observations$value,c(1,3,5))
    session$setInputs(vas_times_txt="21, 0, 3")
    testthat::expect_identical(imported$VAS()$settings$times,c(13,16,19))
    testthat::expect_identical(imported$VAS()$settings$vas_max,200)
  })
  html<-as.character(ui)
  for(label in c("Suggest calibration","vas_target","kss_ui_subjects","state_variable","BioMas 1.3.5"))
    testthat::expect_match(html,label,fixed=TRUE)
})

testthat::test_that("screenshot KSS headers map across midnight after one confirmation", {
  shiny::testServer(server,{
    setup_inputs(session)
    columns<-paste0("KSS_",c(seq(10,24,2),seq(2,8,2),paste0(seq(10,20,2),"_SD")))
    upload(session,"kss",columns)
    testthat::expect_false(sim_inputs_live()$ok)
    testthat::expect_warning({
      testthat::expect_identical(output$index_loaded,"false")
      testthat::expect_null(output$cross_modality_section)
    },NA)
    testthat::expect_match(output$kss_numeric_time_proposal$html,"Use column numbers as clock hours",fixed=TRUE)
    testthat::expect_match(output$kss_numeric_time_proposal$html,"24 is midnight",fixed=TRUE)
    session$setInputs(kss_use_numeric_clock=1)
    r<-imported$KSS()
    testthat::expect_identical(r$settings$columns,columns)
    testthat::expect_equal(r$settings$recording_hours,seq(10,44,2))
    testthat::expect_equal(r$settings$times,seq(2,36,2))
    testthat::expect_match(r$settings$time_mapping,"User-confirmed numeric")
    testthat::expect_identical(output$index_loaded,"true")
    testthat::expect_true(any(grepl("Midnight adjustment",notice_log)))
    session$setInputs(btn_run=1)
    testthat::expect_equal(valid_bundle()$results$comparisons$KSS$cells$time_h,seq(2,36,2))
    # An equal pair is a quiet validation failure, not a guessed extra day.
    session$setInputs(kss_times_txt=paste(rep(10,length(columns)),collapse=", "))
    e<-tryCatch(imported$KSS(),error=function(e)e)
    testthat::expect_s3_class(e,"shiny.silent.error")
    testthat::expect_match(conditionMessage(e),"equal")
    testthat::expect_identical(output$index_loaded,"false")
    testthat::expect_null(output$cross_modality_section)
    testthat::expect_false(sim_inputs_live()$ok)
    testthat::expect_error(valid_bundle())
  })
})

testthat::test_that("all panels preserve edits and recover from incomplete numeric clock mapping", {
  shiny::testServer(server,{
    setup_inputs(session)
    for(prefix in c("kss","vas","pvt")){
      upload(session,prefix,paste0(toupper(prefix),c("_22","_24","_2")))
      do.call(session$setInputs,setNames(list(1),paste0(prefix,"_use_numeric_clock")))
      testthat::expect_identical(imported[[toupper(prefix)]]()$settings$times,c(14,16,18))
      do.call(session$setInputs,setNames(list("21, 0, 3",2),paste0(prefix,c("_times_txt","_use_numeric_clock"))))
      testthat::expect_identical(imported[[toupper(prefix)]]()$settings$clock_text,"21, 0, 3")
      do.call(session$setInputs,setNames(list("21, 0"),paste0(prefix,"_times_txt")))
      testthat::expect_false(sim_inputs_live()$ok)
      testthat::expect_warning({output$index_loaded;output$cross_modality_section},NA)
      do.call(session$setInputs,setNames(list("21, 0, 3"),paste0(prefix,"_times_txt")))
    }
    session$setInputs(btn_run=1)
    b<-valid_bundle()
    testthat::expect_true(all(c("KSS","VAS","PVT")%in%names(b$results$comparisons)))
    writeLines(fr_export_script(b),file.path(biomas_test_root,"work/biomas/reproduce_numeric.R"))
  })
})

testthat::test_that("composite plot always exposes aligned English KSS anchors", {
  shiny::testServer(server,{
    setup_inputs(session)
    session$setInputs(btn_run=1)
    for(mode in c("off","kss","risk")) {
      session$setInputs(risk_overlay_sleep=mode!="off",risk_overlay_scale=if(mode=="off")"kss"else mode)
      z<-jsonlite::fromJSON(output$plt_risk,simplifyVector=FALSE)$x
      testthat::expect_identical(z$layout$yaxis2$side,"right")
      testthat::expect_equal(unlist(z$layout$yaxis2$tickvals),1:9)
      labels<-unlist(z$layout$yaxis2$ticktext)
      testthat::expect_length(labels,9)
      testthat::expect_match(labels[1],"Extremely alert")
      testthat::expect_match(labels[9],"fighting sleep")
      testthat::expect_equal(unlist(z$layout$yaxis2$range),1+8*unlist(z$layout$yaxis$range))
      testthat::expect_true(z$layout$yaxis$fixedrange)
      testthat::expect_true(any(vapply(z$data,function(tr)identical(tr$yaxis,"y2"),logical(1))))
    }
    testthat::expect_match(as.character(ui),"not a validated conversion",fixed=TRUE)
  })
})
