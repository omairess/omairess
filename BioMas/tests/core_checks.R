options(biomas.launch=FALSE,warn=1)
source("app.R")
checks<-0L
ok<-function(value,label){if(!isTRUE(value))stop("FAILED: ",label,call.=FALSE);checks<<-checks+1L;cat("PASS ",label,"\n",sep="")}
near<-function(x,y,tol=1e-7) isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=tol,check.attributes=FALSE))
fails<-function(expr,pattern=NULL){e<-tryCatch({force(expr);NULL},error=function(e)e);inherits(e,"error")&&(is.null(pattern)||grepl(pattern,conditionMessage(e),ignore.case=TRUE))}

fixture<-function(){list(origin=fr_origin(),end=4.07,schedule=data.frame(onset_h=1.013,offset_h=2.017),
  duty=data.frame(onset_h=.25,offset_h=3.5),version="2024fe",overrides=list(PVT=list(),KSS=list()),
  baseline=list(onset=23,offset=7),initial=list(PVT=NULL,KSS=NULL),spacing=1/60,
  index=list(weights=c(.5,.5),pvt_ceiling=20),exposure=list(metric="Index",threshold=.3,scope="duty"),empirical=list())}

cat("\nMODEL AND INITIALIZATION\n")
old<-new.env(parent=globalenv())
for(e in parse("tests/fixtures/legacy_model_reference.R")) if(is.call(e)&&identical(e[[1]],as.name("<-"))&&is.symbol(e[[2]])&&
 as.character(e[[2]])%in%c("mccauley_params","mc2024_step","circadian_c"))eval(e,old)
set.seed(42);maxerr<-0
for(v in c("2021","2024jtb","2024fe"))for(m in c("PVT","KSS")){
  P<-fr_parameters(m,v);ref<-old$mccauley_params(m,v)
  ok(all(vapply(names(ref),function(nm)identical(P[[nm]],ref[[nm]]),logical(1))),paste("Published parameters:",v,m))
  for(w in c(TRUE,FALSE))for(i in 1:20){
    st<-list(p=runif(1,1,30),u=runif(1,1,40),h=if(P$use_inertia)runif(1,0,10)else 0,kappa=runif(1,.1,.9))
    y<-c(st$p-P$p_f-P$xi*st$h,st$u,st$h,st$kappa);t<-runif(1,0,24)
    dy<-fr_rhs(t,y,P,w)[[1]];dp<-c(dy[1]+P$xi*dy[3],dy[2:4])
    actual<-(unlist(old$mc2024_step(st,w,t,P,1e-4))-unlist(st))/1e-4
    maxerr<-max(maxerr,max(abs(dp-actual)))
  }
  b<-fr_baseline(P,t0_clock=7)
  ok(b$residual<=b$tolerance && b$cycle_multiplier<1,paste("Verified stable baseline:",v,m))
}
ok(maxerr<1e-8,"240 independent derivative comparisons")
cat("Maximum derivative difference:",maxerr,"\n")
for(m in c("PVT","KSS")){
  P<-mccauley_params(m);b<-fr_baseline(P,t0_clock=8);p<-P$p_f+b$state[1]+P$xi*b$state[3]
  expected<-if(m=="PVT")4.405841 else 3.908006
  ok(abs(p-expected)<1e-5,paste("08:00 initial phase corrected:",m))
  from7<-fr_simulate_arm(P,data.frame(onset_h=16,offset_h=24),26,7,extra_times=1)$predictions
  from8<-fr_simulate_arm(P,data.frame(onset_h=15,offset_h=23),25,8)$predictions
  ok(abs(from7$p[match(1,from7$time_h)]-from8$p[1])<1e-7,paste("Shifted study start matches continuous history:",m))
}
ok(fails(fr_baseline(mccauley_params("PVT"),onset=3.18,offset=7),"not a stable"),"Feedback stability guard catches W=20.18 below nominal Wc")
nearcritical<-fr_baseline(mccauley_params("PVT","2021"),onset=3.19,offset=7)
ok(nearcritical$residual<=nearcritical$tolerance,"Near-critical continuous baseline is verified, never silently truncated")
for(nm in c("alpha_s","eta_w","xi","lambda"))ok(fails(fr_parameters("PVT","2024fe",setNames(list(0),nm))),paste("Reject zero",nm))
ok(fails(fr_parameters("PVT","2024fe",list(W_c=24))),"Reject singular Wc=24")
ok(fails(fr_parameters("PVT","2024fe",list(tau=25))),"Published period cannot become a discontinuous 25-hour clock")
p<-fr_parameters("PVT","2021",list(gamma=.71,zeta=1.31));ok(p$gamma==0&&p$zeta==0&&!p$use_inertia,"2021 structural exclusions survive edits")
ok(fr_parameters("KSS","2024jtb",list(zeta=2))$zeta==0,"JTB structural feedback exclusion survives edits")

cat("\nEVENTS, DURATION AND EXPOSURE\n")
cfg<-fixture();b<-fr_capture_run(cfg);pred<-b$results$predictions
ok(all(c(0,1.013,2.017,.25,3.5,4.07)%in%pred$time_h),"All sleep/duty/end boundaries appear exactly")
ok(!pred$wake[match(1.013,pred$time_h)]&&pred$wake[match(2.017,pred$time_h)],"Sleep interval uses onset-inclusive/end-exclusive states")
ok(near(max(pred$time_h),4.07,1e-12),"Explicit study end is included")
empty<-cfg;empty$schedule<-fr_empty_intervals();empty$duty<-fr_empty_intervals();empty$exposure$scope<-"wake"
ok(all(fr_run(empty)$results$predictions$wake),"Continuous-wake studies need no dummy sleep interval")
ok(fails(fr_intervals(data.frame(onset_h=c(1,1.5),offset_h=c(2,3)),4)),"Reject overlapping sleep intervals")
ok(fails(fr_intervals(data.frame(onset_h=2,offset_h=1),4)),"Reject reversed intervals")
ok(fails(fr_intervals(data.frame(onset_h=0,offset_h=5),4)),"Reject intervals beyond explicit end")
ok(near(fr_positive_segment(0,10,1,5),c(.5,1.25)),"Analytic rising threshold crossing")
ok(near(fr_positive_segment(10,0,1,5),c(.5,1.25)),"Analytic falling threshold crossing")
ok(near(fr_positive_segment(2,2,1,1),c(1,1)),"One-hour constant exposure has no endpoint overcount")
ok(near(fr_positive_segment(1,1,1,1),c(0,0)),"Exactly at threshold is not strictly above")
I<-fr_index_settings(list(weights=c(0,1),pvt_ceiling=20))
lin<-data.frame(time_h=c(0,1),PVT_raw=c(-10,30),KSS_raw=c(1,1))
e<-fr_exposure(lin,fr_empty_intervals(),fr_empty_intervals(),list(metric="Index",threshold=.5,scope="wake"),I)$summary
ok(near(e$hours_above_threshold,.5)&&near(e$excess_area,.1875),"Composite AUC splits clipping corners before threshold integration")
lin<-data.frame(time_h=c(0,.25,.75,1),PVT_raw=c(0,2.5,7.5,10),KSS_raw=1)
e<-fr_exposure(lin,fr_empty_intervals(),data.frame(onset_h=.25,offset_h=.75),list(metric="PVT_raw",threshold=5,scope="duty"),I)$summary
ok(near(e$exposure_hours,.5)&&near(e$hours_above_threshold,.25)&&near(e$excess_area,.3125),"Duty mask restricts both exposure duration and severity")
e<-fr_exposure(lin,data.frame(onset_h=0,offset_h=1),fr_empty_intervals(),list(metric="Index",threshold=.5,scope="wake"),I)$summary
ok(e$exposure_hours==0&&e$excess_area==0&&is.na(e$mean_index),"Zero awake exposure has defined area and undefined mean")
ok(fails(fr_index_settings(list(weights=c(0,0),pvt_ceiling=20))),"Reject zero-sum weights")
ok(fails(fr_index_settings(list(weights=c(-1,2),pvt_ceiling=20))),"Reject negative weights")
ok(fails(fr_index_settings(list(weights=c(1,1),pvt_ceiling=0))),"Reject zero index ceiling")
ok(near(fr_index(10,5,fr_index_settings(list(weights=c(1,1),pvt_ceiling=20))),.5),"Composite normalization and units")
ok(near(fr_index_settings(list(weights=c(1e308,1e308),pvt_ceiling=20))$normalized_weights,c(.5,.5)),"Weight normalization remains finite for large positive weights")

cat("\nINPUTS AND OBSERVATIONS\n")
tf<-tempfile(fileext=".csv");writeLines(c("id;kss1;kss2","1;4,5;5,5","2;3,5;4,5"),tf)
d<-fr_read_table(tf,"values.csv","semicolon")
ok(identical(fr_parse_numbers(d$kss1,","),c(4.5,3.5)),"Decimal-comma files preserve 4.5 rather than producing 45")
ok(fails(fr_parse_numbers(d$kss1,".")),"Wrong numeric locale is rejected, not silently converted")
ok(fails(fr_parse_numbers(c("4","four"))),"Reject invalid numeric text")
ok(fails(fr_parse_numbers(c(1,Inf))),"Reject infinite observations")
ok(fails(fr_parse_numbers("1,000",".")),"Reject ambiguous grouping separators")
bad_csv<-tempfile(fileext=".csv");writeLines(c("id,value","1,2,3"),bad_csv)
ok(fails(suppressWarnings(fr_read_table(bad_csv,"broken.csv","comma"))),"Reject malformed CSV field counts")
header_csv<-tempfile(fileext=".csv")
writeLines(c(",vas1,vas1,vas1__1,column_1", "s1,20,30,40,50", "s2,60,70,80,90"),header_csv)
repaired<-fr_read_table(header_csv,"vas.csv","comma")
ok(identical(names(repaired),c("column_1","vas1","vas1__2","vas1__1","column_1__1")),"Blank/duplicate headers repaired without colliding with existing suffixes")
ok(ncol(repaired)==5L&&nrow(repaired)==2L&&identical(unname(unlist(repaired[1,])),c("s1","20","30","40","50")),"Header repair preserves every column and value in original order")
header_map<-attr(repaired,"column_name_map")
ok(identical(header_map$original,c("","vas1","vas1","vas1__1","column_1"))&&identical(header_map$position,1:5),"Header mapping retains original labels and column positions")
ok(identical(header_map$changed,c(TRUE,FALSE,TRUE,FALSE,TRUE)),"Only changed headers are marked")
ok(identical(fr_read_table(header_csv,"vas.csv","comma"),repaired),"Header repair is deterministic across imports")
edge<-data.frame(a=1,b=2,c=3,d=4);names(edge)<-c("\ufeffid"," ",NA,"id")
edge<-fr_repair_headers(edge)
ok(!anyNA(names(edge))&&all(nzchar(names(edge)))&&!anyDuplicated(names(edge)),"BOM, whitespace and missing headers become usable names")
ok(identical(names(fr_repair_headers(data.frame(id=1,vas1=2))),c("id","vas1")),"Valid headers stay unchanged")
header_rds<-tempfile();saveRDS(repaired,header_rds)
ok(identical(readRDS(header_rds),repaired),"Repaired headers and their original mapping survive saved exports")
ok(fails(fr_repair_headers(data.frame())),"Files with no columns still fail clearly")
excel_fixture<-system.file("extdata","datasets.xlsx",package="readxl")
excel<-fr_read_table(excel_fixture,"datasets.xlsx")
numeric_column<-which(vapply(excel,is.numeric,logical(1)))[1]
ok(!is.na(numeric_column)&&identical(fr_parse_numbers(excel[[numeric_column]],","),as.numeric(excel[[numeric_column]])),"Native Excel numbers are independent of text decimal locale")
sav<-tempfile(fileext=".sav");haven::write_sav(data.frame(id=c("s1","s2"),kss1=c(4.5,5.5)),sav)
ok(near(fr_read_table(sav,"test.sav")$kss1,c(4.5,5.5)),"SPSS numeric observations retain their values")
o<-fr_origin("dated",start_iso="2026-03-29T01:30:00+01:00")
ok(near(fr_elapsed("2026-03-29T03:30:00+02:00","iso",o),1),"DST change preserves elapsed UTC duration")
ok(near(o$t0_clock,1.5),"Datetime origin carries its actual local starting phase")
ok(fails(fr_parse_iso("2026-10-25 02:30:00")),"Reject ambiguous local datetime without offset")
ok(fails(fr_parse_iso("2026-15-16T08:00:00Z")),"Reject invalid calendar date")
ok(fails(fr_parse_iso("2026-09-16T08:00:00+99:00")),"Reject invalid UTC offsets")
ok(fails(fr_parse_iso("2026-09-16T25:00:00Z")),"Reject invalid ISO clock hours")
ok(near(fr_elapsed("2026-03-29T03:30+02:00","iso",o),1),"ISO datetimes can omit seconds without changing the origin")
act<-fr_actigraphy(c(0,1,2),c(1,0,0),1,3)
ok(identical(act,data.frame(onset_h=1,offset_h=3)),"Terminal sleep includes the confirmed final epoch")
ok(nrow(fr_actigraphy(c(0,1,2),c(1,1,1),1,3))==0,"All-wake actigraphy yields an empty interval table")
ok(fails(fr_actigraphy(c(0,2),c(1,0),1,3)),"Reject gaps in actigraphy")
ok(fails(fr_actigraphy(c(0,1),c(1,NA),1,2)),"Reject unknown actigraphy states")
ok(fails(fr_actigraphy(c(0,1),c(1,0),1,3)),"Reject unobserved trailing study time")
settings<-list(columns=c("kss1","kss2"),times=c(.5,1.5),id="id",group="__none__",decimal=",")
obs<-fr_prepare_observations(d,settings,"KSS")
ok(near(obs$value,c(4.5,3.5,5.5,4.5)),"Strict wide-to-long conversion preserves values")
bad<-settings;bad$times<-c(1,1);ok(fails(fr_prepare_observations(d,bad,"KSS")),"Reject duplicate subject/occasion times")
badraw<-d;badraw$id<-c(1,1);ok(fails(fr_prepare_observations(badraw,settings,"KSS")),"Reject duplicated participant rows")
badraw<-d;badraw$group<-c("A",NA);bad<-settings;bad$group<-"group"
kept<-fr_prepare_observations(badraw,bad,"KSS");ok(nrow(kept)==2 && all(kept$group=="A"),"Preserve supplied app: omit rows with missing group labels")
badraw<-d;badraw$kss1[1]<-"10";ok(fails(fr_prepare_observations(badraw,settings,"KSS")),"Reject out-of-range KSS")
vas<-list(columns="vas1",times=.5,id="id",group="__none__",decimal=".",vas_min=0,vas_max=100,direction="direct")
v<-fr_prepare_observations(data.frame(id=1:3,vas1=c(0,50,100)),vas,"VAS")
ok(near(v$value,c(1,5,9)),"VAS endpoint rescaling")
vas$direction<-"reverse";v<-fr_prepare_observations(data.frame(id=1:3,vas1=c(0,50,100)),vas,"VAS")
ok(near(v$value,c(9,5,1)),"Reverse VAS endpoint rescaling")
ok(fails(fr_prepare_observations(data.frame(id=1,vas1=101),vas,"VAS")),"VAS out-of-range values are not clipped")
pvt<-list(columns="pvt1",times=0,id="id",group="__none__",decimal=".",pvt_units_confirmed=FALSE,lapse_rule=">=500 ms")
ok(fails(fr_prepare_observations(data.frame(id=1,pvt1=250),pvt,"PVT")),"PVT comparison requires explicit lapse-count unit confirmation")

cat("\nCOMPARISONS, UI STATE AND EXPORT\n")
observations<-data.frame(subject=rep(c("s1","s2"),5),group="A",occasion=rep(paste0("t",1:5),each=2),
  time_h=rep(c(-1,.5,1.5,2.5,5),each=2),value_raw=rep(c(4,5),5),value=rep(c(4,5),5),metric="KSS")
cfg<-fixture();cfg$empirical<-list(KSS=list(observations=observations,selected_groups="A",settings=list(),source=list(raw_table=observations)))
b<-fr_capture_run(cfg);cmp<-b$results$comparisons$KSS
ok(nrow(cmp$cells)==2&&nrow(cmp$exclusions)==6,"Comparison excludes sleep and out-of-coverage observations")
ok(identical(sort(unique(cmp$exclusions$exclusion)),sort(c("Outside study coverage","During sleep"))),"Exclusion reasons retained")
ok(all(cmp$cells$n==2)&&near(cmp$cells$SEM,rep(.5,2)),"SD/SEM use valid independent subject rows")
none<-fr_compare(b$results$predictions,observations,character(),cfg$schedule,cfg$end,"KSS")
ok(nrow(none$cells)==0&&all(none$statistics$pairs==0),"Empty group selection includes no observations")
ok(near(fr_fit_statistics(3,1)[c("RMSE","MAE")],c(2,2)),"Single-pair RMSE and MAE are available")
ok(is.na(fr_fit_statistics(c(2,2),c(1,2))["Pearson_r"]),"Constant predictions have undefined correlation")
ok(near(fr_fit_statistics(c(2,3),c(1,2))["Pearson_r"],1),"Correlation uses the standard formula")


cat("ALL",checks,"CORE CHECKS PASSED\n")
