# ORCHID | IJMI reproducibility analysis | cleaned from verified V3.6
# Analytical calculations below are retained from IJMI_Final_Analysis_V3_6(2).R.
# R >= 4.3. Required packages: data.table, dplyr, ggplot2, pROC, xgboost, Matrix.
# Run in RStudio with Source. Change INPUT and ROOT_OUT here if needed.
# IMPORTANT: "primary" in legacy code means EXPANDED predictors; "reduced" is
# the manuscript's primary specification (age, gender, race, OPO).
# The script does not download or distribute the restricted ORCHID dataset.
options(stringsAsFactors=FALSE)
required_packages <- c("data.table", "dplyr", "ggplot2", "pROC", "xgboost", "Matrix")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing_packages)) stop("Install missing R packages: ", paste(missing_packages, collapse=", "))
suppressPackageStartupMessages(invisible(lapply(required_packages, library, character.only=TRUE)))
# A local OPOReferrals.csv in the working directory is used by default.
# If absent, select the source file interactively in RStudio.
INPUT <- "OPOReferrals.csv"
if(!file.exists(INPUT) && interactive()) INPUT <- file.choose()
# A new output folder is created under the current working directory.
# Set an absolute path here to write somewhere else.
ROOT_OUT <- file.path(getwd(), "ORCHID_V3_6_Outputs")
BOOT_B <- 2000L # manuscript bootstrap repetitions
MIN_GROUP_N <- 30L
SEED <- 123L
RUN_MODES <- c("random", "temporal")
RUN_VARIANTS <- c("primary", "reduced")
# The reduced model excludes mechanism/circumstances/height/weight and missingness flags.
# Primary model retains them provisionally for comparison, NOT proof of timing safety.
# Complete-case sensitivity is an optional third run; it changes the target population.
RUN_COMPLETE_CASE <- TRUE
# Exclude referrals with logically inconsistent stage outcomes (27 in original data).
EXCLUDE_INCONSISTENT <- TRUE
# Set FALSE for a fast diagnostic test, TRUE for final model runs.
SAVE_MODELS <- TRUE
set.seed(SEED)
if(!file.exists(INPUT)) stop("Place OPOReferrals.csv beside this script or set INPUT to its path.")
dir.create(ROOT_OUT,showWarnings=FALSE,recursive=TRUE)
run_one <- function(SPLIT_MODE, VARIANT) {
  OUT <- file.path(ROOT_OUT,paste0(SPLIT_MODE,"_",VARIANT))
  dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
  set.seed(SEED)
  save_csv <- function(x,name) fwrite(as.data.table(x),file.path(OUT,name))
  logline <- function(...) cat(...,"\n",file=file.path(OUT,"run_log.txt"),append=TRUE)
  cat("Running:",SPLIT_MODE,VARIANT,"\n")
d <- fread(INPUT,na.strings=c('','NA','NULL','NaN'),showProgress=FALSE)
required <- c('approached','authorized','procured','age','gender','race','opo','referral_year')
if(length(setdiff(required,names(d)))) stop('Missing columns: ',paste(setdiff(required,names(d)),collapse=', '))
raw_n <- nrow(d); d[,row_id:=seq_len(.N)]
# Explicit binary encoding; unexpected values become NA and are audited
as01 <- function(z){z<-toupper(trimws(as.character(z))); fifelse(z %in% c('1','TRUE','T','YES','Y'),1L,fifelse(z %in% c('0','FALSE','F','NO','N'),0L,NA_integer_))}
for(v in c('approached','authorized','procured')) d[,(v):=as01(get(v))]
# Data audit BEFORE exclusions or imputation
original_columns <- setdiff(names(d),'row_id')
miss <- data.table(variable=original_columns,missing_n=sapply(d[,..original_columns],function(x) sum(is.na(x))),total_n=raw_n)
miss[,missing_pct:=round(100*missing_n/total_n,3)];save_csv(miss,'01_missingness_all_variables.csv')
contradiction <- d[(authorized==1 & approached!=1) | (procured==1 & (authorized!=1 | approached!=1))]
save_csv(contradiction[,.(row_id,patient_id,opo,referral_year,approached,authorized,procured)],'02_stage_inconsistencies.csv')
# Plausibility thresholds are SCREENING RULES, not verified clinical definitions.
# Flag values and replace with NA; sensitivity analysis should test alternatives.
for(v in intersect(c('age','height_in','weight_kg'),names(d))) d[,(v):=suppressWarnings(as.numeric(get(v)))]
flags <- data.table(row_id=d$row_id,age_flag=!is.na(d$age)&(d$age<0|d$age>120),height_flag=!is.na(d$height_in)&(d$height_in<12|d$height_in>96),weight_flag=!is.na(d$weight_kg)&(d$weight_kg<1|d$weight_kg>500))
save_csv(flags[age_flag|height_flag|weight_flag],'03_implausible_measurements.csv')
d[!is.na(age)&(age<0|age>120),age:=NA_real_]
d[!is.na(height_in)&(height_in<12|height_in>96),height_in:=NA_real_]
d[!is.na(weight_kg)&(weight_kg<1|weight_kg>500),weight_kg:=NA_real_]
if(EXCLUDE_INCONSISTENT) d <- d[!row_id %in% contradiction$row_id]
# Keep referral cohort intact; eligibility for stage 2 and 3 enforced below
# Split on referral-level records. Duplicate patient_id grouped if present.
# The random split is not external validation; temporal split is a robustness test.
# Patient-level split: duplicate patient IDs are never split across random partitions.
yr <- suppressWarnings(as.integer(as.character(d$referral_year)))
if(anyNA(yr) || any(!yr %in% 2015:2021)) stop("Unexpected referral years")
if(SPLIT_MODE=="temporal") {
  d[,partition:=fifelse(yr<=2019,"train",fifelse(yr==2020,"calibration","test"))]
} else if(SPLIT_MODE=="random") {
  ids <- if("patient_id" %in% names(d)) as.character(d$patient_id) else as.character(d$row_id)
  ids[is.na(ids)|ids==""] <- paste0("missing_id_",d$row_id[is.na(ids)|ids==""])
  u <- unique(ids); set.seed(SEED); u <- sample(u)
  ntr <- floor(.60*length(u)); nca <- floor(.20*length(u))
  tr <- u[seq_len(ntr)]; ca <- u[seq.int(ntr+1L,ntr+nca)]
  d[,partition:=fifelse(ids %in% tr,"train",fifelse(ids %in% ca,"calibration","test"))]
} else stop("Unknown split")
# Check year-by-OPO coverage: OPO2 is absent until 2018, which may affect transportability.
year_opo <- d[,.(n=.N, approached_n=sum(approached==1,na.rm=TRUE), authorized_n=sum(authorized==1 & approached==1,na.rm=TRUE), procured_n=sum(procured==1 & authorized==1 & approached==1,na.rm=TRUE)),by=.(referral_year,opo)]
save_csv(year_opo,'00_year_by_opo_counts.csv')
year_cohort <- d[,.(n=.N, approached_n=sum(approached==1,na.rm=TRUE), authorized_n=sum(authorized==1 & approached==1,na.rm=TRUE), procured_n=sum(procured==1 & authorized==1 & approached==1,na.rm=TRUE)),by=.(referral_year,partition)]
save_csv(year_cohort,'00_year_stage_counts.csv')
# Check for duplicated patient IDs across chronological partitions, if patient_id exists.
if('patient_id' %in% names(d)) {
  cross <- d[!is.na(patient_id),.(partitions=uniqueN(partition),n=.N),by=patient_id][partitions>1]
  save_csv(cross,'00_patient_ids_across_partitions.csv')
  if(nrow(cross)) warning('Patient IDs occur across temporal partitions; evaluate patient-level dependence.')
}
# Features are intentionally limited to referral descriptors; no stage outcomes/timestamps.
# For temporal tests, referral_year is EXCLUDED to avoid extrapolation to unseen year 2021.
base_features <- c("age","gender","race","opo")
questionable <- c("height_in","weight_kg","mechanism_of_death","circumstances_of_death")
if(VARIANT=="primary" || VARIANT=="complete_case") features <- c(base_features,questionable)
else if(VARIANT=="reduced") features <- base_features
else stop("Unknown variant")
if(SPLIT_MODE=="random") features <- c(features,"referral_year")
features <- intersect(features,names(d))
if(VARIANT=="complete_case") {
  keep <- complete.cases(d[,..questionable])
  d <- d[keep]
  # IMPORTANT: complete-case selection is outcome-associated; diagnostic only.
  save_csv(data.table(complete_case_n=nrow(d),excluded_n=raw_n-nrow(d)),"00_complete_case_selection.csv")
}
numerics <- intersect(c("age","height_in","weight_kg","referral_year"),features)
categoricals <- setdiff(features,numerics)
# Original (unimputed) age and gender retained for subgroup reporting.
d[,age_for_group:=age]
d[,gender_for_group:=as.character(gender)]
d[,race_for_group:=as.character(race)]
save_csv(data.table(feature=features,variant=VARIANT,split=SPLIT_MODE),"03_model_feature_manifest.csv")
# Audit whether clinical-variable missingness predicts the downstream outcome.
for(v in intersect(c("mechanism_of_death","circumstances_of_death","height_in","weight_kg"),names(d))) {
  z <- data.table(
    missing = is.na(d[[v]]),
    approached = d[["approached"]],
    authorized = d[["authorized"]],
    procured = d[["procured"]],
    partition = d[["partition"]]
  )
  stopifnot(nrow(z) == nrow(d))
  a <- z[,.(n=.N,approach_rate=mean(approached==1,na.rm=TRUE),
    procurement_rate=mean(procured==1,na.rm=TRUE)),by=.(partition,missing)]
  a[,variable:=v]
  save_csv(a,paste0("03b_missingness_outcome_",v,".csv"))
}
# Training-only imputation and shared factor levels; model.matrix on all partitions uses same columns.
for (v in numerics) {
  # Work only with the post-exclusion data.table, never the original 133,101-row vector.
  values <- d[[v]]
  if (length(values) != nrow(d)) stop('Length mismatch before imputation: ', v)
  missing_idx <- which(is.na(values))
  set(d, j=paste0(v, '_missing'), value=as.integer(is.na(values)))
  train_values <- values[d[['partition']] == 'train']
  med <- suppressWarnings(median(train_values, na.rm=TRUE))
  if (!is.finite(med)) med <- 0
  if (length(missing_idx)) set(d, i=missing_idx, j=v, value=med)
  stopifnot(nrow(d) == length(d[[paste0(v, '_missing')]]))
}
for(v in categoricals){tmp<-as.character(d[[v]]);tmp[is.na(tmp)|trimws(tmp)=='']<-'Missing';
  # unseen levels are represented as OtherUnseen to prevent accidental model.matrix mismatch
  seen<-unique(tmp[d$partition=='train']);tmp[!tmp%in%seen]<-'OtherUnseen';
  d[,(v):=factor(tmp,levels=unique(c(seen,'OtherUnseen')))]}
features<-c(features,paste0(numerics,'_missing'))
# Avoid automatic target leakage from stage indicators or downstream timestamps: not in feature list.
# Stage eligibility is OBSERVED history, not model-predicted eligibility.
stage_specs <- list(approach=list(y='approached',eligible=rep(TRUE,nrow(d))),authorization=list(y='authorized',eligible=d$approached==1),procurement=list(y='procured',eligible=d$approached==1 & d$authorized==1))
cohort <- rbindlist(lapply(names(stage_specs),function(s){a<-stage_specs[[s]];rbindlist(lapply(c('train','calibration','test'),function(p){idx<-which(!is.na(a$eligible)&a$eligible & d$partition==p & !is.na(d[[a$y]]));data.table(stage=s,partition=p,n=length(idx),positive=sum(d[[a$y]][idx]==1),prevalence=if(length(idx))mean(d[[a$y]][idx]==1) else NA_real_)}))}))
save_csv(cohort,'04_stage_cohorts_prevalence.csv')
# Shared sparse design matrix built with training-defined factor levels; no outcome information.
fml <- reformulate(features)
X <- sparse.model.matrix(fml,d)[,-1,drop=FALSE]
# Prediction, calibration and thresholds are always fitted without test outcomes.
clip <- function(x) pmin(pmax(as.numeric(x),1e-6),1-1e-6)
fit_platt <- function(raw,y){fit<-glm(y~qlogis(clip(raw)),family=binomial());function(p) clip(plogis(coef(fit)[1]+coef(fit)[2]*qlogis(clip(p))))}
threshold_youden <- function(y,p){rr<-pROC::roc(y,p,quiet=TRUE,direction='<');as.numeric(pROC::coords(rr,x='best',best.method='youden',ret='threshold',transpose=FALSE)[1])}
score <- function(y,p,threshold){p<-clip(p); pred<-as.integer(p>=threshold); data.table(n=length(y),prevalence=mean(y),auroc=if(length(unique(y))==2)as.numeric(pROC::auc(pROC::roc(y,p,quiet=TRUE,direction='<')))else NA_real_,brier=mean((p-y)^2),calibration_intercept=tryCatch(unname(coef(glm(y~offset(qlogis(p)),family=binomial()))[1]),error=function(e)NA_real_),calibration_slope=tryCatch(unname(coef(glm(y~qlogis(p),family=binomial()))[2]),error=function(e)NA_real_),sensitivity=if(sum(y==1)>0)mean(pred[y==1]==1)else NA_real_,specificity=if(sum(y==0)>0)mean(pred[y==0]==0)else NA_real_,ppr=mean(pred==1))}
# PR AUC by precision-recall step integration (average precision)
ap <- function(y,p){ord<-order(p,decreasing=TRUE);y<-y[ord];if(sum(y)==0)return(NA_real_);sum((cumsum(y)/seq_along(y))[y==1])/sum(y)}
models<-list(); pred_rows<-list();metrics<-list(); curves<-list(); threshold_table<-list()
for(s in names(stage_specs)){
  spec<-stage_specs[[s]];elig<-which(!is.na(spec$eligible)&spec$eligible & !is.na(d[[spec$y]]))
  it<-elig[d$partition[elig]=='train'];ic<-elig[d$partition[elig]=='calibration'];ie<-elig[d$partition[elig]=='test']
  if(min(length(it),length(ic),length(ie))<50)stop('Insufficient cohort for ',s)
  yt<-d[[spec$y]][it];yc<-d[[spec$y]][ic];ye<-d[[spec$y]][ie]
  if(any(sapply(list(yt,yc,ye),function(y)length(unique(y))<2)))stop('Single-class partition for ',s)
  # Unweighted LR: generally more directly interpretable as probabilities
  # Fit LR on the SHARED design matrix: prevents unseen factor-level errors.
  # Non-estimable coefficients (rank-deficient columns) are set to zero for prediction.
  lr_xtrain <- cbind('(Intercept)'=1,as.matrix(X[it,,drop=FALSE]))
  glm_fit <- glm.fit(x=lr_xtrain,y=yt,family=binomial(),control=glm.control(maxit=100))
  lr_beta <- coef(glm_fit); lr_beta[!is.finite(lr_beta)] <- 0
  lr_predict <- function(indices) clip(plogis(as.vector(cbind(1,as.matrix(X[indices,,drop=FALSE])) %*% lr_beta)))
  lr_cal <- lr_predict(ic)
  lr_test <- lr_predict(ie)
  # XGB early stopping uses an INTERNAL subset of train, not calibration/test.
  set.seed(100+match(s,names(stage_specs)));va<-sample(seq_along(it),size=max(50,floor(.2*length(it))))
  fitidx<-setdiff(seq_along(it),va)
  xparams<-list(objective='binary:logistic',eval_metric='logloss',eta=.05,max_depth=4,min_child_weight=10,subsample=.8,colsample_bytree=.8,nthread=2)
  xfit<-xgb.train(params=xparams,data=xgb.DMatrix(X[it[fitidx],,drop=FALSE],label=yt[fitidx]),nrounds=1000,evals=list(validation=xgb.DMatrix(X[it[va],,drop=FALSE],label=yt[va])),early_stopping_rounds=35,verbose=0)
  # refit on all training data for best iteration, no calibration/test exposure
  # XGBoost 3.1.x stores best iteration in early_stop (1-based), not at top level.
  # In xgboost 3.1.2, early_stop is an R OBJECT ATTRIBUTE, not xfit$early_stop.
  es <- attr(xfit, 'early_stop', exact=TRUE)
  rounds <- if (!is.null(es)) es$best_iteration else NULL
  xa <- xgboost::xgb.attributes(xfit)
  if (is.null(rounds) || length(rounds) != 1L || !is.finite(as.numeric(rounds)) || as.numeric(rounds) < 1) {
    # Booster attribute is zero-based; training-round count is one-based.
    if (!is.null(xa$best_iteration) && length(xa$best_iteration)==1L)
      rounds <- as.integer(xa$best_iteration) + 1L
  }
  if (is.null(rounds) || length(rounds)!=1L || !is.finite(as.numeric(rounds)) || as.numeric(rounds)<1L)
    stop('Cannot retrieve valid XGBoost best iteration for stage ',s)
  rounds <- as.integer(rounds)
  best_score <- if (!is.null(es) && length(es$best_score)) as.numeric(es$best_score) else as.numeric(xa$best_score)
  stopped_max <- if (!is.null(es) && length(es$stopped_by_max_rounds)) as.logical(es$stopped_by_max_rounds) else NA
  save_csv(data.table(stage=s, best_iteration=rounds,
                      best_score=best_score,
                      stopped_by_max_rounds=stopped_max),
           paste0('04b_xgb_early_stopping_',s,'.csv'))
  cat('XGBoost:',s,'best iteration:',rounds,'best validation logloss:',
      as.numeric(best_score),'\n')
  xfinal<-xgb.train(params=xparams,data=xgb.DMatrix(X[it,,drop=FALSE],label=yt),nrounds=rounds,verbose=0)
  xcal_raw<-clip(predict(xfinal,X[ic,,drop=FALSE])); xtest_raw<-clip(predict(xfinal,X[ie,,drop=FALSE]))
  platt<-fit_platt(xcal_raw,yc); xcal<-platt(xcal_raw);xtest<-platt(xtest_raw)
  models[[s]]<-list(lr=glm_fit,xgb=xfinal,platt=platt,features=features)
  for(m in c('LR','XGB_raw','XGB_platt')){
    pc<-switch(m,LR=lr_cal,XGB_raw=xcal_raw,XGB_platt=xcal)
    pe<-switch(m,LR=lr_test,XGB_raw=xtest_raw,XGB_platt=xtest)
    th<-threshold_youden(yc,pc);threshold_table[[paste(s,m)]]<-data.table(stage=s,model=m,threshold=th)
    met<-score(ye,pe,th);met[,`:=`(stage=s,model=m,auprc_ap=ap(ye,pe),threshold=th)]
    metrics[[paste(s,m)]]<-met
    pred_rows[[paste(s,m)]]<-data.table(row_id=d$row_id[ie],stage=s,model=m,y=ye,probability=pe,threshold=th,prediction=as.integer(pe>=th),age=d$age_for_group[ie],race=d$race_for_group[ie],gender=d$gender_for_group[ie],opo=as.character(d$opo[ie]))
    caldf<-data.table(p=pe,y=ye);caldf[,bin:=pmin(10L,as.integer(cut(p,breaks=unique(quantile(p,probs=seq(0,1,.1))),include.lowest=TRUE,labels=FALSE)))]
    curves[[paste(s,m)]]<-caldf[!is.na(bin),.(mean_pred=mean(p),observed=mean(y),n=.N),by=bin][,`:=`(stage=s,model=m)]
  }
}
preds<-rbindlist(pred_rows); save_csv(preds,'05_test_predictions.csv')
met<-rbindlist(metrics);setcolorder(met,c('stage','model',setdiff(names(met),c('stage','model'))));save_csv(met,'06_test_model_performance.csv')
save_csv(rbindlist(threshold_table),'07_calibration_selected_thresholds.csv')
calcurves<-rbindlist(curves);save_csv(calcurves,'08_calibration_bins.csv')
p<-ggplot(calcurves,aes(mean_pred,observed,color=model))+geom_point()+geom_line(aes(group=interaction(stage,model)))+geom_abline(slope=1,intercept=0,linetype=2)+facet_wrap(~stage)+coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_bw()+labs(title='Held-out calibration',x='Predicted probability',y='Observed event rate');ggsave(file.path(OUT,'09_calibration.png'),p,width=11,height=4,dpi=200)
# Subgroup fairness, on stage-eligible test samples; confidence intervals below
# Group assignment uses explicit table columns, not data.table get() scoping.
preds[,age_group:=cut(age,c(-Inf,39,64,Inf),labels=c('<40','40-64','65+'))]
fair <- rbindlist(lapply(c("age_group","gender","race"),function(attr) {
  zz <- copy(preds)
  set(zz, j="group", value=as.character(zz[[attr]]))
  zz <- zz[!is.na(group) & trimws(group)!="" & !tolower(trimws(group)) %in% c("missing","na","n/a","not reported")]
  if(attr=="gender") zz <- zz[group %in% c("F","M")]
  a <- zz[,.(n=.N,observed_rate=mean(y),mean_pred=mean(probability),
    auroc=if(uniqueN(y)==2) as.numeric(pROC::auc(pROC::roc(y,probability,quiet=TRUE,direction="<"))) else NA_real_,
    ppr=mean(prediction==1),brier=mean((probability-y)^2),
    sensitivity=if(sum(y==1)>0)mean(prediction[y==1]==1) else NA_real_,
    specificity=if(sum(y==0)>0)mean(prediction[y==0]==0) else NA_real_),by=.(stage,model,group)]
  a[,attribute:=attr]
  a[n>=MIN_GROUP_N]
}))
save_csv(fair,'10_subgroup_fairness.csv')
# Stage-wise subgroup disparity CIs (conditional on fixed fitted model).
set.seed(SEED+1000L)
stage_boot <- rbindlist(lapply(c("age_group","race","gender"),function(attr) {
  rbindlist(lapply(c("LR","XGB_platt"),function(model_name) {
    rbindlist(lapply(names(stage_specs),function(stage_name) {
      z <- preds[stage==stage_name & model==model_name]
      set(z, j="group", value=as.character(z[[attr]]))
      z <- z[!is.na(group) & !tolower(trimws(group)) %in% c("missing","na","n/a","not reported")]
      if(attr=="gender") z <- z[group %in% c("F","M")]
      sz <- z[, .N, by=group][N>=MIN_GROUP_N,group]
      z <- z[group %in% sz]
      if(uniqueN(z$group)<2) return(NULL)
      gv <- split(z$probability,z$group)
      if(attr=="gender" && !all(c("F","M") %in% names(gv))) return(NULL)
      contrast <- function(x) { mu <- vapply(x,mean,numeric(1)); if(attr=="gender") unname(mu["M"]-mu["F"]) else diff(range(mu)) }
      est <- contrast(gv)
      draws <- replicate(BOOT_B,contrast(lapply(gv,function(x) sample(x,length(x),replace=TRUE))))
      data.table(stage=stage_name,attribute=attr,model=model_name,estimate=est,
        lo=unname(quantile(draws,.025)),hi=unname(quantile(draws,.975)),B=BOOT_B)
    }))
  }))
}))
save_csv(stage_boot,"10b_stagewise_gap_bootstrap.csv")
# Subgroup reliability bins (descriptive, n>=10 per bin).
rel <- copy(preds)
set(rel, j="group", value=rep(NA_character_, nrow(rel)))
reliability <- rbindlist(lapply(c("age_group","gender","race"),function(attr){
  zz<-copy(preds);set(zz, j="group", value=as.character(zz[[attr]]))
  zz<-zz[!is.na(group) & !tolower(trimws(group)) %in% c("missing","na","n/a","not reported")]
  if(attr=="gender")zz<-zz[group %in% c("F","M")]
  zz[,bin:=pmin(5L,as.integer(floor(probability*5)+1L))]
  out<-zz[,.(n=.N,mean_pred=mean(probability),observed=mean(y)),by=.(stage,model,group,bin)]
  out[,attribute:=attr];out[n>=10]
}))
save_csv(reliability,"10c_subgroup_reliability_bins.csv")
if(nrow(reliability)) {
  pp<-ggplot(reliability,aes(mean_pred,observed,color=group))+geom_point()+geom_line(aes(group=interaction(attribute,stage,model,group)))+geom_abline(slope=1,intercept=0,linetype=2)+facet_grid(attribute+stage~model)+theme_bw()+labs(x="Predicted probability",y="Observed frequency",title="Subgroup reliability: fixed-width bins")
  ggsave(file.path(OUT,"10d_subgroup_reliability.png"),pp,width=13,height=13,dpi=160)
}
# Pipeline probabilities must be generated for ALL held-out referrals, not only stage-eligible rows.
# The stage models are trained conditional on eligibility but scored counterfactually for all referrals.
# This is a model-based extrapolation for downstream-ineligible referrals; flag as limitation.
testidx<-which(d$partition=='test');r<-data.table(row_id=d$row_id[testidx],age=d$age_for_group[testidx],gender=d$gender_for_group[testidx],race=d$race_for_group[testidx],opo=as.character(d$opo[testidx]),approached=d$approached[testidx],authorized=d$authorized[testidx],procured=d$procured[testidx])
for(s in names(stage_specs)){
  mo<-models[[s]]
  lr_beta <- coef(mo$lr); lr_beta[!is.finite(lr_beta)] <- 0
  lp <- clip(plogis(as.vector(cbind(1,as.matrix(X[testidx,,drop=FALSE])) %*% lr_beta)))
  xp<-mo$platt(clip(predict(mo$xgb,X[testidx,,drop=FALSE])))
  if(length(lp)!=nrow(r) || length(xp)!=nrow(r)) stop('Stage prediction length mismatch: ',s)
  set(r,j=paste0(s,'_LR'),value=lp);set(r,j=paste0(s,'_XGB'),value=xp)
}
for (m in c('LR','XGB')) {
  required_cols <- paste0(c('approach_','authorization_','procurement_'),m)
  if (!all(required_cols %in% names(r))) stop('Missing reachability columns: ',paste(setdiff(required_cols,names(r)),collapse=', '))
  pieces <- lapply(required_cols, function(nm) r[[nm]])
  if (!all(vapply(pieces,length,integer(1)) == nrow(r))) stop('Reachability column length mismatch')
  set(r, j=paste0('reach_',m), value=pieces[[1]]*pieces[[2]]*pieces[[3]])
}
r[,observed_procurement:=as.integer(approached==1 & authorized==1 & procured==1)]
r[,age_group:=cut(age,c(-Inf,39,64,Inf),labels=c('<40','40-64','65+'))]
save_csv(r,'11_test_pipeline_reachability_individual.csv')
reach<-rbindlist(lapply(c('age_group','gender','race'),function(attr){ rr<-copy(r); set(rr,j='group',value=as.character(rr[[attr]])); rr[!is.na(group) & !tolower(group) %in% c('missing','na','n/a','not reported'),.(n=.N,observed_procurement=mean(observed_procurement,na.rm=TRUE),mean_reach_LR=mean(reach_LR),mean_reach_XGB=mean(reach_XGB),product_of_group_means_LR=mean(approach_LR)*mean(authorization_LR)*mean(procurement_LR),product_of_group_means_XGB=mean(approach_XGB)*mean(authorization_XGB)*mean(procurement_XGB)),by=group][,attribute:=attr]}))
save_csv(reach,'12_reachability_vs_observed.csv')
# Subgroup bootstrap: exclude missing demographic categories from comparisons.
# Gender: prespecified Male minus Female contrast (signed).
# Age/race: maximum minus minimum among eligible, sufficiently sized groups.
# Stratified bootstrap preserves subgroup sample sizes and avoids empty resamples.
set.seed(999)
valid_group <- function(z) {
  z <- trimws(as.character(z))
  !is.na(z) & nzchar(z) & !(tolower(z) %in% c('missing','na','n/a','unknown','not reported'))
}
bootstrap_gap <- function(attr, model, B=BOOT_B) {
  x <- as.character(r[[attr]])
  v <- r[[paste0('reach_',model)]]
  ok <- valid_group(x) & is.finite(v)
  x <- x[ok]; v <- v[ok]
  if (attr == 'gender') {
    keep <- x %in% c('F','M')
    x <- x[keep]; v <- v[keep]
    groups <- c('F','M')
  } else {
    sizes <- table(x)
    groups <- names(sizes)[sizes >= MIN_GROUP_N]
    keep <- x %in% groups
    x <- x[keep]; v <- v[keep]
  }
  if (length(groups) < 2L || any(table(factor(x,levels=groups))==0L))
    stop('Insufficient valid groups for ',attr,' / ',model)
  group_values <- lapply(groups,function(g) v[x==g])
  names(group_values) <- groups
  contrast <- function(z) {
    means <- vapply(z,mean,numeric(1))
    if (attr == 'gender') unname(means['M']-means['F']) else diff(range(means))
  }
  estimate <- contrast(group_values)
  bb <- replicate(B, {
    sampled <- lapply(group_values,function(z) sample(z,length(z),replace=TRUE))
    contrast(sampled)
  })
  data.table(attribute=attr,model=model,comparison=if(attr=='gender') 'M minus F' else 'max minus min',
             groups_included=paste(groups,collapse=';'),n_included=length(v),
             estimate=estimate,lo=unname(quantile(bb,.025)),
             hi=unname(quantile(bb,.975)),B=B)
}
boot <- rbindlist(lapply(c('age_group','race','gender'),function(a)
  rbindlist(lapply(c('LR','XGB'),function(m)bootstrap_gap(a,m)))))
save_csv(boot,'13_bootstrap_reachability_gaps.csv')
# PAT: illustrative validation-set threshold adjustment requires separate calibration predictions
# and explicit objective/constraints. It is intentionally NOT optimized on test outcomes here.
# We provide a transparent baseline deterministic pipeline using calibration-derived Youden thresholds.
thr<-rbindlist(threshold_table)
for(m in c('LR','XGB')){
  for(s in names(stage_specs)){
    model_label<-if(m=='LR')'LR' else 'XGB_platt'; t<-thr[stage==s & model==model_label,threshold][1]
    set(r,j=paste0('pass_',s,'_',m),value=as.integer(r[[paste0(s,'_',m)]]>=t))
  }
  set(r,j=paste0('simulated_',m),value=r[[paste0('pass_approach_',m)]]*r[[paste0('pass_authorization_',m)]]*r[[paste0('pass_procurement_',m)]])
}
sim<-rbindlist(lapply(c('age_group','gender','race'),function(attr){ rr<-copy(r); set(rr,j='group',value=as.character(rr[[attr]])); rr[!is.na(group) & !tolower(group) %in% c('missing','na','n/a','not reported'),.(n=.N,observed_procurement=mean(observed_procurement,na.rm=TRUE),simulated_LR=mean(simulated_LR),simulated_XGB=mean(simulated_XGB)),by=group][,attribute:=attr]}))
save_csv(sim,'14_deterministic_simulation.csv')
# Clinical adjustment: compare within available referral-time clinical strata only.
# No causal or 'unjustified discrimination' claims from observational data.
# Simple age-stratified summaries for documented clinical factors, if available.
clinical<-intersect(c('mechanism_of_death','circumstances_of_death'),names(d))
if(length(clinical)){
  ci<-d[partition=='test',c('row_id',clinical),with=FALSE]
  cmerge<-merge(r[,.(row_id,age_group,observed_procurement,reach_LR,reach_XGB)],ci,by='row_id',all.x=TRUE)
  for(v in clinical){a<-cmerge[,.(n=.N,observed=mean(observed_procurement,na.rm=TRUE),reach_LR=mean(reach_LR),reach_XGB=mean(reach_XGB)),by=c('age_group',v)];save_csv(a,paste0('15_clinical_strata_',v,'.csv'))}
}
if(SAVE_MODELS) saveRDS(models,file.path(OUT,'16_fitted_models.rds'))
writeLines(c('Analysis notes:',paste('Input records:',raw_n),paste('Records after contradiction exclusions:',nrow(d)),paste('Split:',SPLIT_MODE),'No final PAT 42% claim is produced: requires a prespecified validation-only optimization and separate evaluation.','IPW is NOT produced: identification and weighting model require separate review.','Product-of-means is descriptive and not an estimator of absolute reachability.','Counterfactual downstream scoring on ineligible referrals is extrapolation; assess overlap and compare with observed transitions.','Check timing of ALL predictors, including cause-of-death and body measurements.',paste('Split mode:',SPLIT_MODE),paste('Feature variant:',VARIANT),'Bootstrap intervals condition on fixed fitted models; they exclude model-refitting uncertainty.','Complete-case analysis (when selected) changes population and may be biased.','Revised results may differ from original manuscript.'),file.path(OUT,'README_interpretation.txt'))
cat('Complete. Outputs in:',normalizePath(OUT),'\n')
  invisible(TRUE)
}
run_status <- data.table(mode=character(),variant=character(),status=character(),message=character())
for(mode in RUN_MODES) {
  for(variant in c(RUN_VARIANTS,if(RUN_COMPLETE_CASE) 'complete_case')) {
    tryCatch({
      run_one(mode,variant)
      run_status <- rbind(run_status,data.table(mode=mode,variant=variant,status='SUCCESS',message=''))
    },error=function(e) {
      run_status <- rbind(run_status,data.table(mode=mode,variant=variant,status='FAILED',message=conditionMessage(e)))
      fwrite(run_status,file.path(ROOT_OUT,'00_run_status.csv'))
      stop(sprintf('V3 stopped at %s / %s: %s',mode,variant,conditionMessage(e)),call.=FALSE)
    })
    fwrite(run_status,file.path(ROOT_OUT,'00_run_status.csv'))
  }
}
writeLines(c('V3.6 analysis completed; compare expanded (code: primary) versus reduced-feature metrics.',
'No IPW or validated PAT effect is claimed. Do not reuse the historical 42% PAT claim.',
'Temporal analysis excludes referral_year; random analysis retains it.',
'Predictor availability is not verified; reduced-feature results are critical.',
'Stagewise and cumulative reachability bootstrap CIs condition on fitted predictions.',
'No leave-one-OPO-out analysis included; requires separate design.'),file.path(ROOT_OUT,'README_MASTER.txt'))
if(nrow(run_status)!=length(RUN_MODES)*(length(RUN_VARIANTS)+as.integer(RUN_COMPLETE_CASE)) || any(run_status$status!='SUCCESS')) stop('Not all V3 configurations completed successfully.')
cat('All V3.6 runs completed successfully in ',normalizePath(ROOT_OUT),'\n')
