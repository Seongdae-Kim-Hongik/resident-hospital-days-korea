# Part of the public analysis code for resident inpatient days in Korean districts, 2011-2024 (MIT licence). Inputs are read from environment variables; see README.md.
suppressPackageStartupMessages({library(INLA);library(sf);library(spdep);library(jsonlite)})
TAG<-Sys.getenv("TAG");VER<-Sys.getenv("ST_VER",unset="std3");SUF<-paste0("_",VER);EXG<-Sys.getenv("EXG",unset="A")
ROOT<-normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),".."))
SHP<-Sys.getenv("ST_SHP",unset=file.path(ROOT,"00_shp","final.shp"))
OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,paste0("results",SUF));if(nzchar(Sys.getenv("ST_RES_OUT"))){RIN<-RES;RES<-Sys.getenv("ST_RES_OUT");dir.create(RES,showWarnings=FALSE,recursive=TRUE);file.copy(file.path(RIN,"S1_likelihood.csv"),RES)}
LOG<-function(...)cat(sprintf("[88X %s %s] ",EXG,format(Sys.time(),"%H:%M:%S")),sprintf(...),"\n",sep="")
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE)
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);DD$fbd_region<-enc2utf8(DD$fbd_region)
best<-sub("best_set=","",grep("^best_set=",readLines(file.path(OUTD,sprintf("HM_best_%s%s.txt",TAG,SUF))),value=TRUE))
HYP<-as.character(S$sets_after_vif[[best]]);FORCED<-as.character(S$forced);MAINV<-unique(c(HYP,FORCED));MAINFAM<-"nbinomial"
lk<-read.csv(file.path(RES,"S1_likelihood.csv"));MAINFAM<-lk$family[which.min(lk$WAIC)]
G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);n<-nrow(G)
g<-inla.read.graph(Sys.getenv("ST_GRAPH"));stopifnot(g$n==n)
nb<-poly2nb(G,queen=TRUE);for(i in seq_len(g$n))nb[[i]]<-as.integer(sort(g$nbs[[i]]));attr(nb,"sym")<-TRUE;lw<-nb2listw(nb,style="W",zero.policy=TRUE)
DD$area<-match(DD$fbd_region,G$region);stopifnot(!anyNA(DD$area));Y0<-min(DD$year);DD$time<-DD$year-Y0+1;DD$yiid<-DD$time
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001);bt<-function(v)sprintf("`%s`",v)
SPB<-"f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)";SPI<-"f(area,model='besag',graph=g,scale.model=TRUE,constr=TRUE,hyper=pcp)"
TI<-"f(time,model='rw1',scale.model=TRUE,hyper=pcp,values=1:14)";IT<-"f(yiid,model='iid',hyper=pcp,values=1:14)"
fit<-function(vars,d=DD,sp=SPB){form<-as.formula(paste("cases ~",paste(c("1",bt(vars),if(nzchar(sp))sp,TI,IT),collapse="+")))
  tryCatch(inla(form,family=MAINFAM,data=d,E=d$E,control.fixed=CF,control.compute=list(waic=TRUE,dic=TRUE),control.inla=list(control.vb=list(emergency=50))),error=function(e){LOG("fit failed %s",conditionMessage(e));NULL})}
FF<-file.path(RES,sprintf("S60_extra_fit_%s.csv",EXG))
summ<-function(an,m,vars,note){if(is.null(m))return(invisible());r<-m$summary.fixed[vars,,drop=FALSE]
  o<-data.frame(analysis=an,cov=vars,beta=r$mean,sd=r$sd,IRR=exp(r$mean),lo=exp(r$`0.025quant`),hi=exp(r$`0.975quant`),cred=(r$`0.025quant`>0|r$`0.975quant`<0),role=ifelse(vars%in%FORCED,"forced","hypothesis"),stringsAsFactors=FALSE)
  write.csv(o,file.path(RES,sprintf("S60_extra_%s.csv",an)),row.names=FALSE)
  f<-data.frame(analysis=an,note=note,n_obs=sum(!is.na(m$.args$data$cases)),WAIC=round(m$waic$waic,1),n_cred_hyp=sum(o$cred&o$role=="hypothesis"),cred_hyp=paste(o$cov[o$cred&o$role=="hypothesis"],collapse=";"))
  write.table(f,FF,append=file.exists(FF),col.names=!file.exists(FF),row.names=FALSE,sep=",");LOG("%s done WAIC %.1f cred %s",an,m$waic$waic,f$cred_hyp)}
dn<-function(an)file.exists(file.path(RES,sprintf("S60_extra_%s.csv",an)))
zf<-function(x)(x-mean(x))/sd(x)*sqrt(length(x)/(length(x)-1))
zf<-function(x){m<-mean(x);s<-sqrt(mean((x-m)^2));(x-m)/s}
if(EXG=="A"){
  if(!dn("X0_main_fullprec"))summ("X0_main_fullprec",fit(MAINV),MAINV,"main model refit, full-precision beta/sd (numbers json source; must match S3 to 3 decimals)")
  if(Sys.getenv("ST_X_ONLY")=="X0"){LOG("ST_X_ONLY=X0 → X0만");quit(save="no",status=0)}
  if(!dn("X1_icar"))summ("X1_icar",fit(MAINV,sp=SPI),MAINV,"ICAR(besag) instead of BYM2")
  if(!dn("X2_nospatial"))summ("X2_nospatial",fit(MAINV,sp=""),MAINV,"no spatial random effect (RW1+iid year only)")
  if(!dn("X3_twoway")){d<-DD;nv<-character(0);for(v in HYP){xi<-ave(d[[v]],d$area);xt<-ave(d[[v]],d$year);d[[paste0(v,"__within")]]<-d[[v]]-xi-xt+mean(d[[v]]);d[[paste0(v,"__between")]]<-xi;nv<-c(nv,paste0(v,"__within"),paste0(v,"__between"))}
    summ("X3_twoway",fit(c(nv,FORCED),d),c(nv,FORCED),"two-way centring: within (x_it - x_i. - x_.t + x..) and between (x_i.) components, per 1 SD of pooled x")}
  if(!dn("X4_area_density")){d<-DD;km2<-d$area_m2/1e6;vv<-MAINV
    for(v in intersect(c("public_sports_facilities_per_100k","village_sports_facilities_per_100k"),vv)){nv<-paste0(v,"__per_km2");d[[nv]]<-zf(d[[paste0(v,"_raw")]]*d$population/1e5/km2);vv[vv==v]<-nv}
    summ("X4_area_density",fit(vv,d),vv,"public and village sports facilities per km2 (count/area, z) instead of per population")}
  if(!dn("X5_boundary_excl")){d<-DD;bx<-grepl("군위|세종|청주|부천",d$fbd_region);d$cases[bx]<-NA;summ("X5_boundary_excl",fit(MAINV,d),MAINV,sprintf("boundary-change districts excluded (Gunwi, Sejong, Cheongju, Bucheon; %d cells)",sum(bx)))}
  if(!dn("X6_interp_year_excl")){d<-DD;bx<-d$year%in%c(2015,2018,2019);d$cases[bx]<-NA;summ("X6_interp_year_excl",fit(MAINV,d),MAINV,sprintf("years with interpolated covariates excluded (2015 public sports facilities, 2018 moderate activity, 2019 obesity; %d cells)",sum(bx)))}
  if(!dn("X9_urban0_excl")){d<-DD;z0<-unique(d$fbd_region[d$urban_population_share_raw==0]);bx<-d$fbd_region%in%z0;d$cases[bx]<-NA;summ("X9_urban0_excl",fit(MAINV,d),MAINV,sprintf("districts with zero urban-area population excluded (%s; %d cells)",paste(z0,collapse=";"),sum(bx)))}
  FP1<-file.path(RES,"S60_P1_moran_constant.csv");if(!file.exists(FP1)){set.seed(2);o<-NULL
    for(y in sort(unique(DD$year))){d<-DD[DD$year==y,];d<-d[match(G$region,d$fbd_region),];a<-moran.mc(log(d$cases/d$E),lw,999,zero.policy=TRUE,alternative="greater");set.seed(2);b<-moran.mc(log((d$cases+0.5)/d$E),lw,999,zero.policy=TRUE,alternative="greater")
      o<-rbind(o,data.frame(year=y,I_no_constant=unname(a$statistic),p_no_constant=a$p.value,I_plus05=unname(b$statistic),p_plus05=b$p.value,dI=unname(b$statistic-a$statistic),zero_cells=sum(d$cases==0)))}
    write.csv(o,FP1,row.names=FALSE);LOG("P1 max |dI| %.5f",max(abs(o$dI)))}
}
if(EXG%in%c("B","C")){
  VC<-read.csv(file.path(OUTD,"gate","tables","G00_village_changed_cells.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);VC$fbd_region<-enc2utf8(VC$fbd_region)
  FL<-list(public_sports_facilities_per_100k=DD$year==2015,moderate_pa_rate=DD$year==2018,obesity_rate_sr=DD$year==2019,village_sports_facilities_per_100k=paste(DD$fbd_region,DD$year)%in%paste(VC$fbd_region,VC$year))
  FL<-FL[names(FL)%in%HYP];stopifnot(length(FL)>=3)
  RR<-if(EXG=="B")1:10 else 11:20
  for(r_ in RR){an<-sprintf("X7_mi_r%02d",r_);if(dn(an))next;set.seed(20261126+r_);d<-DD
    for(v in names(FL)){raw<-d[[paste0(v,"_raw")]];f<-FL[[v]];sdv<-ave(ifelse(f,NA,raw),d$area,FUN=function(s)sd(s,na.rm=TRUE));sdv[is.na(sdv)]<-sd(raw)*0.5
      raw[f]<-pmax(raw[f]+rnorm(sum(f),0,sdv[f]),0);d[[v]]<-zf(raw)}
    summ(an,fit(MAINV,d),MAINV,sprintf("multiple imputation replicate %d: interpolated/corrected cells perturbed (%s)",r_,paste(sprintf("%s %d",names(FL),sapply(FL,sum)),collapse=", ")))}
}
LOG("EXG %s ALL DONE",EXG)
