/******************************************
SAS Code for DunedinPoAm analyses in NAS

By Xu Gao
2020.4.12

*******************************************/
/*data preperation*/
libname nas "D:\OneDrive - cumc.columbia.edu\Cohorts\nas";
libname mort "D:\OneDrive - cumc.columbia.edu\Cohorts\nas\Mortality_incidence";

ods html;
data workdata;
	set nas.meth_biomarker;
run;

proc sort data = workdata;
	by samplename_450k;
run;

proc sort data = Nas.Xg_07oct19;
	by samplename_450k;
run;

data workdata; /*merge phenotype and methylation data, clean work data*/
	merge Nas.Xg_07oct19 (in=in1) workdata(in=in2);
	by samplename_450k;
	if in1 = in2;
	if neduc =>16 then new_educ = 3;
	else if 15>=neduc =>13 then new_educ = 2;
	else if neduc <=12 then new_educ = 1;
	if . < bmi < 25.0 then bmic = 1;
	else if 25.0 <= 25.0 <30.0 then bmic = 2;
	else if bmi => 30.0 then bmic = 3; 
	p_ca = 0;	if yrca1 <= year and yrca1 ne . then p_ca = 1;
	creat_1 = creat/0.9;
	creat_3 = creat/0.9;
	if race ne '2' then egfr = 141* (min(creat_1, 1)**(-0.411)) * (max(creat_3, 1)**(-1.209)) * (0.993**(age)) ;
	else if race = '2' then egfr = 141* (min(creat_1, 1)**(-0.411)) * (max(creat_3, 1)**(-1.209)) * (0.993**(age)) *1.159;
	ckd = 0 ; 	if .< egfr< 60 then ckd = 1;
	if .< FEV1FVC <70 and FEV_93 < 80 then COPD = 1;
	else COPD = 0;
	cvd = 0;
	if stroke = 1 or CHD = 1 then cvd = 1;
	if ht = . then ht = 0;
	if diabete = . then diabete = 0;
	chro_dis = ht + diabete + cvd + copd + ckd + p_ca;
	if packyrs = . then packyrs = 0;
run;

proc sort data = workdata;
	by id date;
run;

data workdata;
set workdata;
countid_check+1;
by id;
if first.id then countid_check=1; 
run;

/************************************************************************/
/*basic description*/
proc sort data = workdata;
	by countid_check;
run;

proc freq data = workdata;
	tables smk twodrink new_educ ht stroke chd diabete chro_dis countid_check;
run;

proc means data = workdata  mean std median q1 q3 qrange n;
	var age chol trig hdl sbp bmi AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelGrim AgeAccelPheno poa_all;
run;

proc univariate data = workdata ;
	var chro_dis;
	histogram  / midpoints =  1 2 3 4 5 6;
	label chro_dis = "number of chronic diseases";
run;

proc univariate data = workdata ;
	var poa_all;
	histogram ;
run;

proc corr data = workdata spearman nomiss  plots=scatter(nvar=2 alpha=.20 .30);
	var poa_all;
	with age;
run;

proc sgscatter data=workdata;                                                                                    
   plot (age)*(poa_all)                                                                               
        / attrid=myid                                                                                                
          reg=(nogroup degree=1 lineattrs=(color=gray));  
	label poa_all = "mPoA"; 
run;  

proc sql; /*get the max number of disease*/
create table want as
select id,max(chro_dis) as max_dis
from workdata
group by id;
quit;

proc freq data = want;
	tables max_dis;
run;

/*correlation*/
proc corr data = workdata spearman plots=scatter(nvar=2 alpha=.20 .30);
	var AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelPheno AgeAccelGrim ;
	with poa_all;
run;

/*fixed effet model for the longitudinal association between age and mPOA*/
proc sort data = workdata;
	by id;
run;

proc glm data = workdata;
 absorb id;
 model poa_all = age / solution noint; random samplename_450k;  run;
quit;

proc glm data = workdata;
 absorb id;
 model poa_all = age NK CD4t CD8t Bcell Mono gran/ solution noint; random samplename_450k;  run;
quit;

proc glm data = workdata;
 absorb id;
 class smk;
 model poa_all = age smk/ solution noint; random samplename_450k;  run;
quit;

proc glm data = workdata;
 absorb id;
 model poa_all = age packyrs/ solution noint; random samplename_450k;  run;
quit;

/************************************************************************/
/*mortality data*/
proc sort data = workdata;
	by samplename_450k;
run;

proc sort data = Mort.Mort_dnam_2019;
	by samplename_450k;
run;

data workdata0; /*merge the longitudinal data with the updated survival data*/
	merge Mort.Mort_dnam_2019 (in=in1) workdata(in=in2);
	by samplename_450k;
	if in1 = in2;
run;

proc means data = workdata0  mean std; /*distributions of DNA methylation markers*/
	var DNAmAge  DNAmAgeHannum DNAmPhenoAge DNAmGrimAge poa_all AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelGrim AgeAccelPheno;
run;

PROC STANDARD DATA=workdata0 MEAN=0 STD=1 OUT=workdata0; /*z-transformation*/
  VAR poa_all AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelGrim AgeAccelPheno;
RUN;

data workdata0; /*categorization of mPoA*/
	set workdata0;
	if poa_all <= -1 then poa_all_c = 0;
	else if poa_all >= 1 then poa_all_c = 2;
	else poa_all_c = 1;
run;

proc corr data = workdata0 pearson; /*correlation of DNAm Ages with chro age*/
	var DNAmAge  DNAmAgeHannum DNAmPhenoAge;
	with age;
run;

proc freq data = workdata0;
	tables poa_all_c;
run;

/*KM for poa group*/
data workdata_km;
	set workdata0;
	time_fu_yr = time_fu/365.25;
	keep poa_all poa_all_c age id dead time_fu_yr;
run;

proc sort data = workdata_km;
	by poa_all_c;
run;

ods graphics on;
   proc lifetest data=workdata_km  method=lt plots=(s,ls,lls,h,p);
      time time_fu_yr*dead(1);
      strata poa_all_c;
   run;
   ods graphics off;

/*Cox regression*/
%let factor = AgeAccelGrim; /* poa_all AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelGrim AgeAccelPheno*/
%let cov = NK CD4t CD8t Bcell Mono gran;
%let time = time_fu_yr; 
%let case = dead; /*cadmx allcvdmx dead*/
%let workdata = workdata0;

proc freq data = &workdata;
	tables &case; 
run;

proc means data = &workdata min max mean median std sum;
	var &factor &time; 
run;

/*continuous********************/
proc phreg data=&workdata plot=survival; /*basic with age*/
	class plate_450k;
	model &time * &case(0)= &factor age;
	random plate_450k;
	hazardratio &factor / diff=REF CL = both; 
run;

proc phreg data=&workdata plot=survival;
	class plate_450k;
	model &time * &case(0)= &factor &cov age;
	random plate_450k;
	hazardratio &factor / diff=REF CL = both; 
run;

proc phreg data=&workdata plot=survival;
	class  smk plate_450k;
	model &time * &case(0)= &factor age smk ; 
	hazardratio &factor / diff=REF CL = both; 
	random plate_450k;
run;

proc phreg data=&workdata plot=survival; 
	class plate_450k;
	model &time * &case(0)= &factor age packyrs; 
	hazardratio &factor / diff=REF CL = both; 
	random plate_450k;
run;

/*with prevalence chronic disease*/
proc freq data = workdata;
	tables chro_dis  ht diabete  cvd  copd  ckd  p_ca;
run;

PROC STANDARD DATA=workdata MEAN=0 STD=1 OUT=workdata_z; /*z-transformation*/
  VAR poa_all poa_all_185 AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelGrim AgeAccelPheno;
RUN;

%let factor = AgeAccelPheno; /* poa_all AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelGrim AgeAccelPheno*/

proc genmod data = workdata_z;
  class plate_450k/param=glm;
  model chro_dis = &factor age plate_450k/ type3 dist=poisson;
run;

proc genmod data = workdata_z;
  class plate_450k/param=glm;
  model chro_dis = &factor age plate_450k NK CD4t CD8t Bcell Mono gran/ type3 dist=poisson;
run;

proc genmod data = workdata_z;
  class plate_450k smk/param=glm;
  model chro_dis = &factor age plate_450k smk/ type3 dist=poisson;
run;

proc genmod data = workdata_z;
  class plate_450k/param=glm;
  model chro_dis = &factor age plate_450k packyrs/ type3 dist=poisson;
run;

/************************************************************************/
/*incident, multible failure regression model*/
proc sort data = workdata; 
	by samplename_450k;
run;

proc sort data = Mort.inci_chro_dis_2019_meth;
	by samplename_450k;
run;

data workdata_inci;
	merge Mort.inci_chro_dis_2019_meth (in=in1) workdata(in=in2);
	by samplename_450k;
	if in1 = in2;
	time_fu_yr = time_fu/365.25;
run;

PROC STANDARD DATA=workdata_inci MEAN=0 STD=1 OUT=workdata_inci; /*z-transformation*/
  VAR poa_all AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelGrim AgeAccelPheno;
RUN;

%let factor = AgeAccelPheno; /* poa_all AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelGrim AgeAccelPheno*/
%let cov = NK CD4t CD8t Bcell Mono gran;
%let time = time_fu_yr; 
%let case = inci_chro; 
%let workdata = workdata_inci;

proc freq data = &workdata;
	tables &case; 
run;

proc means data = &workdata min max median std;
	var &factor &time; 
run;

/*continuous*/
proc phreg data=&workdata plot=survival; /*basic with age, number of chronic diseases at baseline*/
	class plate_450k;
	model &time * &case(0)= &factor age chro_dis_1;
	random plate_450k;
	hazardratio &factor / diff=REF CL = both; 
run;

proc phreg data=&workdata plot=survival;
	class plate_450k;
	model &time * &case(0)= &factor &cov age chro_dis_1;
	random plate_450k;
	hazardratio &factor / diff=REF CL = both; 
run;

proc phreg data=&workdata plot=survival;
	class  smk plate_450k;
	model &time * &case(0)= &factor age smk chro_dis_1; 
	hazardratio &factor / diff=REF CL = both; 
	random plate_450k;
run;

proc phreg data=&workdata plot=survival; 
	class plate_450k;
	model &time * &case(0)= &factor age packyrs chro_dis_1; 
	hazardratio &factor / diff=REF CL = both; 
	random plate_450k;
run;

/************************************************************************/
/*correlation between visits*/
data tim1;
	set workdata;
	keep id poa_all  AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelPheno;
	where countid_check = 1;
run;

data tim1;
	set tim1;
	rename poa_all= poa_all_1  AgeAccelerationResidual=aahor_1 AgeAccelerationResidualHannum=aahan_1 AgeAccelPheno=aapheno_1;
run;

data tim2;
	set workdata;
	keep id poa_all AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelPheno;
	where countid_check = 2;
run;
data tim2;
	set tim2;
	rename poa_all= poa_all_2  AgeAccelerationResidual=aahor_2 AgeAccelerationResidualHannum=aahan_2 AgeAccelPheno=aapheno_2;
run;
data tim3;
	set workdata;
	keep id poa_all AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelPheno;
	where countid_check = 3;
run;
data tim3;
	set tim3;
	rename poa_all= poa_all_3 AgeAccelerationResidual=aahor_3 AgeAccelerationResidualHannum=aahan_3 AgeAccelPheno=aapheno_3;
run;

data tim4;
	set workdata;
	keep id poa_all AgeAccelerationResidual AgeAccelerationResidualHannum AgeAccelPheno;
	where countid_check = 4;
run;
data tim4;
	set tim4;
	rename poa_all= poa_all_4  AgeAccelerationResidual=aahor_4 AgeAccelerationResidualHannum=aahan_4 AgeAccelPheno=aapheno_4;
run;

data all;	
	merge tim1 tim2 tim3 tim4;
	by id;
run;

proc corr data = all spearman;
	var poa_all_1 poa_all_2 poa_all_3;
run;

proc corr data = all ;
	var aahor_1 aahor_2 aahor_3;
run;

proc corr data = all ;
	var aahan_1 aahan_2 aahan_3;
run;

proc corr data = all ;
	var aapheno_1 aapheno_2 aapheno_3;
run;
