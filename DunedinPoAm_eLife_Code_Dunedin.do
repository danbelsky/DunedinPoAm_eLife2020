//****************************************************************//
//****************************************************************//
//DUNEDIN
//****************************************************************//
//****************************************************************//
global user db3275
use "/Users/$user/Box/Belsky/mPoA/Documentation/CodeToPrepareData/mPoA_DUNEDINpheno.dta", clear
merge 1:1 snum using "/Users/$user/Box/Belsky/mPoA/Documentation/CodeToPrepareData/mPoA_DUNEDINmeth.dta",nogen
merge 1:1 snum using "/Users/$user/Box/Belsky/mPoA/mPOA_revisionsFeb21.dta"
merge 1:1 snum using "/Users/db3275/Box/Belsky/PNAS2015/Documentation/1PaceofAging File.dta", keepusing(Big6craep38 Big6crvep38 VenularRsd38 ArteriolRsd38) nogen

//Grooved Pegboard Time
gen gpdom38=gpdomt38*-1
gen gpdom45=GpDomTim45*-1
gen d_gpdom=(gpdom45-gpdom38)

quietly sum gpdom38
gen zd_gpdom=(d_gpdom/r(sd))
mean gpdom* d_gpdom
//Retinal Vessel Phenotypes
	//Alternate Standardization to show change
foreach v in 38 45{
	capture drop A`v'*
	gen A`v'=.
	foreach s in 1 2{
		quietly sum Big6crvep`v' if ArteriolRsd`v'!=. & sex==`s'
		local S=r(mean)
		reg Big6craep`v' Big6crvep`v' if ArteriolRsd`v'!=.  & sex==`s'
		predict A`v'`s' if e(sample), r
		replace A`v'`s' = A`v'`s'+`S'*_b[Big6crvep`v'] + _b[_cons]
		replace A`v'=A`v'`s' if sex==`s'
		}
	drop A`v'1 A`v'2
	
	capture drop V`v'*
	gen V`v'=.
	foreach s in 1 2{
		quietly sum Big6craep`v' if VenularRsd`v'!=. & sex==`s'
		local S=r(mean)
		reg Big6crvep`v' Big6craep`v' if VenularRsd`v'!=.  & sex==`s'
		predict V`v'`s' if e(sample), r
		replace V`v'`s' = V`v'`s'+`S'*_b[Big6craep`v'] + _b[_cons]
		replace V`v'=V`v'`s' if sex==`s'
		}
	drop V`v'1 V`v'2	
	}
foreach x in A V{
	capture drop D_`x'
	gen D_`x' = `x'45-`x'38
	capture drop zD_`x'
	quietly sum `x'38
	gen zD_`x'= D_`x'/r(sd)	
	}
	//Confirm new standardized variable matches distribution of the raw variable
sum   Big6craep38 A38  Big6craep45 A45 Big6crvep38 V38 Big6crvep45 V45 if A38 !=. & A45!=.
	//Confirm new standardized variable correlates perfectly with HL's standardized variable
corr A38 ArteriolRsd38 Big6craep38

//HL Variables
gen d_Art=ArteriolRsd45-ArteriolRsd38
gen d_Ven=VenularRsd45-VenularRsd38
quietly sum ArteriolRsd38
gen zd_Art=d_Art/r(sd)
quietly sum VenularRsd38
gen zd_Ven=d_Ven/r(sd)
mean ArteriolRsd* d_Art VenularRsd*	d_Ven
corr mpoa ArteriolRsd45 d_Art VenularRsd45 d_Ven

#delimit ;
graph box d_Art zD_A d_Ven zD_V, 
	box(1, color(blue)) marker(1, mcolor(blue))
	box(2, color(red)) marker(2, mcolor(red)) 
	box(3, color(dknavy)) marker(3, mcolor(dknavy)) 
	box(4, color(cranberry)) marker(4, mcolor(cranberry))
	 ; #delimit cr

save "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_Dunedin.dta", replace 	
	
use "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_Dunedin.dta", clear	
	
putexcel set "/Users/$user/Box/Belsky/mPoA/Documentation/mPoATables_R1.xlsx", sheet(Dunedin45) modify
putexcel A1 = "Dunedin Age 45 Analysis"	
count if mpoa!=.
putexcel A3 = "mPoA N="
putexcel B3 = matrix(r(N))
count if paceofaging!=.	
putexcel A4 = "Pace of Aging N="
putexcel B4 = matrix(r(N))
count if paceofaging!=.	& mpoa!=.
putexcel A5 = "mPoA & Pace of Aging N="
putexcel B5 = matrix(r(N))

/*ArteriolRsd45 */
/*VenularRsd45	*/

#delimit ;
global PHENO "
	balClsMax45 
	Velocity_m45 
	StepPlace45 
	ChairStands45 
	GripMax45
	gpdom45
	PhyLimts45 
	pri45
	wmi45
	psi45
	A45 
	V45 
	Health45  
	ZFacialAge45 
	" ; #delimit cr
	
	//Summary Statistics
tabstat $PHENO, s(mean sd n) by(sex) save
matrix A1 = r(StatTotal)' , r(Stat1)', r(Stat2)'
matrix list A1
	//Regression Analysis
foreach y in mpoa paceofaging{
	capture drop z`y'
	egen z`y' =std(`y')
	matrix Fx_`y' = J(1,5,999)
	matrix Fx_`y'_cells = J(1,5,999)
	matrix Fx_`y'_smk = J(1,5,999)
	foreach x in $PHENO{
		capture drDop Z 
		capture drop Z 
		gen Z =.
		foreach s in 1 2{
			quietly sum `x' if sex==`s'
			replace Z = (`x'-r(mean))/r(sd) if sex==`s'
			}
		//Base Model
		quietly reg Z z`y' sex, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix Fx_`y' = Fx_`y' \ A
		//Adjusted for Cells
		quietly reg Z z`y' sex plasmablast38 cd8pcd28ncd45ran38 cd8naive38 cd4naive38 cd8t38 cd4t38 nk38 bcell38 mono38 gran38, robust 		
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)	
		matrix rownames A = `x'
		matrix Fx_`y'_cells = Fx_`y'_cells \ A
		//Adjusted for Smoking (pack years)
		quietly reg Z z`y' sex PackYrLifTm38, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix Fx_`y'_smk = Fx_`y'_smk \ A
		drop Z
		}
		drop z`y'
	foreach q in "" _cells _smk{ 
		matrix Fx_`y'`q' = Fx_`y'`q'[2...,1...]
		matrix colnames Fx_`y'`q' = r lb ub p N
		matrix list Fx_`y'`q' 
		}
	}

//Analysis of cases w/ mPoA & Pace of Aging
keep if mpoa!=. & paceofaging!=.
	//Summary Statistics
tabstat $PHENO, s(mean sd n) by(sex) save
matrix RA1 = r(StatTotal)' , r(Stat1)', r(Stat2)'
matrix list RA1
	//Regression Analysis
foreach y in mpoa paceofaging{
	capture drop z`y'
	egen z`y' =std(`y')
	matrix RFx_`y' = J(1,5,999)
	matrix RFx_`y'_cells = J(1,5,999)
	matrix RFx_`y'_smk = J(1,5,999)
	foreach x in $PHENO{
		capture drop Z 
		gen Z =.
		foreach s in 1 2{
			quietly sum `x' if sex==`s'
			replace Z = (`x'-r(mean))/r(sd) if sex==`s'
			}
		//Base Model
		quietly reg Z z`y' sex, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix RFx_`y' = RFx_`y' \ A
		//Adjusted for Cells
		quietly reg Z z`y' sex plasmablast38 cd8pcd28ncd45ran38 cd8naive38 cd4naive38 cd8t38 cd4t38 nk38 bcell38 mono38 gran38, robust 		
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)	
		matrix rownames A = `x'
		matrix RFx_`y'_cells = RFx_`y'_cells \ A
		//Adjusted for Smoking (pack years)
		quietly reg Z z`y' sex PackYrLifTm38, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix RFx_`y'_smk = RFx_`y'_smk \ A
		drop Z
		}
		drop z`y'
	foreach q in "" _cells _smk{ 
		matrix RFx_`y'`q' = RFx_`y'`q'[2...,1...]
		matrix colnames RFx_`y'`q' = r lb ub p N
		matrix list RFx_`y'`q' 
		}
	}


//Post Tables to Excel Sheet
	//Tables from Analysis of All Available Data
putexcel A8 = "Analysis of All Available Data"	
putexcel B9 = "mPoA Effect-Sizes"
putexcel H9 = "Pace of Aging Effect-Sizes"
putexcel A27 = "Models Adjusted for Estimated Cell Counts at Age 38"	
putexcel A47 = "Models Adjusted for Pack Years at Age 38"	
putexcel A10 = matrix(Fx_mpoa), names
putexcel H10 = matrix(Fx_paceofaging), colnames	
putexcel A28 = matrix(Fx_mpoa_cells), names
putexcel H28 = matrix(Fx_paceofaging_cells), colnames
putexcel A48 = matrix(Fx_mpoa_smk), names
putexcel H48 = matrix(Fx_paceofaging_smk), colnames
	
		//Tables from Analysis of Restricted Dataset with both mPoA and Pace of Aging
putexcel A65 = "Analysis of SMs w/ Data on mPoA and Pace of Aging"	
putexcel A66 = "mPoA Effect-Sizes"
putexcel H66 = "Pace of Aging Effect-Sizes"
putexcel A84 = "Models Adjusted for Estimated Cell Counts at Age 38"
putexcel A104 = "Models Adjusted for Pack Years at Age 38"	
putexcel A67 = matrix(RFx_mpoa), names
putexcel H67 = matrix(RFx_paceofaging), colnames	
putexcel A85 = matrix(RFx_mpoa_cells), names
putexcel H85 = matrix(RFx_paceofaging_cells), colnames
putexcel A105 = matrix(RFx_mpoa_smk), names
putexcel H105 = matrix(RFx_paceofaging_smk), colnames	


//Summary Statistics Tables 
	//All Available Data
putexcel A123 = "Summary Statistics - All Available Data"
putexcel B124 = "All SMs"
putexcel E124 = "Women"
putexcel H124 = "Men"
putexcel A125 = matrix(A1), names

putexcel A141 = "Summary Statistics - SMs with data on both mPoA and Pace of Aging"
putexcel B142 = "All SMs"
putexcel E142 = "Women"
putexcel H142 = "Men"
putexcel A143 = matrix(RA1), names


//****************************************************************************************//
//****************************************************************************************//
// CHANGE SCORE ANALYSIS 38-45
//****************************************************************************************//
//****************************************************************************************//
putexcel set "/Users/$user/Box/Belsky/mPoA/Documentation/mPoATables_R1.xlsx", sheet(Dunedin3845) modify

count if mpoa!=. & (balance38!=. & balance45!=. | grip38!=. & grip45!=. | limitations38!=. & limitations 45!=. | srh38!=. & srh45!=. | face38!=. & face45!=.)
putexcel A3 = "mPoA N="
putexcel B3 = matrix(r(N))
count if paceofaging!=.	& (balance38!=. & balance45!=. | grip38!=. & grip45!=. | limitations38!=. & limitations 45!=. | srh38!=. & srh45!=. | face38!=. & face45!=.)
putexcel A4 = "Pace of Aging N="
putexcel B4 = matrix(r(N))
count if paceofaging!=.	& mpoa!=. & (balance38!=. & balance45!=. | grip38!=. & grip45!=. | limitations38!=. & limitations 45!=. | srh38!=. & srh45!=. | face38!=. & face45!=.)
putexcel A5 = "mPoA & Pace of Aging N="
putexcel B5 = matrix(r(N))

//****************************************************************************************//
//Summary Statistics
//****************************************************************************************//
global DPHENO "DB DG d_gpdom DL Diq D_A D_V iPoor DF" //"DB DG d_gpdom DL Diq d_Art d_Ven iPoor DF"
	//Change Score Summary Statistics
mean $DPHENO if mpoa !=.	
tabstat $DPHENO, s(mean sd n) by(sex) save
matrix A = r(StatTotal)' , r(Stat1)', r(Stat2)'
putexcel E8 = "Women"
putexcel H8 = "Men"
putexcel A9 = "Change Score Summary Statistics - All Available Data"
putexcel A10 = matrix(A), names
	//Summary Statistics for SMs w mpoa
tabstat $DPHENO if mpoa!=., s(mean sd n) by(sex) save
matrix A = r(StatTotal)' , r(Stat1)', r(Stat2)'
putexcel A21 = "Change Score Summary Statistics - SMs w/ mPoA Data"
putexcel A22 = matrix(A), names
	//Summary Statistics for SMs w Pace of Aging
tabstat $DPHENO if paceofaging!=., s(mean sd n) by(sex) save
matrix A = r(StatTotal)' , r(Stat1)', r(Stat2)'
putexcel A33 = "Change Score Summary Statistics - SMs w/ Pace of Aging Data"
putexcel A34 = matrix(A), names

	//Means & 95% CIs of change in terms of baseline SDs for in-text reporting
matrix Fx = J(1,4,999)
foreach y in zDB zDG zd_gpdom zDL Diq zD_A zD_V DF { // zDB zDG zd_gpdom zDL Diq zd_Art zd_Ven DF
	mean `y' if mpoa!=.
	matrix A = _b[`y'] , _b[`y'] - invttail(e(df_r),0.025)*_se[`y'], _b[`y'] + invttail(e(df_r),0.025)*_se[`y'], e(N)
	matrix rownames A = `y'
	matrix Fx = Fx \ A
	}
matrix Fx =Fx[2...,1...]
matrix colnames Fx = M lb ub N
matrix list Fx
putexcel A83 = "Change in Terms of Baseline SD"
putexcel A84 = matrix(Fx), names

	//% with Fair or Poor SRH at 38 and 45 for in-text reporting
tab srh45 if mpoa!=. & DSRH!=.
tab srh38 if mpoa!=. & DSRH!=.

	//Regression Analysis
foreach y in mpoa paceofaging{
	capture drop z`y'
	egen z`y' =std(`y')
	matrix DFx_`y' = J(1,5,999)
	matrix DFx_`y'_cells = J(1,5,999)
	matrix DFx_`y'_smk = J(1,5,999)
	foreach x in $DPHENO{
		if `"`x'"'=="iPoor"{
				//Base Model
			quietly poisson `x' z`y' sex, robust 
			matrix A = exp(_b[z`y']) ,  exp(_b[z`y'] - invnormal(0.975)*_se[z`y']), exp(_b[z`y'] + invnormal(0.975)*_se[z`y']), 2*normal(-abs(_b[z`y']/_se[z`y'])), e(N)
			matrix rownames A = `x'
			matrix DFx_`y' = DFx_`y' \ A			
				//Adjusted for Cells
			quietly poisson `x' z`y' sex plasmablast38 cd8pcd28ncd45ran38 cd8naive38 cd4naive38 cd8t38 cd4t38 nk38 bcell38 mono38 gran38, robust 
			matrix A = exp(_b[z`y']) ,  exp(_b[z`y'] - invnormal(0.975)*_se[z`y']), exp(_b[z`y'] + invnormal(0.975)*_se[z`y']), 2*normal(-abs(_b[z`y']/_se[z`y'])), e(N)
			matrix rownames A = `x'
			matrix DFx_`y'_cells = DFx_`y'_cells \ A
				//Adjusted for Smoking
			quietly poisson `x' z`y' sex PackYrLifTm38, robust 
			matrix A = exp(_b[z`y']) ,  exp(_b[z`y'] - invnormal(0.975)*_se[z`y']), exp(_b[z`y'] + invnormal(0.975)*_se[z`y']), 2*normal(-abs(_b[z`y']/_se[z`y'])), e(N)
			matrix rownames A = `x'
			matrix DFx_`y'_smk = DFx_`y'_smk \ A
			}	
		else{
		capture drop Z 
		gen Z =.
		foreach s in 1 2{
			quietly sum `x' if sex==`s'
			replace Z = (`x'-r(mean))/r(sd) if sex==`s'
			}
		//Base Model
		quietly reg Z z`y' sex, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix DFx_`y' = DFx_`y' \ A
		//Adjusted for Cells
		quietly reg Z z`y' sex plasmablast38 cd8pcd28ncd45ran38 cd8naive38 cd4naive38 cd8t38 cd4t38 nk38 bcell38 mono38 gran38, robust 		
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)	
		matrix rownames A = `x'
		matrix DFx_`y'_cells = DFx_`y'_cells \ A
		//Adjusted for Smoking (pack years)
		quietly reg Z z`y' sex PackYrLifTm38, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix DFx_`y'_smk = DFx_`y'_smk \ A
		drop Z
		}
		}
		drop z`y'
	foreach q in "" _cells _smk{ 
		matrix DFx_`y'`q' = DFx_`y'`q'[2...,1...]
		matrix colnames DFx_`y'`q' = r lb ub p N
		matrix list DFx_`y'`q' 
		}
	}

//Post Tables to Excel
	putexcel A45 = "mPoA Regression Analysis of Change Scores"
	putexcel A46 = matrix(DFx_mpoa), names
	putexcel H46 = matrix(DFx_paceofaging), colnames
	
	putexcel A57 = "mPoA Regression Analysis of Change Scores - Adjustd for Estimated Cell Counts at Age 38"
	putexcel A58 = matrix(DFx_mpoa_cells), names
	putexcel H58 = matrix(DFx_paceofaging_cells), colnames	

	putexcel A69 = "mPoA Regression Analysis of Change Scores - Adjustd for Smoking at Age 38 (pack years)"	
	putexcel A70 = matrix(DFx_mpoa_smk), names
	putexcel H70 = matrix(DFx_paceofaging_smk), colnames	
	
//Sensitivity Analysis - in-text reporting of residualized change analysis of IQ
preserve
	capture drop zmpoa
	egen zmpoa=std(mpoa)
	egen ziq45=std(fsiq45a)
	egen zDiq=std(Diq)
	reg zDiq zmpoa sex, robust  
	reg zDiq zmpoa wfsiq713 sex, robust 
restore

 	

//****************************************************************************************//
//****************************************************************************************//	
	

	
	
//****************************************************************************************//
//****************************************************************************************//
// COMPARISON WITH EPIGENETIC CLOCKS -- see code at end of this document beginning line 647
//****************************************************************************************//
//****************************************************************************************//	
do "/Users/$user/Box/Belsky/mPoA/eLife/R1/Analysis_Dunedin_R1_ClockComparison.do"	
	
//****************************************************************************************//
//****************************************************************************************//	
	


//FIGURES	
global user db3275
use "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_Dunedin.dta", clear	

//******************************************************************************************//
//SF1A. mPoA by Pace of Aging
//******************************************************************************************//
#delimit ;
scatter mpoa paceofaging, mfcolor(blue%30) mlcolor(blue%5) 
	|| lfit mpoa paceofaging, lcolor(cranberry) lwidth(thick) lpattern(dash)
	scheme(s1mono)
	legend(off)
	ylabel(, angle(horiz) nogrid labsize(large) format(%9.2f)) 
	xlabel(, angle(horiz) nogrid labsize(large))
	xtitle(12 year Longitudinal Pace of Aging, size(large) margin(small))
	ytitle(DunedinPoAm at age 38, size(medlarge) margin(small))
	name(mpoabypoa,replace)
	; #delimit cr
//******************************************************************************************//
//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/SuppFigure1A.pdf", replace	
//******************************************************************************************//


//******************************************************************************************//
//SF1B. Comparison of mPoA & Pace of Aging Effect-Sizes
//******************************************************************************************//
//Effect-size Comparison - SMs with data on both mPoA and Pace of Aging
clear
svmat2 RFx_paceofaging, names(col) rnames(DV)
gen X = 1 
gen Pheno=_n
save temp, replace
clear
svmat2 RFx_mpoa, names(col) rnames(DV)
gen X = 2 
gen Pheno=_n
append using temp
save temp, replace 

use temp, clear
drop if inlist(Pheno,11,12)
recode Pheno (13=11) (14=12) 
capture label drop Pheno
#delimit ;
	label define Pheno 
		1 "Balance"
		2 "Gait Speed"
		3 "Steps in Place"
		4 "Chair Stands"
		5 "Grip Strength"
		6 "Motor Coordination"
		7 "Physical Limitations"
		8 "Perceptual Reasoning"
		9 "Working Memory"
		10 "Proceessing Speed"
		11 "Self-Rated Health"
		12 "Facial Aging"
		; #delimit cr
label values Pheno Pheno
preserve 
	foreach x in r ub lb{
	replace `x' = `x'*-1 if inlist(Pheno,7,12)
	replace `x' = `x'*-1
	}
#delimit ;	
		graph hbar r, over(X) over(Pheno, label(valuelabel angle(horiz)) ) 
		asyvars bar(2,color(dknavy)) bar(1,color(midblue) lwidth(medthick))
		plotregion(color(white)) graphregion(color(white)) note("")
		legend(ring(0) pos(3) order(1 2) cols(1) lab(1 "Pace of Aging") lab(2 "DunedinPoAm") region(lcolor(white)) symxsize(5))
		ylabel(, angle(horiz) nogrid)
		yline(0, lcolor(gs10))
		yscale(range(0 .5))
		ytitle(Effect-Size (Pearson r equivalent))
		title("")	
		name(compfx, replace)
		; #delimit cr		
restore
//******************************************************************************************//
//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/SuppFigure1B.pdf", replace	
//******************************************************************************************//


//******************************************************************************************//
//Figure1A. mPoA Effect-Sizes for Age 45 Levels
//******************************************************************************************//
//Create Dataset of Effect-Sizes		
clear
svmat2 Fx_mpoa, names(col) rnames(DV)
gen X = 2 
gen Pheno=_n
save temp, replace 
drop if inlist(Pheno,11,12)
recode Pheno (13=11) (14=12) 

capture label drop Pheno
#delimit ;
	label define Pheno 
		1 "Balance"
		2 "Gait Speed"
		3 "Steps in Place"
		4 "Chair Stands"
		5 "Grip Strength"
		6 "Motor Coordination"
		7 "Physical Limitations"
		8 "Perceptual Reasoning"
		9 "Working Memory"
		10 "Proceessing Speed"
		11 "Self-Rated Health"
		12 "Facial Aging"
		; #delimit cr
label values Pheno Pheno

preserve 
foreach x in r ub lb{
replace `x' = `x'*-1 if inlist(Pheno,7,12)
replace `x' = `x'*-1
}
#delimit ; 
twoway bar r Pheno if X==2 & inlist(Pheno,1,2,3,4,5,6,7), horiz color(dknavy)
	|| bar r Pheno if X==2 & inlist(Pheno,8,9,10), horiz color(midblue)
	|| bar r Pheno if X==2 & inlist(Pheno,11,12), horiz color(orange)
	|| rcap lb ub Pheno if X==2, horiz lwidth(medthick) lcolor(gs12)
	yscale(reverse)
	ylabel(1(1)12,valuelabels angle(horiz) labsize(large) nogrid)
	ytitle("")
	xlabel(,labsize(large))
	xtitle(Effect-size (Pearson r), size(large))
	graphregion(color(white))
	plotregion(color(white))
	name(age45fx, replace)	
	legend(off)
	; #delimit cr
restore
//******************************************************************************************//
//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/Figure1A.pdf", replace	
//******************************************************************************************//


	
//******************************************************************************************//
//Figure1B. mPoA Effect-Sizes for Change Scores
//******************************************************************************************//
global user db3275
use "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_Dunedin.dta", clear	

global ytitlesize "medium"
egen zmpoa = std(mpoa)
#delimit ;
binscatter zDB zmpoa , controls(sex) nq(50) lcolor(red) mcolor(dknavy)
	ylabel(,angle(horiz) labsize(medlarge) format(%9.1f) nogrid)
	ytitle("Balance 38-45", size($ytitlesize))
	xtitle(Age 38 DunedinPoAm (z-score units))
	graphregion(color(white))
	plotregion(color(white))
	name(balance, replace)
	; #delimit cr

#delimit ;
binscatter zDG zmpoa , controls(sex) nq(50) lcolor(red) mcolor(dknavy)
	ylabel(,angle(horiz) labsize(medlarge) format(%9.1f) nogrid)
	ytitle("Grip Strength 38-45" , size($ytitlesize))
	xtitle(Age 38 DunedinPoAm (z-score units))
	graphregion(color(white))
	plotregion(color(white))
	name(grip, replace)
	; #delimit cr
	
#delimit ;
binscatter zd_gpdom zmpoa , controls(sex) nq(50) lcolor(red) mcolor(dknavy)
	ylabel(,angle(horiz) labsize(medlarge) format(%9.1f) nogrid)
	ytitle("Motor Coordination 38-45" , size($ytitlesize))
	xtitle(Age 38 DunedinPoAm (z-score units))
	graphregion(color(white))
	plotregion(color(white))
	name(motor, replace)
	; #delimit cr	

#delimit ;
binscatter zDL zmpoa , controls(sex) nq(50) lcolor(red) mcolor(dknavy)
	ylabel(,angle(horiz) labsize(medlarge) format(%9.1f) nogrid)
	ytitle("Limitations 38-45", size($ytitlesize))
	xtitle(Age 38 DunedinPoAm (z-score units))
	graphregion(color(white))
	plotregion(color(white))
	name(limitations, replace)
	; #delimit cr	
	
#delimit ;
binscatter Diq zmpoa , controls(sex) nq(50) lcolor(red) mcolor(midblue)
	ylabel(,angle(horiz) labsize(medlarge) format(%9.0f) nogrid)
	ytitle("Cognitive Decline 13-45 (IQ points)" , size($ytitlesize))
	xtitle(Age 38 DunedinPoAm (z-score units))
	graphregion(color(white))
	plotregion(color(white))
	name(diq, replace)
	; #delimit cr
	
#delimit ;
binscatter DF zmpoa , controls(sex) nq(50) lcolor(cranberry) mcolor(orange)
	ylabel(,angle(horiz) labsize(medlarge) format(%9.1f) nogrid)
	ytitle("Facial Aging 38-45" , size($ytitlesize))
	xtitle(Age 38 DunedinPoAm (z-score units))
	graphregion(color(white))
	plotregion(color(white))
	name(face, replace)
	; #delimit cr	
	
logit iPoor zmpoa sex
margins,at(zmpoa=(-2(.25)2))
#delimit ;
marginsplot, recastci(rarea) 
	plot1opts(lcolor(cranberry) mcolor(none)) 
	ci1opts(lcolor(white) fcolor(orange*.5))
	ylabel(,angle(horiz) labsize(medlarge) format(%9.2f) nogrid)
	xlabel(-2(1)2, labsize(medlarge))
	ytitle("Incident Fair/Poor Health 38-45", size($ytitlesize))
	xtitle(Age 38 DunedinPoAm (z-score units))
	graphregion(color(white))
	plotregion(color(white))
	title("")
	name(srh, replace)
	; #delimit cr
	
//******************************************************************************************//
//Produce combined graph
#delimit ;
graph combine balance grip motor limitations diq srh face , cols(2) holes(6) scheme(s1mono) name(F1B,replace) xsize(4) ysize(7)
; #delimit cr
	//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/Figure1B.pdf", replace	
//******************************************************************************************//





//****************************************************************************************//
//DUNEDIN COMPARISON WITH EPIGENETIC CLOCKS
//****************************************************************************************//
global user db3275
use "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_Dunedin.dta", clear	
	
gen Horvath=horvathage38-ExactAge38 
gen Hannum=hannumage38-ExactAge38
gen Levine=LevineAge38-ExactAge38

//****************************************************************************************//
//****************************************************************************************//
#delimit ;
global PHENO "
	balClsMax45 
	Velocity_m45 
	StepPlace45 
	ChairStands45 
	GripMax45
	gpdom45
	PhyLimts45 
	pri45
	wmi45
	psi45
	A45 
	V45 
	Health45  
	ZFacialAge45 
	" ; #delimit cr
	
//Regression Analysis
foreach y in mpoa Horvath Hannum Levine {
	capture drop z`y'
	egen z`y' =std(`y')
	matrix Fx_`y' = J(1,5,999)
	matrix Fx_`y'_cells = J(1,5,999)
	matrix Fx_`y'_smk = J(1,5,999)
	foreach x in $PHENO{
		capture drop Z 
		gen Z =.
		foreach s in 1 2{
			quietly sum `x' if sex==`s'
			replace Z = (`x'-r(mean))/r(sd) if sex==`s'
			}
		//Base Model
		quietly reg Z z`y' sex, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix Fx_`y' = Fx_`y' \ A
		//Adjusted for Cells
		quietly reg Z z`y' sex plasmablast38 cd8pcd28ncd45ran38 cd8naive38 cd4naive38 cd8t38 cd4t38 nk38 bcell38 mono38 gran38, robust 		
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)	
		matrix rownames A = `x'
		matrix Fx_`y'_cells = Fx_`y'_cells \ A
		//Adjusted for Smoking (pack years)
		quietly reg Z z`y' sex PackYrLifTm38, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix Fx_`y'_smk = Fx_`y'_smk \ A
		drop Z
		}
		drop z`y'
	foreach q in "" _cells _smk{ 
		matrix Fx_`y'`q' = Fx_`y'`q'[2...,1...]
		matrix colnames Fx_`y'`q' = r lb ub p N
		matrix list Fx_`y'`q' 
		}
	}
//****************************************************************************************//	
	
//****************************************************************************************//	
//Post Tables to Excel Sheet
//****************************************************************************************//
putexcel set "/Users/$user/Box/Belsky/mPoA/Documentation/mPoATables_R1.xlsx", sheet(D45Comp) modify

putexcel B9 = "mPoA Effect-Sizes"
putexcel H9 = "Horvath Clock Effect-Sizes"
putexcel N9 = "Hannum Clock Effect-Sizes"
putexcel T9 = "Levine Clock Effect-Sizes"

putexcel A27 = "Models Adjusted for Estimated Cell Counts at Age 38"	
putexcel A47 = "Models Adjusted for Pack Years at Age 38"	
putexcel A10 = matrix(Fx_mpoa), names
putexcel H10 = matrix(Fx_Horvath), colnames	
putexcel N10 = matrix(Fx_Hannum), colnames	
putexcel T10 = matrix(Fx_Levine), colnames	

putexcel A28 = matrix(Fx_mpoa_cells), names
putexcel H28 = matrix(Fx_Horvath_cells), colnames
putexcel N28 = matrix(Fx_Hannum_cells), colnames
putexcel T28 = matrix(Fx_Levine_cells), colnames

putexcel A48 = matrix(Fx_mpoa_smk), names
putexcel H48 = matrix(Fx_Horvath_smk), colnames
putexcel N48 = matrix(Fx_Hannum_smk), colnames
putexcel T48 = matrix(Fx_Levine_smk), colnames
//****************************************************************************************//
//****************************************************************************************//


//****************************************************************************************//
//****************************************************************************************//
//FIGURE S1C. Comparison of Effect-sizes between DunedinPoAm & Epigenetic Clocks
//****************************************************************************************//
clear
svmat2 Fx_mpoa, names(col) rnames(DV)
gen X = 1 
gen Pheno=_n
save temp, replace
clear
svmat2 Fx_Horvath, names(col) rnames(DV)
gen X = 2 
gen Pheno=_n
append using temp
save temp, replace 
clear
svmat2 Fx_Hannum, names(col) rnames(DV)
gen X = 3
gen Pheno=_n
append using temp
save temp, replace 
clear
svmat2 Fx_Levine, names(col) rnames(DV)
gen X = 4
gen Pheno=_n
append using temp
save temp, replace 

use temp, clear
drop if inlist(Pheno,11,12)
recode Pheno (13=11) (14=12) 
capture label drop Pheno
#delimit ;
	label define Pheno 
		1 "Balance"
		2 "Gait Speed"
		3 "Steps in Place"
		4 "Chair Stands"
		5 "Grip Strength"
		6 "Motor Coordination"
		7 "Physical Limitations"
		8 "Perceptual Reasoning"
		9 "Working Memory"
		10 "Proceessing Speed"
		11 "Self-Rated Health"
		12 "Facial Aging"
		; #delimit cr
foreach x in r lb ub{
	replace `x' = `x'*-1 if inlist(Pheno,7,12)
	replace `x' = `x'*-1	
	}

label values Pheno Pheno
capture label drop X
label define X 1 "DunedinPoAm" 2 "Horvath Clock" 3 "Hannum Clock" 4 "Levine Clock"
label values X X
sort Pheno X
gen Z =  _n
replace Z = Z+3 if Z>4
replace Z = Z+3 if Z>11
replace Z = Z+3 if Z>18
replace Z = Z+3 if Z>25
replace Z = Z+3 if Z>32
replace Z = Z+3 if Z>39
replace Z = Z+3 if Z>46
replace Z = Z+3 if Z>53
replace Z = Z+3 if Z>60
replace Z = Z+3 if Z>67
replace Z = Z+3 if Z>74

#delimit ; 
twoway scatter Z r if X ==1, msize(medium) mcolor(midblue)
	|| scatter Z r if X ==2, msize(medium) mcolor(black)
	|| scatter Z r if X ==3, msize(medium) mcolor(gs6)
	|| scatter Z r if X ==4, msize(medium) mcolor(purple)
	|| rcap lb ub Z if X==1, horiz lwidth(medthick) lcolor(midblue)
	|| rcap lb ub Z if X==2, horiz lwidth(medthick) lcolor(black)
	|| rcap lb ub Z if X==3, horiz lwidth(medthick) lcolor(gs6)
	|| rcap lb ub Z if X==4, horiz lwidth(medthick) lcolor(purple)
	graphregion(color(white)) plotregion(color(white))
	ylabel(
1 "Balance"
8 "Gait Speed"
15 "Steps in Place"
22 "Chair Stands"
29 "Grip Strength"
36 "Motor Coordination"
43 "Physical Limitations"
50 "Perceptual Reasoning"
57 "Working Memory"
64 "Proceessing Speed"
71 "Self-Rated Health"
78 "Facial Aging"
		,nogrid angle(horiz) noticks labsize(medium)) 
	yscale(reverse)
	ytitle("")
	legend(pos(3) cols(1) order(1 2 3 4) region(lcolor(white)) 
		lab(1 "DunedinPoAm") lab(2 "Horvath Clock") lab(3 "Hannum Clock") lab(4 "Levine Clock")
		size(medsmall)
		)
	xline(0, lcolor(gs8))
	xsize(6) ysize(4)
	xtitle("Effect-size (Pearson r)", size(medlarge) margin(small))
; #delimit cr
 
//****************************************************************************************//
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/SuppFigure1C.pdf", replace	
//****************************************************************************************//
//****************************************************************************************//


//****************************************************************************************//
//ANALYSIS OF CHANGE
//****************************************************************************************//
global user db3275
use "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_Dunedin.dta", clear	
	
gen Horvath=horvathage38-ExactAge38 
gen Hannum=hannumage38-ExactAge38
gen Levine=LevineAge38-ExactAge38

global DPHENO "DB DG d_gpdom DL Diq D_A D_V iPoor DF" //"DB DG d_gpdom DL Diq d_Art d_Ven iPoor DF"

	//Regression Analysis
foreach y in mpoa Horvath Hannum Levine{
	capture drop z`y'
	egen z`y' =std(`y')
	matrix DFx_`y' = J(1,5,999)
	matrix DFx_`y'_cells = J(1,5,999)
	matrix DFx_`y'_smk = J(1,5,999)
	foreach x in $DPHENO{
		if `"`x'"'=="iPoor"{
				//Base Model
			quietly poisson `x' z`y' sex, robust 
			matrix A = exp(_b[z`y']) ,  exp(_b[z`y'] - invnormal(0.975)*_se[z`y']), exp(_b[z`y'] + invnormal(0.975)*_se[z`y']), 2*normal(-abs(_b[z`y']/_se[z`y'])), e(N)
			matrix rownames A = `x'
			matrix DFx_`y' = DFx_`y' \ A			
				//Adjusted for Cells
			quietly poisson `x' z`y' sex plasmablast38 cd8pcd28ncd45ran38 cd8naive38 cd4naive38 cd8t38 cd4t38 nk38 bcell38 mono38 gran38, robust 
			matrix A = exp(_b[z`y']) ,  exp(_b[z`y'] - invnormal(0.975)*_se[z`y']), exp(_b[z`y'] + invnormal(0.975)*_se[z`y']), 2*normal(-abs(_b[z`y']/_se[z`y'])), e(N)
			matrix rownames A = `x'
			matrix DFx_`y'_cells = DFx_`y'_cells \ A
				//Adjusted for Smoking
			quietly poisson `x' z`y' sex PackYrLifTm38, robust 
			matrix A = exp(_b[z`y']) ,  exp(_b[z`y'] - invnormal(0.975)*_se[z`y']), exp(_b[z`y'] + invnormal(0.975)*_se[z`y']), 2*normal(-abs(_b[z`y']/_se[z`y'])), e(N)
			matrix rownames A = `x'
			matrix DFx_`y'_smk = DFx_`y'_smk \ A
			}	
		else{
		capture drop Z 
		gen Z =.
		foreach s in 1 2{
			quietly sum `x' if sex==`s'
			replace Z = (`x'-r(mean))/r(sd) if sex==`s'
			}
		//Base Model
		quietly reg Z z`y' sex, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix DFx_`y' = DFx_`y' \ A
		//Adjusted for Cells
		quietly reg Z z`y' sex plasmablast38 cd8pcd28ncd45ran38 cd8naive38 cd4naive38 cd8t38 cd4t38 nk38 bcell38 mono38 gran38, robust 		
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)	
		matrix rownames A = `x'
		matrix DFx_`y'_cells = DFx_`y'_cells \ A
		//Adjusted for Smoking (pack years)
		quietly reg Z z`y' sex PackYrLifTm38, robust 
		matrix A = _b[z`y'] , _b[z`y'] - invttail(e(df_r),0.025)*_se[z`y'], _b[z`y'] + invttail(e(df_r),0.025)*_se[z`y'], 2*ttail(e(df_r),abs(_b[z`y']/_se[z`y'])), e(N)
		matrix rownames A = `x'
		matrix DFx_`y'_smk = DFx_`y'_smk \ A
		drop Z
		}
		}
		drop z`y'
	foreach q in "" _cells _smk{ 
		matrix DFx_`y'`q' = DFx_`y'`q'[2...,1...]
		matrix colnames DFx_`y'`q' = r lb ub p N
		matrix list DFx_`y'`q' 
		}
	}

//Post Tables to Excel
putexcel set "/Users/$user/Box/Belsky/mPoA/Documentation/mPoATables_R1.xlsx", sheet(D3845Comp) modify
	putexcel A45 = "mPoA Regression Analysis of Change Scores"
	putexcel A46 = matrix(DFx_mpoa), names
	putexcel H46 = matrix(DFx_Horvath), colnames
	putexcel N46 = matrix(DFx_Hannum), colnames
	putexcel T46 = matrix(DFx_Levine), colnames
	
	putexcel A57 = "mPoA Regression Analysis of Change Scores - Adjusted for Estimated Cell Counts at Age 38"
	putexcel A58 = matrix(DFx_mpoa_cells), names
	putexcel H58 = matrix(DFx_Horvath_cells), colnames
	putexcel N58 = matrix(DFx_Hannum_cells), colnames
	putexcel T58 = matrix(DFx_Levine_cells), colnames

	putexcel A69 = "mPoA Regression Analysis of Change Scores - Adjusted for Smoking at Age 38 (pack years)"	
	putexcel A70 = matrix(DFx_mpoa_smk), names
	putexcel H70 = matrix(DFx_Horvath_smk), colnames	
	putexcel N70 = matrix(DFx_Hannum_smk), colnames	
	putexcel T70 = matrix(DFx_Levine_smk), colnames	
	


//****************************************************************************************//
//FIGURE S1D. Comparison of Change Score Effect-sizes between DunedinPoAm & Epigenetic Clocks
//****************************************************************************************//
clear
svmat2 DFx_mpoa, names(col) rnames(DV)
gen X = 1 
gen Pheno=_n
save temp, replace
clear
svmat2 DFx_Horvath, names(col) rnames(DV)
gen X = 2 
gen Pheno=_n
append using temp
save temp, replace 
clear
svmat2 DFx_Hannum, names(col) rnames(DV)
gen X = 3
gen Pheno=_n
append using temp
save temp, replace 
clear
svmat2 DFx_Levine, names(col) rnames(DV)
gen X = 4
gen Pheno=_n
append using temp
save temp, replace 

use temp, clear
drop if inlist(Pheno,6,7)
recode Pheno (8=6) (9=7) 
capture label drop Pheno
#delimit ;
	label define Pheno 
		1 "Change in Balance 38-45"
		2 "Change in Grip Strength 38-45"
		3 "Change in Motor Coordination 38-45"
		4 "Change in Physical Limitations 38-45"
		5 "Change in Cognitive Decline 13-45"
		6 "Incident Fair/Poor Health 38-45"
		7 "Change in Facial Aging 38-45"
		; #delimit cr
label values Pheno Pheno
capture label drop X
label define X 1 "DunedinPoAm" 2 "Horvath Clock" 3 "Hannum Clock" 4 "Levine Clock"
label values X X

#delimit ; 
	twoway scatter X r if X == 1, msize(large) mcolor(midblue)
		|| scatter X r if X == 2, msize(large) mcolor(black)
		|| scatter X r if X == 3, msize(large) mcolor(gs6)
		|| scatter X r if X == 4, msize(large) mcolor(purple)
		|| rcap lb ub X if X == 1, horiz lcolor(midblue) lwidth(thick)
		|| rcap lb ub X if X == 2, horiz lcolor(black) lwidth(thick)
		|| rcap lb ub X if X == 3, horiz lcolor(gs6) lwidth(thick)
		|| rcap lb ub X if X == 4, horiz lcolor(purple) lwidth(thick)
		by(Pheno, rescale plotregion(color(white)) graphregion(color(white)) note("") holes(6)
			) 
	ylabel( 1 " " 2 " " 3 " " 4 " ", nolabels nogrid angle(horiz) noticks labsize(medium)) 
	yscale(reverse)
	ytitle("")
	legend(pos(3) cols(1) order(1 2 3 4) region(lcolor(white))  
		lab(1 "DunedinPoAm") lab(2 "Horvath Clock") lab(3 "Hannum Clock") lab(4 "Levine Clock")
		size(medsmall)
		)
	xline(0, lcolor(gs8))
	xsize(6) ysize(4)
	xtitle("Effect-size", size(medlarge) margin(small))
	subtitle(,fcolor(white) lcolor(white) size(medsmall))	
; #delimit cr
		

//****************************************************************************************//
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/SuppFigure1D.pdf", replace	
//****************************************************************************************//
//****************************************************************************************//


//****************************************************************************************//
//FIGURE S1E. Comparison of Change Score Effect-sizes between DunedinPoAm & Epigenetic Clocks
//****************************************************************************************//

use "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_Dunedin.dta", clear	
merge 1:1 snum using "/Users/$user/Box/Belsky/PNAS2015/Documentation/RenatePoA.dta", nogen
corr mpoa paceofaging poaslp*

matrix A = r(C)
matrix A = A[3...,1..2]
matrix list A
clear
svmat2 A, names(col) rnames(biomarker)
split biomarker, parse(_) gen(bm)
sort bm2


replace bm2="totchol" if bm2=="chol"
replace bm2="creat" if bm2=="crcl"
replace bm2="gum" if bm2=="cal"
replace bm2="bun" if bm2=="urea"
replace bm2="leuktl" if bm2=="tel"
replace bm2="cardrespfit" if bm2=="vox"
replace bm2="hscrp" if bm2=="crp"
replace bm2="hba1c" if bm2=="gly"
replace bm2="waisthip" if bm2=="wh"
sort bm2
gen marker=_n
capture label drop marker

#delimit ;
twoway lfit mpoa paceofaging 
	|| scatter mpoa paceofaging, scheme(s1mono)
		mcolor(dknavy)
		msymbol(O)
		msize(vlarge)
		mlabel(marker)
		mlabpos(0) mlabcolor(white) mlabsize(small)
	ylabel(,angle(horiz) nogrid)
	legend(off)
	xtitle("Correlation with Age 26-38 Pace of Aging", size(medlarge) margin(small))
	ytitle("Correlation with Age 38 DunedinPoAm", size(medlarge) margin(small))
; #delimit cr

#delimit ;
label define marker
1 "ApoB100/ApoA1"
2 "BMI"
3 "BUN"
4 "Cardiorespiratory Fitness"
5 "Creatinine Clearance"
6 "FEV1"
7 "FEV1/FVC"
8 "Gum health"
9 "HbA1C"
10 "HDL cholesterol"
11 "hsCRP"
12 "Leukocyte Telomere Length"
13 "Lipoprotein (a)"
14 "Mean Arterial Pressure"
15 "Total Cholesterol"
16 "Triglycerides"
17 "Waist-hip ratio"
18 "White blood cell count"
; #delimit cr
label values marker marker
tab marker

label var mpoa DunedinPoAm
label var paceofaging "Pace of Aging"
label var marker "Slope of Change 26-38"
drop biomarker bm1 bm2
order marker mpoa paceofaging
list

graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/FigureS1E.pdf", replace	

