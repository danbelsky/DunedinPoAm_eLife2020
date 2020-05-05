
//****************************************************************//
//E-Risk
//****************************************************************//	
global user db3275
use "/Users/$user/Box/Belsky/mPoA/Documentation/CodeToPrepareData/mPoA_ERISKpheno.dta", clear
merge 1:1 atwinid using "/Users/$user/Box/Belsky/mPoA/Documentation/CodeToPrepareData/mPoA_ERISKmeth.dta", nogen

putexcel set "/Users/$user/Box/Belsky/mPoA/Documentation/mPoATables_R1.xlsx", sheet(ERiskComp) modify

gen Horvath = dnamage-age
gen Hannum = bioage4hastatic - age
gen Levine = LevineAge - age

global cells "plasmablast cd8pcd28ncd45ran cd8naive cd4naive cd8t cd4t nk bcell mono gran"
 
foreach x in Horvath Hannum Levine mpoa{
	if `"`x'"' == "mpoa" {
		local Y ""
		}
	else{
		local Y "Age Diff"
		}
	egen z`x' = std(`x')
	label var z`x' `"Z-score of `x'"'
	reg `x' mpc*
	predict br`x', r
	label var br`x' `"`x' `Y' resid for Batch Effects (tech PCs)"'
	egen zbr`x'=std(br`x')
	label var zbr`x' `"Z-score of `x' `Y' resid for Batch Effects (tech PCs)"'
	reg `x' mpc** $cells
	predict bcr`x', r
	label var bcr`x' `"`x' `Y'resid for Batch Effects (tech PCs) & Cell Counts"'
	egen zbcr`x' = std(bcr`x')
	label var zbcr`x' `"Z-score of `x' `Y' resid for Batch Effects (tech PCs) & Cell Counts"'
	}	
	

label var age "Chronological Age"
label var dnamage "Horvath Clock"
label var bioage4hastatic "Hannum Clock"
label var LevineAge "Levine Clock"
label var mpoa "mPoA"

label var brHorvath "Horvath Clock"
label var brHannum "Hannum Clock"
label var brLevine "Levine Clock"
label var brmpoa "mPoA"

label var bcrHorvath "Horvath Clock (adj cells)"
label var bcrHannum "Hannum Clock (adj cells)"
label var bcrLevine "Levine Clock (adj cells)"
label var bcrmpoa "mPoA (adj cells)"

label var Horvath "Horvath Clock Age Difference"
label var Hannum "Hannum Clock Age Difference"
label var Levine "Levine Clock Age Difference"


save "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_ERisk.dta", replace 	

use "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_ERisk.dta", clear
	
//Summary Statistics	
putexcel A1="ERISK methylation N"
sum mpoa
putexcel B1 = matrix(r(N))
putexcel A3 = "Complete Pairs"
matrix Fx = (0,0)
preserve
	egen mdata= count(mpoa), by(familyid)
	foreach x in 1 2{
		unique familyid if Zyg==`x' & mdata==2
		matrix A = r(sum)
		unique familyid if Zyg==`x' & mdata==2 & mpoa!=. 
		matrix A = A, r(sum)
		matrix Fx = Fx \ A
		}
	matrix Fx = Fx[2...,1...]
	matrix rownames Fx = MZ DZ
	matrix colnames Fx = All wMethylation
	matrix list Fx
	putexcel A4=matrix(Fx), names

	gen male = sampsex==1
	tabstat age male mpoa dnamage bioage4hastatic LevineAge brmpoa brHorvath brHannum brLevine if mpoa!=., s(mean sd min max n ) save
	matrix X = r(StatTotal)'
	matrix list X
	putexcel A10=matrix(X), names
restore
	//SES
matrix Fx = 999,999
foreach x in 1 2 3{
	count if ses==`x'
	matrix A = r(N)
	count if ses==`x' & mpoa!=.
	matrix A = A,r(N)
	matrix Fx = Fx\A
	}
matrix Fx = Fx[2...,1...]
matrix colnames Fx = Nall Nmeth
matrix rownames Fx = Low Middle High
matrix list Fx
putexcel A22="SES"
putexcel A23=matrix(Fx), names
	//Victimization
matrix Fx = 999,999
forvalues x=0(1)5{
	count if polyve512 ==`x'
	matrix A = r(N)
	count if polyve512==`x' & mpoa!=.
	matrix A = A,r(N)
	matrix Fx = Fx\A
	}
matrix Fx = Fx[2...,1...]
matrix colnames Fx = Nall Nmeth
matrix rownames Fx = Low Middle High
matrix list Fx
putexcel A30="Victimization"
putexcel A31=matrix(Fx), names


//Summary Statistics
tabstat mpoa brmpoa bcrmpoa	, s(mean sd n) save
matrix A = r(StatTotal)
putexcel A59="mPoA"
putexcel A60=matrix(A), names


//SES & VIC Fx
preserve
egen zses=std(seswq) if mpoa!=.
egen zvic=std(polyve512c) if mpoa!=.
foreach x in mpoa{
foreach z in "" _cells _smk _nsmk { 
	matrix Fx`z'=J(1,5,999)
	matrix colnames Fx`z' = b lb ub p N
	if `"`z'"'=="" {
		global C "sampsex"
		}
	if `"`z'"'=="_cells" {
		global C "sampsex plasmablast cd8pcd28ncd45ran cd8naive cd4naive cd8t cd4t nk bcell mono gran"
		}
	if `"`z'"'=="_smk" {
		global C "sampsex smkcure18 smkdlye18"
		}
	if `"`z'"'=="_nsmk" {
		global C "sampsex if smkcure==0"
		}		
	foreach y in zses zvic{
		quietly reg zbr`x' `y' sampsex $C, cluster(familyid) robust 
		matrix A = _b[`y'] , _b[`y'] - invttail(e(df_r),0.025)*_se[`y'], _b[`y'] + invttail(e(df_r),0.025)*_se[`y'], 2*ttail(e(df_r),abs(_b[`y']/_se[`y'])), e(N)
		matrix rownames A = `y'
		matrix Fx`z' = Fx`z' \ A
		}
	quietly reg zbr`x' zses zvic sampsex $C, cluster(familyid) robust 
		matrix B = ( _b[zses] , _b[zses] - invttail(e(df_r),0.025)*_se[zses], _b[zses] + invttail(e(df_r),0.025)*_se[zses], 2*ttail(e(df_r),abs(_b[zses]/_se[zses])), e(N) \ _b[zvic] , _b[zvic] - invttail(e(df_r),0.025)*_se[zvic], _b[zvic] + invttail(e(df_r),0.025)*_se[zvic], 2*ttail(e(df_r),abs(_b[zvic]/_se[zvic])), e(N) )
		matrix rownames B = Zses_mv Zvic_mv
		matrix Fx`z' = Fx`z' \ B	
	foreach y in "seswq35"  { 
		quietly reg zbr`x' i.`y' $C, cluster(familyid) robust 
		#delimit ;
		matrix A = (_b[2.`y'] , _b[2.`y'] - invttail(e(df_r),0.025)*_se[2.`y'], _b[2.`y'] + invttail(e(df_r),0.025)*_se[2.`y'], 2*ttail(e(df_r),abs(_b[2.`y']/_se[2.`y'])), e(N) \ 
		_b[3.`y'] , _b[3.`y'] - invttail(e(df_r),0.025)*_se[3.`y'], _b[3.`y'] + invttail(e(df_r),0.025)*_se[3.`y'], 2*ttail(e(df_r),abs(_b[3.`y']/_se[3.`y'])), e(N) ); #delimit cr
		matrix rownames A = `y'_2 `y'_3
		matrix Fx`z' = Fx`z' \ A			
		}
	foreach y in "polyve512c" {
		quietly reg zbr`x' i.`y' $C, cluster(familyid) robust 
		#delimit ;
		matrix A = (_b[1.`y'] , _b[1.`y'] - invttail(e(df_r),0.025)*_se[1.`y'], _b[1.`y'] + invttail(e(df_r),0.025)*_se[1.`y'], 2*ttail(e(df_r),abs(_b[1.`y']/_se[1.`y'])), e(N) \
		_b[2.`y'] , _b[2.`y'] - invttail(e(df_r),0.025)*_se[2.`y'], _b[2.`y'] + invttail(e(df_r),0.025)*_se[2.`y'], 2*ttail(e(df_r),abs(_b[2.`y']/_se[2.`y'])), e(N) \ 
		_b[3.`y'] , _b[3.`y'] - invttail(e(df_r),0.025)*_se[3.`y'], _b[3.`y'] + invttail(e(df_r),0.025)*_se[3.`y'], 2*ttail(e(df_r),abs(_b[3.`y']/_se[3.`y'])), e(N) ); #delimit cr
		matrix rownames A = `y'_1 `y'_2 `y'_3
		matrix Fx`z' = Fx`z' \ A	
		}
matrix Fx`z' = Fx`z'[2...,1...]
matrix list Fx`z'
		}
		}

putexcel A75 = "REGRESSIONS"
putexcel A76 = "Base Model"
putexcel A77 = matrix(Fx), names
putexcel A88 = "Adjusted for Cells"
putexcel A89 = matrix(Fx_cells), names
putexcel A100 = "Adjusted for Smoking"
putexcel A101 = matrix(Fx_smk), names
putexcel A112 = "Non-Smokers"
putexcel A113 = matrix(Fx_nsmk), names


//******************************************************************************************//
//FIGURE SHOWING MPOA BY SES 
	//Compute effect-sizes for strata for graphing
reg zbrmpoa	i.ses sampsex, robust  cluster(familyid)
	margins, over(ses)
	#delimit ;
	marginsplot, graphregion(color(white) ) plotregion(color(white) )
		recast(bar)
		ytitle("DunednPoAm z-score, Age 18", margin(small) size(medlarge))
		xtitle(Family Social Class, margin(small) size(medlarge))
		ylabel(-0.4(.2).6, format(%9.1f) angle(horiz) labsize(medlarge) nogrid)
		xlabel(1 "Low" 2 "Middle" 3 "High", labsize(medlarge))
		plot1opts(color(blue) lcolor(white))
		ci1opts(color(gs12) lwidth(thick))
		title("A." " " " ", color(black) size(vlarge) position(11) span ring(10))
		name(sesmargins, replace)
	; #delimit cr
//******************************************************************************************//
//FIGURE SHOWING MPOA BY Victimization 
	//Compute effect-sizes for strata for graphing
reg zbrmpoa	i.polyve512c sampsex, robust  cluster(familyid)
	margins, over(polyve512c)
	#delimit ;
	marginsplot, graphregion(color(white)) plotregion(color(white)) 
		recast(bar)
		ytitle("DunedinPoAm z-score, Age 18", margin(small) size(medlarge))
		xtitle("Extent of Polyvictimization (types)" , margin(small) size(medlarge))
		ylabel(-0.4(.2).6, format(%9.1f) angle(horiz) labsize(medlarge) nogrid)
		xlabel(0 "None" 1 "1" 2 "2" 3 "3+", labsize(medlarge))
		plot1opts(color(blue) lcolor(white))
		ci1opt(color(gs10) lwidth(thick))
		title("B." " " " ", color(black) size(vlarge) position(11) span ring(10))
		name(vicmargins, replace)
	; #delimit cr
//******************************************************************************************//
//Produce combined graph
graph combine sesmargins vicmargins, scheme(s1mono) ycommon 
	//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/Figures/Figure5.pdf", replace	
//******************************************************************************************//		
		
		


		
		
		
//****************************************************************************************//
//E-RISK COMPARISON WITH EPIGENETIC CLOCKS
//****************************************************************************************//

global user db3275
use "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_ERisk.dta", clear
putexcel set "/Users/$user/Box/Belsky/mPoA/Documentation/mPoATables_R1.xlsx", sheet(ERiskComp) modify


//SES & VIC Fx
preserve
egen zses=std(seswq) if mpoa!=.
egen zvic=std(polyve512c) if mpoa!=.
//foreach x in mpoa Horvath Hannum Levine {
foreach x in mpoa Horvath Hannum Levine {
foreach z in "" _cells _smk _nsmk { 
	matrix Fx_`x'`z'=J(1,5,999)
	matrix colnames Fx_`x'`z' = b lb ub p N
	if `"`z'"'=="" {
		global C "sampsex"
		}
	if `"`z'"'=="_cells" {
		global C "sampsex plasmablast cd8pcd28ncd45ran cd8naive cd4naive cd8t cd4t nk bcell mono gran"
		}
	if `"`z'"'=="_smk" {
		global C "sampsex smkcure18 smkdlye18"
		}
	if `"`z'"'=="_nsmk" {
		global C "sampsex if smkcure==0"
		}		
	foreach y in zses zvic{
		quietly reg zbr`x' `y' sampsex $C, cluster(familyid) robust 
		matrix A = _b[`y'] , _b[`y'] - invttail(e(df_r),0.025)*_se[`y'], _b[`y'] + invttail(e(df_r),0.025)*_se[`y'], 2*ttail(e(df_r),abs(_b[`y']/_se[`y'])), e(N)
		matrix rownames A = `y'
		matrix Fx_`x'`z' = Fx_`x'`z' \ A
		}
	quietly reg zbr`x' zses zvic sampsex $C, cluster(familyid) robust 
		matrix B = ( _b[zses] , _b[zses] - invttail(e(df_r),0.025)*_se[zses], _b[zses] + invttail(e(df_r),0.025)*_se[zses], 2*ttail(e(df_r),abs(_b[zses]/_se[zses])), e(N) \ _b[zvic] , _b[zvic] - invttail(e(df_r),0.025)*_se[zvic], _b[zvic] + invttail(e(df_r),0.025)*_se[zvic], 2*ttail(e(df_r),abs(_b[zvic]/_se[zvic])), e(N) )
		matrix rownames B = Zses_mv Zvic_mv
		matrix Fx_`x'`z' = Fx_`x'`z' \ B	
	foreach y in "seswq35"  { 
		quietly reg zbr`x' i.`y' $C, cluster(familyid) robust 
		#delimit ;
		matrix A = (_b[2.`y'] , _b[2.`y'] - invttail(e(df_r),0.025)*_se[2.`y'], _b[2.`y'] + invttail(e(df_r),0.025)*_se[2.`y'], 2*ttail(e(df_r),abs(_b[2.`y']/_se[2.`y'])), e(N) \ 
		_b[3.`y'] , _b[3.`y'] - invttail(e(df_r),0.025)*_se[3.`y'], _b[3.`y'] + invttail(e(df_r),0.025)*_se[3.`y'], 2*ttail(e(df_r),abs(_b[3.`y']/_se[3.`y'])), e(N) ); #delimit cr
		matrix rownames A = `y'_2 `y'_3
		matrix Fx_`x'`z' = Fx_`x'`z' \ A			
		}
	foreach y in "polyve512c" {
		quietly reg zbr`x' i.`y' $C, cluster(familyid) robust 
		#delimit ;
		matrix A = (_b[1.`y'] , _b[1.`y'] - invttail(e(df_r),0.025)*_se[1.`y'], _b[1.`y'] + invttail(e(df_r),0.025)*_se[1.`y'], 2*ttail(e(df_r),abs(_b[1.`y']/_se[1.`y'])), e(N) \
		_b[2.`y'] , _b[2.`y'] - invttail(e(df_r),0.025)*_se[2.`y'], _b[2.`y'] + invttail(e(df_r),0.025)*_se[2.`y'], 2*ttail(e(df_r),abs(_b[2.`y']/_se[2.`y'])), e(N) \ 
		_b[3.`y'] , _b[3.`y'] - invttail(e(df_r),0.025)*_se[3.`y'], _b[3.`y'] + invttail(e(df_r),0.025)*_se[3.`y'], 2*ttail(e(df_r),abs(_b[3.`y']/_se[3.`y'])), e(N) ); #delimit cr
		matrix rownames A = `y'_1 `y'_2 `y'_3
		matrix Fx_`x'`z' = Fx_`x'`z' \ A	
		}
matrix Fx_`x'`z' = Fx_`x'`z'[2..3,1...]
matrix list Fx_`x'`z'
		}
		}
		
# delimit ;
matrix A = Fx_mpoa, Fx_Horvath, Fx_Hannum, Fx_Levine \
	Fx_mpoa_cells, Fx_Horvath_cells, Fx_Hannum_cells, Fx_Levine_cells \
	Fx_mpoa_smk, Fx_Horvath_smk, Fx_Hannum_smk, Fx_Levine_smk \
	Fx_mpoa_nsmk, Fx_Horvath_nsmk, Fx_Hannum_nsmk, Fx_Levine_nsmk 
	; #delimit cr

putexcel B2 = matrix(A)	, names


matrix B = (Fx_mpoa, (1,1\1,2)) \ (Fx_Horvath, (2,1\2,2))  \ (Fx_Hannum, (3,1\3,2)) \ (Fx_Levine, (4,1\4,2))
clear
svmat2 B, names(col) rnames(DV)
rename c6 BA 
rename c7 Y
capture label drop BA
label define BA 1 "DunedinPoAm" 2 "Horvath" 3 "Hannum" 4 "Levine"
label values BA BA
capture label d rop Y
label define Y 1 "Low SES" 2 "Victimization"
label values Y Y
sort Y BA 
gen X = _n
replace X = X + 1 if X>4
foreach x in b lb ub{
	replace `x'=`x'*-1 if Y==1
	}
#delimit ; 
	twoway scatter X b if BA == 1, msize(large) mcolor(midblue)
		|| scatter X b if BA == 2, msize(large) mcolor(black)
		|| scatter X b if BA == 3, msize(large) mcolor(gs6)
		|| scatter X b if BA == 4, msize(large) mcolor(purple)
		|| rcap lb ub X if BA == 1, horiz lcolor(midblue) lwidth(thick)
		|| rcap lb ub X if BA == 2, horiz lcolor(black) lwidth(thick)
		|| rcap lb ub X if BA == 3, horiz lcolor(gs6) lwidth(thick)
		|| rcap lb ub X if BA == 4, horiz lcolor(purple) lwidth(thick)
	graphregion(color(white)) plotregion(color(white))
	text(0 .1 "Low SES", size(large))
	text(5 .1 "Victimization", size(large))
	ylabel(, nolabels nogrid angle(horiz) noticks labsize(medium)) 
	yscale(reverse)
	ytitle("")
	legend(ring(0) pos(5) cols(1) order(1 2 3 4) region(lcolor(white))  
		lab(1 "DunedinPoAm") lab(2 "Horvath Clock") lab(3 "Hannum Clock") lab(4 "Levine Clock")
		size(medsmall)
		)
	xline(0, lcolor(gs8) lpattern(dash))
	xsize(6) ysize(4)
	xtitle("Effect-size (Pearson r)", size(large) margin(small))
	xlabel(,labsize(large))
; #delimit cr
		

//****************************************************************************************//
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/SuppFigure5.pdf", replace	
//****************************************************************************************//
//****************************************************************************************//

		
