
//****************************************************************//
//Understanding Society
//****************************************************************//
global user db3275
use "/Users/$user/Box/Belsky/mPoA/Documentation/CodeToPrepareData/mPoA_UNDSOCpheno.dta", clear
merge 1:1 id using "/Users/$user/Box/Belsky/mPoA/Documentation/CodeToPrepareData/mPoA_UNDSOCmeth.dta",nogen

putexcel set "/Users/$user/Box/Belsky/mPoA/Documentation/mPoATables_R1.xlsx", sheet(UndSoc) modify
putexcel A1 = "UNDERSTANDING SOCIETY"

	//Batch Controls
tab rackbar, gen(RB)
capture drop plate
gen plate = .
forvalues x= 1(1)14{
	replace plate = `x' if RB`x'==1
	}
//Batch Residualized mPoA
foreach x in mpoa DNAmAge BioAge4HAStatic LevineAge{
	reg `x' RB*
	predict br`x', r
	replace br`x'=br`x'+_b[_cons]
	label var br`x' `"`x' resid for Plate"'
	egen zbr`x'=std(br`x')
	label var zbr`x' `"Z-score of `x' resid for Plate"'
	}
	
	//Age Residualized Clocks
gen Horvath=DNAmAge
gen Hannum=BioAge4HAStatic
gen Levine=LevineAge 
foreach x in Horvath Hannum Levine{
	local Y "Age Resid"
	reg `x' RB* age
	predict rbr`x', r
	label var rbr`x' `"`x' `Y' resid for Plate"'
	egen zrbr`x'=std(rbr`x')
	label var zrbr`x' `"Z-score of `x' `Y' resid for Plate"'
	}	
egen zmpoa=std(mpoa)	
gen zrbrmpoa=zbrmpoa
drop Horvath Hannum Levine



//Summary Statistics
	//Raw Variables
gen male=sex==2
tabstat age male mpoa DNAmAge BioAge4HAStatic LevineAge BAnochol srh nonsmoker, s(mean sd n) save
return list
matrix A = r(StatTotal)
putexcel A3="Raw Variables"
putexcel A4=matrix(A), names
	//Analysis Values
tabstat brmpoa rbrHorvath rbrHannum rbrLevine baanochol, s(mean sd n) save
return list
matrix A = r(StatTotal)
putexcel A9="Analysis Variables"
putexcel A10=matrix(A), names


//1. mPoA mean & assoc w/ age and sex
sum mpoa 
matrix A = r(mean) , r(sd), 999, 999, r(N)
preserve
	keep if mpoa!=.
	egen zage=std(age)
	reg zbrmpoa zage sex, robust 
	matrix A = A \ _b[zage], _b[zage] - invttail(e(df_r),0.025)*_se[zage], _b[zage] + invttail(e(df_r),0.025)*_se[zage], 2*ttail(e(df_r),abs(_b[zage]/_se[zage])), e(N)
	matrix A = A \ _b[sex], _b[sex] - invttail(e(df_r),0.025)*_se[sex], _b[sex] + invttail(e(df_r),0.025)*_se[sex], 2*ttail(e(df_r),abs(_b[sex]/_se[sex])), e(N)
	matrix colnames A = "r/d" lb ub p N 
	matrix rownames A = mPoA Age Sex
	matrix list A
	putexcel A15= "mPoA, Age, & Sex"
	putexcel A16=matrix(A), names
restore

//2. Correlations ofpub'd clocks with chrono age
corr age brDNAmAge brBioAge4HAStatic brLevineAge
matrix A = r(C)
putexcel A21="Clock Correlations with Chronological Age"
putexcel A22=matrix(A), names

//3. Correlations of age-residuals from pub'd clocks with mPoA
corr zrbrmpoa rbrHorvath rbrHannum rbrLevine 
matrix A = r(C)
putexcel A28="Clock Age-Accel Residual Correlations with mPoA"
putexcel A29=matrix(A), names

//4. Associations of mPoA & Epi Clocks w/ KDM and SRH 
capture drop _BAA_ 
gen _BAA_ = baanochol
//Assoc w/ BA & SRH
foreach x in _BAA_ srh {
	capture drop z`x'
	egen z`x' = std(`x')
	}
foreach z in "" _cells _smk _nsmk{
	matrix Fx`z' = J(1,10,000)
	matrix colnames Fx`z' = b_BA lb ub P N b_SRH lb ub P N
	if `"`z'"'=="" {
		global C "age sex"
		}
	if `"`z'"'=="_cells" {
		global C "age sex PlasmaBlast-Gran"
		}
	if `"`z'"'=="_smk" {
		global C "age sex nonsmoker SM*"
		}
	if `"`z'"'=="_nsmk" {
		global C "age sex if nonsmoker ==1"
		}		
	foreach x in mpoa Horvath Hannum Levine { 
		matrix A`z' = 0
		foreach y in z_BAA_ zsrh{
			quietly reg `y' zrbr`x' $C, robust
			matrix A`z' = A`z', _b[zrbr`x'], _b[zrbr`x'] - (invnormal(.975)*_se[zrbr`x']), _b[zrbr`x'] + (invnormal(.975)*_se[zrbr`x']), 2*ttail(e(df_r),abs(_b[zrbr`x']/ _se[zrbr`x'])), e(N)
			}
		matrix A`z' = A`z'[1...,2...]
		matrix rownames A`z' = `x'
		matrix Fx`z' = Fx`z' \ A`z'
		}
	matrix Fx`z' = Fx`z'[2...,1...]
	matrix list Fx`z'
	}
putexcel A36="REGRESSION MODELS of KDM BA and SRH"
putexcel A37="Effect-sizes for Association with KDM SRH: Standardized Regression Coefs from Models Adjusted for Age & Sex"
putexcel A38=matrix(Fx), names
putexcel A44="Effect-sizes for Association with BA & SRH: Standardized Regression Coefs from Models Adjusted for Age, Sex, and Cells"
putexcel A45=matrix(Fx_cells), names
putexcel A51="Effect-sizes for Association with BA & HD: Standardized Regression Coefs from Models Adjusted for Age, Sex, and smoking"
putexcel A52=matrix(Fx_smk), names
putexcel A58="Effect-sizes for Association with BA & SRH: Standardized Regression Coefs from Models Adjusted for Age, Sex -- NonSmokers Only"
putexcel A59=matrix(Fx_nsmk), names


//******************************************************************************************//
//Supplemental Figure 3. Comparison of KDM BA and SRH associations with methylation variables
//******************************************************************************************//
clear
matrix A = ( Fx[1...,1..5], (1\2\3\4),(1\1\1\1) ) \ (Fx[1...,6..10], (1\2\3\4), (2\2\2\2) )
svmat2 A, names(col) rnames(ba)
rename c6 BA 
rename c7 Y
capture label drop BA
label define BA 1 "DunedinPoAm" 2 "Horvath" 3 "Hannum" 4 "Levine"
label values BA BA
capture label d rop Y
label define Y 1 "KDM BA" 2 "Self-rated Health"
label values Y Y
sort Y BA 
gen X = _n
replace X = X + 1 if X>4
rename b_BA b
foreach x in b lb ub{
	replace `x'=`x'*-1 if Y==2
	}
drop ba
rename BA ba
rename b fx
#delimit ;
graph twoway  scatter X fx if ba == 1, msize(large) mcolor(midblue)
		|| scatter X fx if ba == 2, msize(large) mcolor(black)
		|| scatter X fx if ba == 3, msize(large) mcolor(gs6)
		|| scatter X fx if ba == 4, msize(large) mcolor(purple)
		|| rcap lb ub X if ba == 1, horiz lcolor(midblue) lwidth(thick)
		|| rcap lb ub X if ba == 2, horiz lcolor(black) lwidth(thick)
		|| rcap lb ub X if ba == 3, horiz lcolor(gs6) lwidth(thick)
		|| rcap lb ub X if ba == 4, horiz lcolor(purple) lwidth(thick)
	graphregion(color(white)) plotregion(color(white))
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
	text(0 .1 "KDM Biological Age", size(large) color(black))
	text(5.5 .1 "Self-rated Health (reversed)", size(large) color(black))
	; #delimit cr

//******************************************************************************************//
//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/SuppFig3.pdf", replace	
//******************************************************************************************//



//FIGURES

label var age "Chronological Age"
label var DNAmAge "Horvath Clock"
label var BioAge4HAStatic "Hannum Clock"
label var LevineAge "Levine Clock"
label var mpoa "mPoA"

preserve
label var brmpoa "mPoA"
label var rbrHorvath "Horvath Clock Accel" 
label var rbrHannum "Hannum Clock Accel"
label var rbrLevine "Levine Clock Accel"

//******************************************************************************************//
//Figure 3A. mPoA by Age and Sex
//******************************************************************************************//
	//Linear Slopes w/ CIs
#delimit ;
twoway scatter mpoa age if sex==1,  msymbol(O) mfcolor(gold%50) mlcolor(gold%1) yline(1, lcolor(dknavy))
	|| scatter mpoa age if sex==2,  msymbol(+) mcolor(blue%30) mlwidth(medthick)
	|| lfitci mpoa age if sex==2, lpattern(dash) lcolor(dknavy) lwidth(thick) ciplot(rarea) fcolor(dknavy%20) acolor(none)
	|| lfitci mpoa age if sex==1, lpattern(dash) lcolor(orange) lwidth(thick) ciplot(rarea) fcolor(orange%20) acolor(none)
	$X $Y 
	ylabel(0.8(.2)1.4, format(%9.1f) nogrid angle(horiz)) 
	ytitle(DunedinPoAm, size(medlarge) margin(small))
	graphregion(color(white)) plotregion(color(white))
	legend(pos(3) cols(1) 
		lab(1 "Women") lab(2 "Men") 
		lab(4 "Fitted Slope") lab(6 "Fitted Slope")
		order(2 4 1 6 ) symxsize(5) region(lcolor(white)) )
		name(mpoabyAge_Lin, replace)
; #delimit cr	
//******************************************************************************************//
//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/Figure3A.pdf", replace	
//******************************************************************************************//


//******************************************************************************************//
//Figure 3b. Matrix Plot of Methylation Clocks
//******************************************************************************************//
//Output data for matrix plot figure creation in R
export delimited brmpoa rbrHorvath rbrHannum rbrLevine baanochol age sex using "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/UndSocMatrixPlotData_200306.txt", delim(tab) replace
//******************************************************************************************//


//******************************************************************************************//
//Figure 3C. Associations with KDM BA
//******************************************************************************************//
foreach y in brmpoa rbrHorvath rbrHannum rbrLevine {
	#delimit ;
	if `"`y'"'=="brmpoa"{; local color "midblue" ; local ylabopt ""; local ytitle "DunedinPoAm" ; } ;
	if `"`y'"'=="rbrHorvath"{; local color "black" ; local ylabopt "labcolor(white) noticks"; local ytitle "Horvath Clock" ;} ;
	if `"`y'"'=="rbrHannum"{; local color "gs6" ; local ylabopt ""; local ytitle "Hannum Clock" ;} ;
	if `"`y'"'=="rbrLevine"{; local color "purple" ; local ylabopt "labcolor(white) noticks"; local ytitle "Levine Clock" ;} ; #delimit cr
	#delimit ;
	binscatter baanochol z`y' if abs(z`y')<5, controls(sex) nq(100) lcolor(red) mcolor(`color') 
		ylabel(,angle(horiz) labsize(medium) format(%9.1f) nogrid `ylabopt')
		xlabel(-3(1)3, labsize(medium)) 
		ytitle("")
		xtitle(`ytitle' (z-score units))
		xscale(range(-3 3.5))
		graphregion(color(white))
		plotregion(color(white))
		name(`y', replace)
		; #delimit cr
	}

graph combine brmpoa rbrHorvath rbrHannum rbrLevine, ycommon  scheme(s1mono) name(srh,replace) rows(2) imargin(vsmall) //xsize(4) ysize(4) 
//******************************************************************************************//
//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/Figure3C.pdf", replace	
//******************************************************************************************//

//******************************************************************************************//
//Figure 3D. Associations with Self-rated Health
//******************************************************************************************//
foreach y in brmpoa rbrHorvath rbrHannum rbrLevine {
	#delimit ;
	if `"`y'"'=="brmpoa"{; local color "midblue" ; local ylabopt ""; local ytitle "DunedinPoAm" ; } ;
	if `"`y'"'=="rbrHorvath"{; local color "black" ; local ylabopt "labcolor(white) noticks"; local ytitle "Horvath Clock" ;} ;
	if `"`y'"'=="rbrHannum"{; local color "gs6" ; local ylabopt ""; local ytitle "Hannum Clock" ;} ;
	if `"`y'"'=="rbrLevine"{; local color "purple" ; local ylabopt "labcolor(white) noticks"; local ytitle "Levine Clock" ;} ; #delimit cr
	reg z`y' i.srh sex age, robust
	margins, at(sex=1.5 age=60) over(srh)
	matrix `y'=r(table)
	#delimit ; 
	marginsplot, graphregion(color(white)) plotregion(color(white))
		plot1opts(connect(none) msize(large) mcolor(`color'))
		ci1opts(lwidth(medthick) lcolor(`color'))
		title(`"`ytitle'"', size(medlarge) color(`color') ring(0) pos(1))
		ytitle("", size(medlarge))
		xtitle("")
		ylabel(-.5(.5)1,angle(horiz) format(%9.1f) nogrid labsize(medium) `ylabopt')
		xlabel(,valuelabels angle(50) labsize(medium))
		yline(0, lcolor(gs8))
		name(`y',replace)
		nodraw
		; #delimit cr
}
graph combine brmpoa rbrHorvath rbrHannum rbrLevine, scheme(s1mono) name(srh,replace) rows(2) imargin(vsmall) xsize(4) ysize(4)
//******************************************************************************************//
//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/Figure3D.pdf", replace	
//******************************************************************************************//





