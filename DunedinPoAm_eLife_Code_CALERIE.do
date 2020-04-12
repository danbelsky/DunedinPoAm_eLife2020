//****************************************************************//
//CALERIE
//****************************************************************//		
global user db3275
use "/Users/$user/Box/Belsky/mPoA/Documentation/CodeToPrepareData/mPoA_CALERIEpheno.dta", clear
merge m:1 deidnum using "/Users/$user/Box/Belsky/mPoA/Documentation/CodeToPrepareData/mPoA_CALERIEmeth.dta", nogen

putexcel set "/Users/$user/Box/Belsky/mPoA/Documentation/mPoATables_R1.xlsx", sheet(CALERIE) modify

tab plate, gen(Plate)
	//Z-scores of mPOA Variables
reg mpoa Plate*
predict rmpoa, r
quietly sum rmpoa if fu==0
gen zrmpoa = (rmpoa-r(mean))/r(sd)

gen Horvath = DNAmAge
gen Hannum = BioAge4HAStatic
gen Levine = LevineAge

foreach x in Horvath Hannum Levine{
	quietly reg `x' age Plate*
	predict X, r
	egen Z = std(X)
	replace `x'=Z
	drop X Z
	}

gen b1 = BA if fu==0
egen baseba=max(b1),by(deidnum)
gen BAgain = BA-baseba
drop b1 baseba
gen baa = BA - age

gen CA=agevis
gen b1 = CA if fu==0
egen baseCA=max(b1),by(deidnum)
gen CAgain = CA-baseCA
drop b1 baseCA

capture label drop CR
label define CR 0 "AL" 1 "CR"
label values CR CR

save "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_CALERIE.dta", replace

//Summary Statistics
capture drop male
gen male = sex==2
tabstat age male BA BAgainslope HD HDgainslope mpoa if fu==0, by(CR) s(mean sd n) save
matrix X = r(Stat1)'  , r(Stat2)' 
matrix list X
putexcel A1="CALERIE"
putexcel A2=matrix(X) , names

tabstat age male BA BAgainslope HD HDgainslope mpoa if fu==0 & mpoa!=., by(CR) s(mean sd n) save
matrix X = r(Stat1)'  , r(Stat2)' 
putexcel A12="CALERIE w/ Methylation"
putexcel A13=matrix(X) , names

//Test of diff in mpoa between CR and AL
reg mpoa CR if fu==0, robust 

//Baseline associations
preserve
	keep if fu==0
	local C
	foreach x in baa HD baseage {
		egen z`x' = std(`x') 
		}
	matrix Fx = J(1,3,999)
	foreach x in baa HD{
		reg mpoa `x' sex baseage Plate* , robust
		matrix A = _b[`x'], _se[`x'], e(N)
		matrix rownames A = `x'
		matrix Fx = Fx \ A
		}
	reg mpoa male baseage Plate* , robust 
	matrix A = (_b[baseage], _se[baseage], e(N)) \ (_b[male],_se[male], e(N))
	matrix rownames A = Age Male
	matrix Fx = Fx \ A
	matrix Fx = Fx[2...,1...]
	matrix colnames Fx = b SE N
	matrix list Fx
restore
putexcel A22 = "Baseline Assoc w mPoA"
putexcel A23 = matrix(Fx), names

//mPoA Association with future rate of aging
preserve
	//Set ref to 38yo female
	capture drop OM*
	gen OMbaseage = (baseage-38)/10
	gen OMsex = sex-1
reg BAgainslope c.zrmpoa##CR OMbaseage OMsex if fu==0, robust 
margins, dydx(zrmpoa) over(CR) 




//Confirmation of substantive result from re-fitting growth model
mixed BAgain CR##c.CAgain##c.zrmpoa  OMsex OMbaseage, nocons || deidnum: CAgain , robust cov(unstr)	
restore

global user db3275
use "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_CALERIE.dta", clear

//FIGURE
	//Original Model from JGBS Paper
preserve
	//Set ref to 38yo female
	capture drop OM*
	gen OMbaseage = (baseage-38)/10
	gen OMsex = sex-1
	mixed BAgain CR##c.CAgain##c.zrmpoa OMsex OMbaseage,  || deidnum: CAgain , robust cov(unstr)	
	margins, over(CR) at(CAgain=(0 1 2) zrmpoa=(-1 1)) 			
	#delimit ;
	marginsplot, 
		by(CR,)
		subtitle(,bcolor(white) size(large) color(black))
		ytitle(Change in KDM Biological Age, margin(small) size(medlarge))
		ylabel(,angle(horiz) labsize(medlarge) format(%9.1f) nogrid)
		xscale(range(-.25 2.25))
		xtitle(Follow-up, margin(small) size(medlarge))
		xlabel(0 "Baseline" 1 "12mo" 2 "24mo",labsize(medlarge))
		ciopt(recast(rarea))
		plot1opt(color(dknavy) msymbol(O) msize(vlarge) lwidth(thick))
		plot2opt(color(midblue) msymbol(T) msize(vlarge) lpattern(dash) lwidth(thick) )
		ci1opt(color(dknavy%20) lcolor(white%5) lwidth(vvthin))
		ci2opt(color(midblue%20) lcolor(white%5) lwidth(vvthin))
		byopts(
			graphregion(color(white)) plotregion(color(white)) 
			legend(cols(1) pos(3) 
				lab(3 "Slow") lab(4 "Fast ") 
				order(4 3) symxsize(5) region(lcolor(white)) size(medlarge))
				title(""))
		legend(cols(1) order(4 3) symxsize(5) lab(3 "Slow") lab(4 "Fast") region(lcolor(white)) size(medlarge) title("Baseline" "DunedinPoAm", color(black) size(medium)))
		name(growthd3, replace)
		; #delimit cr
//******************************************************************************************//
//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/Figure4.pdf", replace	
//******************************************************************************************//




//******************************************************************************************//
//Comparison with Clocks
//******************************************************************************************//

use "/Users/$user/Box/Belsky/mPoA/eLife/R1/mPoAR1_CALERIE.dta", clear

preserve
	//Set ref to 38yo female
	capture drop OM*
	gen OMbaseage = (baseage-38)/10
	gen OMsex = sex-1
	foreach x in zrmpoa Horvath Hannum Levine{	
		foreach Z in 1 2 {
			if `Z'==1 {
				local C ""
				}
			if `Z'==2{
				local C "PlasmaBlast-Gran"
				}
			quietly reg BAgainslope c.`x'##CR OMbaseage OMsex `C' if fu==0, robust 
			quietly margins, dydx(`x') over(CR) 
			matrix A = r(table)
			matrix A = (A[1,1...] \ A[5..6,1...] \ A[4,1...]) 
			matrix Fx_`x'`Z' = A'
			matrix list Fx_`x'`Z'
			}
			}
restore

putexcel A37=matrix(Fx_zrmpoa1), names
putexcel A41=matrix(Fx_Horvath1), names
putexcel A45=matrix(Fx_Hannum1), names
putexcel A49=matrix(Fx_Levine1), names


putexcel A53=matrix(Fx_zrmpoa2), names
putexcel A57=matrix(Fx_Horvath2), names
putexcel A61=matrix(Fx_Hannum2), names
putexcel A65=matrix(Fx_Levine2), names
//******************************************************************************************//

//******************************************************************************************//
//Figure for Clock Comparison
//******************************************************************************************//
matrix A = Fx_zrmpoa, (1,1\1,2) \ Fx_Horvath, (2,1\2,2) \ Fx_Hannum, (3,1\3,2) \ Fx_Levine, (4,1\4,2) 
matrix list A

clear
svmat2 A, names(col) rnames(ba)
rename c5 BA 
rename c6 Y
capture label drop BA
label define BA 1 "DunedinPoAm" 2 "Horvath" 3 "Hannum" 4 "Levine"
label values BA BA
capture label drop Y
label define Y 1 "AL" 2 "CR"
label values Y Y
sort BA Y 
gen X = _n
replace X = X + 1 if X>2
replace X = X + 1 if X>5
replace X = X + 1 if X>8
drop ba
rename BA ba
rename b fx
rename ll lb
rename ul ub
#delimit ;
graph twoway  
		   rcap lb ub X if ba == 1, horiz lcolor(midblue) lwidth(thick)
		|| rcap lb ub X if ba == 2, horiz lcolor(black) lwidth(thick)
		|| rcap lb ub X if ba == 3, horiz lcolor(gs6) lwidth(thick)
		|| rcap lb ub X if ba == 4, horiz lcolor(purple) lwidth(thick)

		|| scatter X fx if ba == 1 & Y==1, msize(large) mcolor(midblue)
		|| scatter X fx if ba == 2 & Y==1, msize(large) mcolor(black)
		|| scatter X fx if ba == 3 & Y==1, msize(large) mcolor(gs6)
		|| scatter X fx if ba == 4 & Y==1, msize(large) mcolor(purple)
		
		|| scatter X fx if ba == 1 & Y==2, msize(large) mlcolor(midblue) mlwidth(medthick) mfcolor(white)
 		|| scatter X fx if ba == 2 & Y==2, msize(large) mlcolor(black) mlwidth(medthick) mfcolor(white)
		|| scatter X fx if ba == 3 & Y==2, msize(large) mlcolor(gs6) mlwidth(medthick) mfcolor(white)
		|| scatter X fx if ba == 4 & Y==2, msize(large) mlcolor(purple)	mlwidth(medthick) mfcolor(white)	
		
	graphregion(color(white)) plotregion(color(white))
	ylabel(1 "AL" 2 "CR" 4 "AL" 5 "CR" 7 "AL" 8 "CR" 10 "AL" 11 "CR"  ,  nogrid angle(horiz) noticks labsize(medium)) 
	yscale(reverse)
	ytitle("")
	legend(ring(0) pos(5) cols(1) order(5 6 7 8) region(lcolor(white))  
		lab(5 "DunedinPoAm") lab(6 "Horvath Clock") lab(7 "Hannum Clock") lab(8 "Levine Clock")
		size(medsmall)
		)
	xline(0, lcolor(gs8) lpattern(dash))
	xsize(6) ysize(4)
	xtitle("Effect-size", size(large) margin(small))
	xlabel(-.5(.1).5,labsize(large))
	; #delimit cr

//******************************************************************************************//
//Export PDF
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/SuppFigure6.pdf", replace	
//******************************************************************************************//





