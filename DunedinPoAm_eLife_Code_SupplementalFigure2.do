global user db3275

//****************************************************************//
//Training Analysis
//****************************************************************//
import delimited using "/Users/$user/Box/DPPP_genomics/Genomics_OUT/ElasticNetProject/FromJoey/poa-hemi-all-20180322/paceofaging-dunedin38-allnox-base-alpha0.5/paceofaging-dunedin38-allnox-base-alpha0.5-bootstats.csv", clear varn(1)
save "/Users/$user/Box/Belsky/mPoA/Documentation/CodeToPrepareData/mPoABootStrapStats.dta", replace

use "/Users/$user/Box/Belsky/mPoA/Documentation/CodeToPrepareData/mPoABootStrapStats.dta", clear
//Number of Training Models
count 
matrix A = r(N)
//Number of Probes used 
sum probes_used
matrix A = A \ r(mean)
//Training Sample Size
sum sample_used
matrix A = A \ r(mean)
matrix rownames A = Reps Probes N
matrix Training1=A

//Probes Selected
sum selected_probes
matrix A = r(mean), r(sd), r(min) , r(max)
//Lamda
sum lambda
matrix A = A \ (r(mean), r(sd), r(min) , r(max))
//Test Corr
sum testing_cor
matrix A = A \ (r(mean), r(sd), r(min) , r(max))
matrix colnames A = M SD Min Max
matrix rownames A = Probes Lambda TestCorr
matrix Training2=A

matrix list Training1
matrix list Training2

putexcel set "/Users/db3275/Box/Belsky/mPoA/Documentation/mPoATables_Dec19.xlsx", sheet(Training) modify
putexcel A1 = matrix(Training1), names
putexcel A5=matrix(Training2), names

cd "/Users/db3275/Box/Belsky/mPoA/Documentation/Figures"
//# of Probes
quietly sum selected_probes
local B = r(mean)
local A = round(r(mean),1)
#delimit ;
hist selected_probes,  freq bin(20) scheme(s1mono) 
	color(purple%20)
	ylabel(,angle(horiz) labsize(medlarge))
	ytitle(Frequency, margin(small) size(medlarge))
	xlabel(,labsize(medlarge) format(%9.0f))
	xtitle(Number of Probes Selected , size(medlarge) margin(small))
	xscale(alt)
	xline(`B', lcolor(red) lwidth(thick))
	text(25 55 `"Mean Number of Probes =`A'"', color(red) placement(e))
	name(num_probes, replace)
	; #delimit cr
graph export "BootstrapProbeNumHist.pdf", replace
	
//Lambda
quietly sum lambda_min
local B= r(mean)
local A= round(r(mean),.01)
#delimit ;
hist lambda_min,  freq bin(20) scheme(s1mono) 
	color(dkorange%30)
	ylabel(,angle(horiz) labsize(medlarge))
	ytitle(Frequency, margin(small) size(medlarge))
	xlabel(,labsize(medlarge))
	xtitle(Lamda , size(medlarge) margin(small))
	xscale(alt)
	xline(`B', lcolor(red) lwidth(thick))
	text(14 .085 `"Mean Lambda = `A'"', placement(e) color(red))
	name(lambda, replace)
	; #delimit cr	
graph export "BootstrapLambdaHist.pdf", replace
	
	
//Corr w PoA
sum testing_cor
local B= r(mean)
#delimit ;
hist testing_cor,  freq bin(20) scheme(s1mono) 
	color(blue%30)
	ylabel(,angle(horiz) labsize(medlarge))
	ytitle(Frequency, margin(small) size(medlarge))
	xlabel(,labsize(medlarge))
	xtitle(mPoA - Pace of Aging Correlation in Test Data , size(medlarge) margin(small))
	xscale(alt)
	xline(`B', lcolor(red) lwidth(thick))
	text(14 .34 "Mean Corr w/" "Pace of Aging=0.33", justification(left) placement(e) color(red))
	name(testcor, replace)
	; #delimit cr	
graph export "BootstrapFitHist.pdf", replace
	
//#Probes - Corr
#delimit ;
scatter testing_cor  selected_probes, scheme(s1mono) mcolor(blue%30)
	ylabel(,angle(horiz) labsize(medlarge))
	ytitle(Test Sample mPoA - Pace of Aging Correlation, margin(small) size(medlarge))
	xlabel(,labsize(medlarge) format(%9.0f))
	xtitle(Number of Probes Selected , size(medlarge) margin(small))
	|| fpfit testing_cor  selected_probes, lcolor(red%30) lwidth(medthick) 
	legend(off)
	name(numprobes_testcor, replace)
	; #delimit cr
graph export "Fit_by_ProbeNum.pdf", replace
	
//Lamda - #Probes 
#delimit ;
scatter selected_probes lambda_min, scheme(s1mono) mcolor(blue%30)
	ylabel(,angle(horiz) labsize(medlarge))
	ytitle(Number of Probes Selected, margin(small) size(medlarge))
	xlabel(,labsize(medlarge) format(%9.2f))
	xtitle(Lamda, size(medlarge) margin(small))
	|| lfit selected_probes lambda_min, lcolor(red%30) lwidth(medthick) lpattern(dash)
	legend(off)
	name(lambda_numprobes, replace)
	; #delimit cr	
graph export "Lamda_by_ProbeNum.pdf", replace

//Correlations among Different Bootsrap mPoA Models
cd "/Users/db3275/Box/DPPP_genomics/Genomics_OUT/ElasticNetProject/FromJoey/bootstrap-results/paceofaging-dunedin38-allnox-base-alpha0.5"
//Assemble dataset of all holdout sample predicted mPoA values for all 100 models
import delimited "r00-scores.tsv", delim(tab) varn(1) clear
rename score mpoa0
forvalues v=1(1)9{
	preserve
	import delimited `"r0`v'-scores.tsv"', delim(tab) varn(1) clear
	rename score mpoa`v'
	save temp, replace
	restore
	merge 1:1 id using temp, nogen
	}
forvalues v=10(1)99{
	preserve
	import delimited `"r`v'-scores.tsv"', delim(tab) varn(1) clear
	rename score mpoa`v'
	save temp, replace
	restore
	merge 1:1 id using temp, nogen
	}	
//Assemble dataset of correlations among the different models
set matsize 11000
matrix Fx = J(1,4,999)
forvalues v1=0(1)99{
	forvalues v2=0(1)99{
		count if mpoa`v1'!=. & mpoa`v2'!=.
		matrix A = `v1', `v2', r(N)
		if r(N)>=3{
			corr mpoa`v1' mpoa`v2'
			matrix A=A,r(rho)
			}
		if r(N)<3{
			matrix A = A,999
			}
		matrix Fx = Fx \ A
		}
		}
clear
matrix colnames Fx = m1 m2 N r
svmat Fx, names(col)
drop if m1==999
drop if m1==m2
count
replace r = . if r==999 // 68 pairs with no match = 34 unique combinations
gen U =.
local a = 0 
local b = 1
forvalues v = 1(1)99{
	replace U = 1 if m1==`a' & m2>=`b'
	local a = `a' + 1
	local b = `b' + 1
	}
keep if U == 1
drop U
save 
use "/Users/db3275/Box/Belsky/mPoA/Documentation/CodeToPrepareData/Dunedin_mpoabootsrap_Corr.dta", clear

sum r, det
#delimit ;
twoway kdensity r , bwidth(.01) scheme(s1mono) kerne(gaus)
	color(blue) lwidth(thick)
	ylabel(,angle(horiz) labsize(medlarge))
	ytitle(Density for Correlation, margin(small) size(medlarge) axis(1))
	xlabel(,labsize(medlarge))
	xtitle(Pairwise Correlation Among Bootstrap mPoA (blue), size(medlarge) margin(small))
	xscale(alt axis(1)) yscale(alt axis(1))
	name(bootstrapmpoa_r, replace)
	note(n=34 pairs with no overlap)
	|| kdensity N, kernel(gaus) bwidth(2) xaxis(2) yaxis(2) 
		color(dknavy) xscale(alt axis(2)) yscale(alt axis(2)) 
		ytitle(Density for Sample Size, margin(small) size(medlarge) axis(2))
		xtitle(Sample Size, axis(2))
		ylabel(,angle(horiz) labsize(medium) format(%9.2f) axis(2))
		xlabel(,labsize(medium) axis(2))
	legend(ring(0) pos(11) cols(1) symxsize(5) lab(1 "Correlation") lab(2 "Sample Size") region(lcolor(white)) )
	; #delimit cr
restore
cd "/Users/db3275/Box/Belsky/mPoA/Documentation/Figures"
graph export "Bootstrap_ePoA_intercorr.pdf", replace

//Combined Figure
graph combine lambda num_probes bootstrapmpoa_r testcor , scheme(s1mono) altshrink name(epoaboot, replace)
graph export "PoA_BootCombined.pdf", replace
