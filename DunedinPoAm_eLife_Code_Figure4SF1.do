clear
input N	ba	fx	lb	ub	dv
1					
2					
3	1	1.29	1.16	1.45	1
4	1	1.19	1.03	1.38	2
5	1	1.15	1.11	1.2	3
6					
7					
8	2	1.07	0.94	1.22	1
9	2	1.02	0.86	1.2	2
10	2	1.07	1.03	1.12	3
11					
12					
13	3	1.05	0.9	1.23	1
14	3	1.17	0.99	1.39	2
15	3	1.04	1	1.09	3
16					
17					
18	4	1.18	1.05	1.34	1
19	4	1.1	0.93	1.31	2
20	4	1.13	1.08	1.18	3
end 
capture label drop ba
label define ba 1 "DunedinPoAm" 2 "Horvath" 3 "Hannum" 4 "PhenoAge"
label values ba ba
capture label drop dv
label define dv 1 "Mortality" 2 "Incident Chronic Disease" 3 "Prevalent Chronic Disease"
label values dv dv

sort dv ba
gen X = _n
replace X = X +2 if X>4
replace X = X + 2 if X>10

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
	xline(1, lcolor(gs8) lpattern(dash))
	xsize(6) ysize(4)
	xtitle("Effect-size", size(large) margin(small))
	xlabel(.9(.1)1.5,labsize(large))
	text(0 1.2 "Mortality (HR)", size(large) color(black))
	text(5.5 1.2 "Incident Chronic Disease (HR)", size(large) color(black))
	text(11.2 1.2 "Prevalent Chronic Disease (IRR)", size(large) color(black))
	; #delimit cr
	
//****************************************************************************************//
graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/SuppFigure4.pdf", replace	
//****************************************************************************************//


	


