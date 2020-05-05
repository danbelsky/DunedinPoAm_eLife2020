

import delimited using "/Users/db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/NAS/eLife/survival.csv", delim(comma) varn(1) clear

stset time_fu_yr, failure(dead==1)

#delimit ;
sts graph, survival by(poa_all_c) censored(single) censopts(lcolor(gs6))
risktable xlabel(0(3)15) 
risktable(, order(1 "Slow DunedinPoAm" 2 "Average DunedinPoAm" 3 "Fast DunedinPoAm") failevents title(Number At Risk (Deaths)))
title("")
xtitle(Analysis Time (years))
ytitle(Survival)
legend(pos(3) cols(1) symxsize(5) lab(1 "Slow") lab(2 "Average") lab(3 "Fast")
	title(DunedinPoAm, size(medsmall)) region(lcolor(white)) )
plot1opts(lwidth(medthick) lcolor(midblue))
plot2opts(lwidth(medthick) lcolor(gs10))
plot3opts(lwidth(medthick) lcolor(red))
ylabel(,angle(horiz) nogrid labsize(medlarge) format(%9.2f))
xlabel(,labsize(medlarge))
graphregion(color(white)) plotregion(color(white))
xsize(7) ysize(4)
name(naswcens,replace)
; #delimit cr

graph export "/Users/$user/Box/Belsky/mPoA/eLife/R1/Figures/Figure4R2.pdf", replace	
