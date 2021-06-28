////////Systematic review: baseline glycemic control and cancer survival////////
//////////////////////////////////*SL 25 Mar 2021*//////////////////////////////

////////////step 1: re-construct time-to-event data from survival curve*////////

{
cd "Z:\NewJob\SR_Apr2020\Data\Figure_data"


/*Kaneda-2012 tumour-free-survival HCC*/
import delimited "Kaneda-2012.csv", clear
ipdfc, surv( hba1c65) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) saving(temp0, replace)
import delimited "Kaneda-2012.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) saving(temp1, replace)
use temp0, clear
gen hba1c65 = 0
append using temp1
replace hba1c65=1 if hba1c65==.
sort time hba1c65
stset time, f(event==1)
strate hba1c65
stcox hba1c65
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>6.5"
gen study ="Kaneda"
gen outcome="tumour-free-survival"
gen cancer="HCC"
save res, replace

/*Kang-2016 cancer-specific mortality Upper Tract Urothelial Carcinoma*/
import delimited "Kang-2016_CSM.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) totevents(7) iso saving(temp0, replace)
import delimited "Kang-2016_CSM.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) totevents(23) iso saving(temp1, replace)
use temp0, clear
gen hba1c7 = 0
append using temp1
replace hba1c7=1 if hba1c7==.
sort time hba1c7
stset time, f(event==1) scale(12)
strate hba1c7
stcox hba1c7
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Kang"
gen outcome="cancer-specific mortality"
gen cancer="Urothelial"
append using res
save res, replace

/*Kang-2016 all-cause mortality Upper Tract Urothelial Carcinoma*/
import delimited "Kang-2016_OS.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) totevents(13) iso saving(temp0, replace)
import delimited "Kang-2016_OS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) totevents(32) iso saving(temp1, replace)
use temp0, clear
gen hba1c7 = 0
append using temp1
replace hba1c7=1 if hba1c7==.
sort time hba1c7
stset time, f(event==1) scale(12)
strate hba1c7
stcox hba1c7
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Kang"
gen outcome="all-cause mortality"
gen cancer="Urothelial"
append using res
save res, replace

import delimited "Kang-2016_RFS.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) totevents(5) iso saving(temp0, replace)
import delimited "Kang-2016_RFS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) totevents(17) iso saving(temp1, replace)
use temp0, clear
gen hba1c7 = 0
append using temp1
replace hba1c7=1 if hba1c7==.
sort time hba1c7
stset time, f(event==1) scale(12)
strate hba1c7
stcox hba1c7
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Kang"
gen outcome="recurrence-free survival"
gen cancer="Urothelial"
append using res
save res, replace

import delimited "Tai-2015-DFS.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) totevents(8) iso  probability saving(temp0, replace)
import delimited "Tai-2015-DFS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) iso totevents(20) probability saving(temp1, replace)
use temp0, clear
gen hba1c7 = 0
append using temp1
replace hba1c7=1 if hba1c7==.
sort time hba1c7
stset time, f(event==1)
strate hba1c7
stcox hba1c7
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Tai"
gen outcome="recurrence-free survival"
gen cancer="Bladder"
append using res
save res, replace

import delimited "Lee-2016-DFS.csv", clear
ipdfc, surv(hba1c9) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) iso saving(temp0, replace)
import delimited "Lee-2016-DFS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) iso saving(temp1, replace)
use temp0, clear
gen hba1c9 = 0
append using temp1
replace hba1c9=1 if hba1c9==.
sort time hba1c9
stset time, f(event==1) scale(12)
strate hba1c9
stcox hba1c9
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=9"
gen study ="Lee"
gen outcome="disease-free survival"
gen cancer="Pancreatic"
append using res
save res, replace

import delimited "Lee-2016-OS.csv", clear
ipdfc, surv(hba1c9) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) iso totevents(34) saving(temp0, replace)
import delimited "Lee-2016-OS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) iso totevents(19) saving(temp1, replace)
use temp0, clear
gen hba1c9 = 0
append using temp1
replace hba1c9=1 if hba1c9==.
sort time hba1c9
stset time, f(event==1) scale(12)
strate hba1c9
stcox hba1c9
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=9"
gen study ="Lee"
gen outcome="Overall survival"
gen cancer="Pancreatic"
append using res
save res, replace

import delimited "Okamura-2017_DFS.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) iso  probability saving(temp0, replace)
import delimited "Okamura-2017_DFS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) iso probability saving(temp1, replace)
use temp0, clear
gen hba1c7 = 0
append using temp1
replace hba1c7=1 if hba1c7==.
sort time hba1c7
stset time, f(event==1)
strate hba1c7
stcox hba1c7
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Okamura"
gen outcome="disease-free survival"
gen cancer="Esophagus"
append using res
save res, replace

import delimited "Okamura-2017-OS.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) iso probability saving(temp0, replace)
import delimited "Okamura-2017-OS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) iso probability saving(temp1, replace)
use temp0, clear
gen hba1c7 = 0
append using temp1
replace hba1c7=1 if hba1c7==.
sort time hba1c7
stset time, f(event==1)
strate hba1c7
stcox hba1c7
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Okamura"
gen outcome="Overall survival"
gen cancer="Esophagus"
append using res
save res, replace

/*analysis of raw data Komatsu-2020 does not match what have been reported in the paper*/
/*I decided to adhere to KM plot in the paper instead of supplemental raw data*/
import delimited "Komatsu-2020-OS.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) iso saving(temp0, replace)
import delimited "Komatsu-2020-OS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) iso saving(temp1, replace)
use temp0, clear
gen hba1c7 = 0
append using temp1
replace hba1c7=1 if hba1c7==.
sort time hba1c7
stset time, f(event==1)
strate hba1c7
stcox hba1c7
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Komatsu"
gen outcome="Overall survival"
gen cancer="Lung"
append using res
save res, replace

import delimited "Komatsu-2020-OS-3levels.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) iso saving(temp0, replace)
import delimited "Komatsu-2020-OS-3levels.csv", clear
ipdfc, surv(hba1c78) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) iso saving(temp1, replace)
import delimited "Komatsu-2020-OS-3levels.csv", clear
ipdfc, surv(hba1c8) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) iso saving(temp2, replace)
use temp0, clear
gen exposure = "hba1c<7"
append using temp1
replace exposure="hba1c7-8" if exposure==""
append using temp2
replace exposure="hba1c>8" if exposure==""
sort time exposure
sencode exposure, replace
stset time, f(event==1)
strate exposure
stcox i.exposure
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
drop if parm == "1b.exposure"
gen exposure="3levels-hba1c7-8"  if parm == "2.exposure"
replace exposure="3levels-hba1c>8" if parm == "3.exposure"
keep estimate min95 max95 exposure
gen study ="Komatsu"
gen outcome="Overall survival"
gen cancer="Lung"
append using res
save res, replace

/*Ahn 2016 progression free survival bladder cancer for person-years*/
import delimited "ahn-2016-PFS.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) totevents(7) iso saving(temp0, replace)
import delimited "ahn-2016-PFS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) totevents(19) iso saving(temp1, replace)
use temp0, clear
gen hba1c7 = 0
append using temp1
replace hba1c7=1 if hba1c7==.
sort time hba1c7
stset time, f(event==1) scale(12)
strate hba1c7
stcox hba1c7
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Ahn"
gen outcome="progression free survival"
gen cancer="bladder"
append using res
save res, replace

/*for person years*/
import delimited "cheon-2014-OS.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) totevents(19) probability iso saving(temp0, replace)
import delimited "cheon-2014-OS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) totevents(37) probability iso saving(temp1, replace)
use temp0, clear
gen hba1c7 = 0
append using temp1
replace hba1c7=1 if hba1c7==.
sort time hba1c7
stset time, f(event==1) scale(365.24)
strate hba1c7
stcox hba1c7
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Cheon"
gen outcome="Overall survival"
gen cancer="pancreatic"
append using res
save res, replace

import delimited "Huang-2020-RFS.csv", clear
ipdfc, surv(hba1c7) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) probability iso saving(temp0, replace)
import delimited "Huang-2020-RFS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) probability iso saving(temp1, replace)
use temp0, clear
gen hba1c7 = 0
append using temp1
replace hba1c7=1 if hba1c7==.
sort time hba1c7
stset time, f(event==1) scale(12)
strate hba1c7
stcox hba1c7
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Cheon"
gen outcome="Overall survival"
gen cancer="pancreatic"
append using res
save res, replace

import delimited "kondo2013-OS.csv", clear
ipdfc, surv(hba1c75) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) probability iso saving(temp0, replace)
import delimited "kondo2013-OS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) probability iso saving(temp1, replace)
use temp0, clear
gen hba1c75 = 0
append using temp1
replace hba1c75=1 if hba1c75==.
sort time hba1c75
stset time, f(event==1) scale(365.24)
strate hba1c75
stcox hba1c75
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Cheon"
gen outcome="Overall survival"
gen cancer="pancreatic"
append using res
save res, replace

import delimited "lee-2017-colon-OS.csv", clear
ipdfc, surv(hba1c8) tstart(x) trisk(trisk) nrisk(nrisk0) gen (time event) probability iso saving(temp0, replace)
import delimited "lee-2017-colon-OS.csv", clear
ipdfc, surv(v3) tstart(x) trisk(trisk) nrisk(nrisk1) gen (time event) probability iso saving(temp1, replace)
use temp0, clear
gen hba1c8 = 0
append using temp1
replace hba1c8=1 if hba1c8==.
sort time hba1c8
stset time, f(event==1) scale(12)
strate hba1c8
stcox hba1c8
parmest, fast format(estimate min95 max95 p %8.2f p %8.1e) list(,) eform
keep estimate min95 max95
gen exposure="hba1c>=7"
gen study ="Cheon"
gen outcome="Overall survival"
gen cancer="pancreatic"
append using res
save res, replace
}

////////////step 2: re-clean database and conduct analysis*////////


/*per 1 unit increment*/

cd "Z:\NewJob\SR_Apr2020\Data"
/*analyse raw data for Komatsu T first to get an coefficient*/

{
import excel "Komatsu-2020-rawdata.xls", firstrow clear
keep if Diabetesmellitus=="Yes"
codebook age
tab Sex, m
sencode Sex, replace
tab smoking, m
sencode smoking, replace
tab CVD, m
tab CAD, m
tab Arrhythmia, m
tab CKD, m
tab COPD, m
tab NewStage, m
gen stage=1 if NewStage=="0" | NewStage=="1A" | NewStage=="1B"
replace stage=2 if NewStage=="2A" | NewStage=="2B"
replace stage=3 if NewStage=="3A" | NewStage=="4"
replace stage=9 if NewStage=="" 
tab stage, m
tab death_censor,m
codebook SurvivalYears
stset SurvivalYears, f(death_censor==1)
stcox HbA1c age Sex i.smoking BMI i.stage
}
/*data conversion from cut-off to continuous*/

{
import excel "data.xlsx", sheet("cut-off to per-1-unit") firstrow clear
gen lnrr=ln(RR)
gen lnlb=ln(LCI)
gen lnub=ln(UCI)
gen se=(lnub-lnlb)/3.92

/*for Cheon YK 2014*/
preserve
keep if stauthor=="Cheon YK"
glst lnrr dose, se(se) cov(N n) ci eform
restore

/*for Kaneda K 2012*/
preserve
keep if stauthor=="Kaneda K"
glst lnrr dose, se(se) cov(N n) ci eform
restore

/*for Lee W 2016*/
preserve
keep if stauthor=="Lee W"
keep if outcome=="overall survival/death"
glst lnrr dose, se(se) cov(N n) ci eform
restore
preserve
keep if stauthor=="Lee W"
keep if outcome=="disease-free survival"
glst lnrr dose, se(se) cov(N n) ci eform
restore

/*for Lee SJ 2017*/
preserve
keep if stauthor=="Lee SJ"
glst lnrr dose, se(se) cov(N n) ci eform
restore

/*for Okamura A 2017*/
preserve
keep if stauthor=="Okamura A"
keep if outcome=="overall survival/death"
glst lnrr dose, se(se) cov(N n) ci eform
restore
preserve
keep if stauthor=="Okamura A"
keep if outcome=="cancer specific mortality"
glst lnrr dose, se(se) cov(N n) ci eform
restore

/*for Siddiqui AA 2008*/
preserve
keep if stauthor=="Siddiqui AA"
keep if outcome=="cancer specific mortality"
glst lnrr dose, se(se) cov(N n) cc eform
restore
preserve
keep if stauthor=="Siddiqui AA"
keep if outcome=="overall survival/death"
glst lnrr dose, se(se) cov(N n) cc eform
restore
}
/*import excel cannot read cd path sometimes, don't know why*/
/*per-1-unit data clean*/
{

/*import excel "data.xlsx", sheet("table1") firstrow clear
tempfile overall
save `overall', replace
*/
import excel "data.xlsx", sheet("per-1-unit") firstrow clear
tostring year, replace
gen Study = stauthor + " et al (" + year + ")"
order Study, first
drop stauthor year journal 
distinct Study
tab outcome, m
replace outcome = "All-cause mortality" if outcome == "overall survival/death"
replace outcome = "Disease-free survival" if outcome == "disease-free survival"
replace outcome = "Cancer-specific mortality" if outcome == "cancer specific mortality"
replace outcome = "Recurrence" if outcome == "recurrence"
bysort outcome: gen studyno=_N
drop if studyno<2
gsort -studyno
tab outcome
label variable N "Participants"
label variable n "Events"
replace site = proper(site)
label variable site "Cancer site"
tostring n, replace
replace n="NR" if n=="."
sort outcome
sencode outcome, replace
tab outcome, m
gen lnrr=ln(RR)
gen lnlb=ln(LCI)
gen lnub=ln(UCI)
gen se=(lnub-lnlb)/3.92
set scheme s2mono
save db2, replace
}

/*binary exposure*/
cd "Z:\NewJob\SR_Apr2020\Data"
{

/// calculate mean of HbA1c by a cut off 7% then calculate HR for >7% vs. <7%
/*
///Boursi bladder use Wang 2014 excel calculator to calculate mean/SD from median/IQR
cutpconv, cutoff(7) mean(7.166666667) sd(1.187563442)

///Boursi breast use Wang 2014 excel calculator to calculate mean/SD from median/IQR
cutpconv, cutoff(7) mean(7.166666667) sd(1.187333863)

///Boursi colorectal use Wang 2014 excel calculator to calculate mean/SD from median/IQR
cutpconv, cutoff(7) mean(7.1) sd(1.112828408)

///Boursi pancreatic use Wang 2014 excel calculator to calculate mean/SD from median/IQR
cutpconv, cutoff(7) mean(8.233333333) sd(2.154723809)

///Boursi prostate use Wang 2014 excel calculator to calculate mean/SD from median/IQR
cutpconv, cutoff(7) mean(6.966666667) sd(1.038580727)
*/


import excel "data.xlsx", sheet("per-1-unit to 7% cut-off") firstrow clear
gen lnrr1=ln(RR)
gen lnlb=ln(LCI)
gen lnub=ln(UCI)
gen se1=(lnub-lnlb)/3.92

gen lnrr = lnrr1 * (mean2 - mean1) /*calculate HR for >7% vs. <7%*/
gen se = se1 * (mean2 - mean1) /*calculate SE for >7% vs. <7%*/
}

{
import excel "data.xlsx", sheet("poorly-controlled") firstrow clear
tostring year, replace
gen Study = stauthor + " et al (" + year + ")"
order Study, first
drop stauthor year journal 
sort outcome site
tab outcome, m
replace outcome = "All-cause mortality" if outcome == "overall survival/death"
replace outcome = "Disease-free survival" if outcome == "disease-free survival"
replace outcome = "Cancer-specific mortality" if outcome == "cancer specific mortality"
replace outcome = "Recurrence" if outcome == "recurrence"
bysort outcome: gen studyno=_N
drop if studyno<2
gsort -studyno
foreach var in N1 N2 n1 n2 {
	replace `var'=0 if `var'==.
}
gen Totpar = N1 + N2
gen totevents = n1+n2
replace Totpar = N if Totpar==0
/*these studies reported total number of participants & events but not in each group*/

replace Totpar=1850 if Study=="Boursi B et al (2016)" & site=="colorectal"
replace totevents=769 if Study=="Boursi B et al (2016)" & site=="colorectal"

replace Totpar=1382 if Study=="Boursi B et al (2016)" & site=="breast"
replace totevents=332 if Study=="Boursi B et al (2016)" & site=="breast"

replace Totpar=1168 if Study=="Boursi B et al (2016)" & site=="bladder"
replace totevents=403 if Study=="Boursi B et al (2016)" & site=="bladder"

replace Totpar=634 if Study=="Boursi B et al (2016)" & site=="pancreatic"
replace totevents=519 if Study=="Boursi B et al (2016)" & site=="pancreatic"

replace Totpar=1994 if Study=="Boursi B et al (2016)" & site=="prostate"
replace totevents=564 if Study=="Boursi B et al (2016)" & site=="prostate"

replace Totpar=417 if Study=="Lee et al (2015)"


label variable Totpar "Participants"
label variable totevents "Events"
replace site = proper(site)
label variable site "Cancer site"
tostring totevents, replace
replace totevents="NR" if totevents=="0"
sort outcome
sencode outcome, replace
tab outcome, m
 
replace lnrr=ln(RR) if lnrr==.
gen lnlb=ln(LCI)
gen lnub=ln(UCI)
replace se=(lnub-lnlb)/3.92 if se==.
set scheme s2mono
save db1, replace
}


/*main analysis: forest plot*/
set scheme s2mono
use db1, clear
sort outcome Study site
metan lnrr se, by(outcome) sortby(lnrr) random eform nooverall lcols(Study Totpar totevents) sgweight boxsca(30) xlabel (0.5, 1.0, 3.0, 5.0, 10.0) texts(100) olineopt(lwidth(thin)) diamopt(lcolor(black)) ciopt(lwidth(vthin)) pointopt(msymbol(O)) effect(RR) nobox title("HbA1c ≥ 7% vs. HbA1c < 7%", size(small))saving("poorly", replace)

use db2, clear
sort outcome Study site
metan lnrr se, by(outcome) sortby(lnrr) random eform nooverall lcols(Study N n) sgweight boxsca(30) xlabel ( 0.8, 1.0, 1.2, 1.5) texts(100) olineopt(lwidth(thin)) diamopt(lcolor(black)) ciopt(lwidth(vthin)) pointopt(msymbol(O)) effect(RR) nobox title("HbA1c per-1-unit increment", size(small))saving ("per", replace)
graph combine "poorly.gph" "per.gph", col(2) iscale(0.5) imargins(zero) xsize(11.7) ysize(8.3) scale(1) graphregion(margin(small)) plotregion(margin(medium))
graph export "Z:/NewJob/SR_Apr2020/manuscript/V2/F2.svg", as(svg) replace

/*main analysis: funnel plots*/
/*poorly controlled vs. well controlled HbA1c*/

use db1, clear
preserve
keep if outcome==1
metabias lnrr se , egger
metafunnel lnrr se , ytitle("SE of log(RR)") subtitle(All-cause mortality, Egger's test (p=0.005)) saving("poorly1", replace) nodraw
meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
meta trimfill, left eform(eformstring)
meta trimfill, funnel(ytitle(SE of log(RR)) xtitle(log(RR)) title(All-cause mortality, size(large))) random(dlaird)
graph save "trim1.gph", replace
graph close
restore

preserve
keep if outcome ==2
metabias lnrr se, egger
metafunnel lnrr se, xtitle("") ytitle("SE of log(RR)") subtitle(Cancer-specific mortality, Egger's test (p=0.017)) saving("poorly2", replace) nodraw
meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
meta trimfill, left eform(eformstring)
meta trimfill, funnel(ytitle(SE of log(RR)) xtitle(log(RR)) title(Cancer-specific mortality, size(large))) random(dlaird)
graph save "trim2.gph", replace
graph close
restore

preserve
keep if outcome ==3
metabias lnrr se, egger
metafunnel lnrr se, xtitle("log(RR)") ytitle("SE of log(RR)") subtitle(Recurrence, Egger's test (p=0.005)) saving("poorly3", replace) nodraw
meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
meta trimfill, left eform(eformstring)
meta trimfill, funnel(ytitle(SE of log(RR)) xtitle(log(RR)) title(Recurrence, size(large))) random(dlaird)
graph save "trim3.gph", replace
graph close
restore
/* 
metabias lnrr se if outcome==4, egger
metafunnel lnrr se if outcome==4, xtitle("") ytitle("SE of log(RR)") subtitle(Disease-free Survival, Egger's test (p=0.038)) saving("poorly4", replace) nodraw
*/
use db2, clear
rename N Totpar
preserve
keep if outcome==1
metabias lnrr se, egger
metafunnel lnrr se, xtitle("") ytitle("") subtitle(All-cause mortality, Egger's test (p=0.040)) saving("per1", replace) nodraw
meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
meta trimfill, left eform(eformstring)
meta trimfill, funnel(ytitle(SE of log(RR)) xtitle(log(RR)) title(All-cause mortality, size(large))) random(dlaird)
graph save "trim1.gph", replace
restore

preserve
keep if outcome ==2
metabias lnrr se , egger
metafunnel lnrr se , xtitle("") ytitle("") subtitle(Cancer-specific mortality, Egger's test (p=0.018)) saving("per2", replace) nodraw
meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
meta trimfill, left eform(eformstring)
meta trimfill, funnel(ytitle(SE of log(RR)) xtitle(log(RR)) title(Cancer-specific mortality, size(large))) random(dlaird)
restore

preserve
keep if outcome ==3
metabias lnrr se , egger
metafunnel lnrr se , xtitle("log(RR)") ytitle("") subtitle(Recurrence, Egger's test (p=0.204)) saving("per3", replace) nodraw
meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
meta trimfill, left eform(eformstring)
meta trimfill, funnel(ytitle(SE of log(RR)) xtitle(log(RR)) title(Recurrence, size(large))) random(dlaird)
restore
/*
metabias lnrr se if outcome==4, egger
metafunnel lnrr se if outcome==4, xtitle("log(RR)") ytitle("") subtitle(Disease-free Survival, Egger's test (p=0.068)) saving("per4", replace) nodraw
*/
graph combine "poorly1.gph" "per1.gph" "poorly2.gph" "per2.gph" "poorly3.gph" "per3.gph", col(2) iscale(0.5) imargins(zero) xsize(8.3) ysize(11.7) scale(1) graphregion(margin(small)) plotregion(margin(medium))
graph export "Z:/NewJob/SR_Apr2020/manuscript/V2/F3.svg", as(svg)
/*
/*keep only HR*/
use db1, clear
keep if estimates=="HR"
metan lnrr se, by(outcome) sortby(lnrr) random eform nooverall lcols(Study Totpar totevents) sgweight boxsca(30) xlabel (0.2, 0.5, 1.0, 3.0, 5.0, 10.0) texts(100) olineopt(lwidth(thin)) diamopt(lcolor(black)) ciopt(lwidth(vthin)) pointopt(msymbol(O)) effect(RR) nobox title("poorly controlled vs. well controlled HbA1c", size(small)) saving("hr1", replace)

use db2, clear
keep if estimates=="HR"
metan lnrr se, by(outcome) sortby(lnrr) random eform nooverall lcols(Study N n) sgweight boxsca(30) xlabel ( 0.8, 1.0, 1.2, 1.5) texts(100) olineopt(lwidth(thin)) diamopt(lcolor(black)) ciopt(lwidth(vthin)) pointopt(msymbol(O)) effect(RR) nobox title("HbA1c per-1-unit increment", size(small))saving ("hr2", replace)

graph combine "hr1.gph" "hr2.gph", imargins(zero)
graph export "hr.svg", as(svg) replace
*/

/*subgroup analyses*/
/*trim and fill analysis*/
use db1, clear
forvalues i=1(1)3 {
	preserve
	keep if outcome == `i'
	display `i'
	meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
	meta trimfill, left eform(eformstring)
	meta trimfill, left funnel
	restore
}

use db2, clear
forvalues i=1(1)3 {
	preserve
	keep if outcome == `i'
	display `i'
	meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
	meta trimfill, left eform(eformstring)
	meta trimfill, left funnel
	restore
}
/*by continent & meta-regression with study sample size*/
use db1, clear
gen continent="Western" if Country=="UK" | Country =="USA"
replace continent="Asia" if continent==""
gsort -continent
sencode continent, replace

/*new stata meta-analysis commands*/
/* cannot combine graphs from offcial stata meta commands*/
keep if outcome==1
replace site="Colorectal" if site=="Colon"
meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
meta forestplot, random(dlaird) subgroup(continent) eform(RR) cibind(parentheses) xline(1) title("All-cause mortality")
meta forestplot, random(dlaird) subgroup(site) eform(RR) cibind(parentheses) xline(1) title("All-cause mortality")

meta regress continent, random(dlaird)
/*
meta regress _meta_studysize, random(dlaird)
estat bubbleplot, title("All-cause mortality, categorised HbA1c (p<0.001)", size(small)) ytitle("log(RR)")
*/

use db1, clear
gen continent="Western" if Country=="UK" | Country =="USA"
replace continent="Asia" if continent==""
gsort -continent
sencode continent, replace
keep if outcome==2
meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
meta forestplot, random(dlaird) subgroup(continent) eform(RR) cibind(parentheses) xline(1) title("Cancer-specific mortality")
meta regress continent, random(dlaird)

/*
meta regress _meta_studysize
estat bubbleplot, title("Recurrence, categorised HbA1c (p=0.315)", size(small)) ytitle("log(RR)")
*/
use db1, clear
gen continent="Western" if Country=="UK" | Country =="USA"
replace continent="Asia" if continent==""
gsort -continent
sencode continent, replace
keep if outcome==3
meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
meta forestplot, random(dlaird) subgroup(continent) eform(RR) cibind(parentheses) xline(1) title("Recurrence")
meta regress _meta_studysize
estat bubbleplot, title("Disease-free survival, categorised HbA1c (p=0.003)", size(small)) ytitle("log(RR)")

use db2, clear
gen continent="Western" if country=="UK" | country =="USA"
replace continent="Asia" if continent==""
gsort -continent
sencode continent, replace

keep if outcome==1
meta set lnrr se, studylabel(Study) studysize(N) random(dlaird)
meta forestplot, random(dlaird) subgroup(continent) eform(RR) cibind(parentheses) xline(1) title("All-cause mortality")
meta regress _meta_studysize, random(dlaird)
estat bubbleplot, title("All-cause mortality, per-1-unit increment in HbA1c (p=0.028)", size(small)) ytitle("log(RR)")

use db2, clear
gen continent="Western" if country=="UK" | country =="USA"
replace continent="Asia" if continent==""
gsort -continent
sencode continent, replace

keep if outcome==2
meta set lnrr se, studylabel(Study) studysize(N) random(dlaird)
meta forestplot, random(dlaird) subgroup(continent) eform(RR) cibind(parentheses) xline(1) title("Cancer-specific mortality")
meta regress _meta_studysize
estat bubbleplot, title("Recurrence, per-1-unit increment in  HbA1c (p=0.590)", size(small)) ytitle("log(RR)")

use db2, clear
gen continent="Western" if country=="UK" | country =="USA"
replace continent="Asia" if continent==""
gsort -continent
sencode continent, replace

keep if outcome==3
meta set lnrr se, studylabel(Study) studysize(N) random(dlaird)
meta regress _meta_studysize
meta forestplot, random(dlaird) subgroup(continent) eform(RR) cibind(parentheses) xline(1) title("Recurrence")
/*
estat bubbleplot, title("Disease-free survival, per-1-unit increment in  HbA1c (p=0.452)", size(small)) ytitle("log(RR)")
*/
/*by cancer site*/
use db1, clear
drop N
replace site="Colorectal" if site=="Colon"
sdecode outcome, replace
gen out = outcome + ": " + site
bysort out: gen N = _N
drop if N==1

metan lnrr se, by(out) sortby(lnrr) random eform nooverall lcols(Study Totpar totevents) sgweight boxsca(30) xlabel (0.2, 0.5, 1.0, 3.0, 5.0, 10.0) texts(100) olineopt(lwidth(thin)) diamopt(lcolor(black)) ciopt(lwidth(vthin)) pointopt(msymbol(O)) effect(RR) nobox title("HbA1c ≥ 7% vs. HbA1c < 7%", size(small)) saving("site1", replace)
keep if site=="Bladder" & outcome=="Recurrence"
	meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
	meta trimfill, left eform(eformstring)
meta trimfill, funnel(ytitle(SE of log(RR)) xtitle(log(RR)) title(Recurrence: Bladder, size(large))) random(dlaird)
keep if outcome==1
sencode site, replace
	meta set lnrr se, studylabel(Study) studysize(Totpar) random(dlaird)
meta regress i.site

use db2, clear
replace site="Colorectal" if site=="Colon"
sdecode outcome, replace
gen out = outcome + ": " + site
bysort out: gen total = _N
drop if total==1
metan lnrr se, by(out) sortby(lnrr) random eform nooverall lcols(Study N n) sgweight boxsca(30) xlabel ( 0.8, 1.0, 1.2, 1.5) texts(100) olineopt(lwidth(thin)) diamopt(lcolor(black)) ciopt(lwidth(vthin)) pointopt(msymbol(O)) effect(RR) nobox title("HbA1c per-1-unit increment", size(small)) saving("site2", replace)
graph combine "site1.gph" "site2.gph", col(2) iscale(0.5) imargins(zero) xsize(11.7) ysize(8.3) scale(1) graphregion(margin(small)) plotregion(margin(medium))
graph export "Z:/NewJob/SR_Apr2020/manuscript/V2/sF3.svg", as(svg) replace