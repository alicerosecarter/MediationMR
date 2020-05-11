

global PCs PC*
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth


********************************************************************************

foreach out in sbp hyperten CVD_inc {

	putexcel set 202002_results_ldl, sheet(`out') modify
	putexcel A1="Method" B1="Scale" C1="Mediator" D1="Outcome" E1="Total Effect" F1="Total - LCI" G1="Total - UCI" H1="Direct Effect" I1="Direct - LCI" J1="Direct - UCI" K1="Indirect - difference" ///
		L1="LCI" M1="UCI" N1="% mediated - difference" O1="LCI" P1="UCI" Q1="Indirect - product" R1="LCI" S1="UCI" T1="% mediated - product" U1="LCI" V1="UCI" ///
		W1="Exposure - Mediator" X1="Exp-Med - LCI" Y1="Exp-Med - UCI" Z1="Mediator - Outcome" AA1="Med-Out - LCI" AB1="Med-Out - UCI" 
}

******************** Total Effects - Observational on RD scale *****************

foreach out in sbp hyperten CVD_inc {
	
		local x=1	
	
		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

		regress `out' eduyears_scaled $max_covar
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="Observational" B`x'="Risk Difference" C`x'="LDL" D`x'="`out_label'" E`x'=_b[eduyears_scaled] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	


******************* Total Effects - Observational on logOR scale ***************

foreach out in hyperten CVD_inc {
	
		local x=3
	
	foreach med in zldl  {

		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

		logit `out' eduyears_scaled $max_covar

	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="Observational" B`x'="Log OR" C`x'="LDL" D`x'="`out_label'" E`x'=_b[eduyears_scaled] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}

******************* Total Effects - Observational on OR scale ***************

foreach out in hyperten CVD_inc {
	
		local x=5
	
	foreach med in zldl  {

		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

		logit `out' eduyears_scaled $max_covar

	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local beta = exp(_b[eduyears_scaled])
	local ll_2sls1 = exp(_b[eduyears_scaled])-1.96*(_se[eduyears_scaled])
	local ul_2sls1 = exp(_b[eduyears_scaled])+1.96*(_se[eduyears_scaled])
	
	putexcel A`x'="Observational" B`x'="Log OR" C`x'="LDL" D`x'="`out_label'" E`x'=`beta' F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}

******************** Total Effects - One-sample MR on RD scale *****************

foreach out in sbp hyperten CVD_inc {
	
		local x=2
		
	foreach med in zldl  {

		local x=`x'+1
		
	putexcel set 202002_results_ldl, sheet(`out') modify

		ivreg2 `out' (eduyears_scaled = ea_weighted) $PCs $max_covar
	
	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="MR - One-sample" B`x'="Risk Difference" C`x'="LDL" D`x'="`out_label'" E`x'=_b[eduyears_scaled] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}


******************* Total Effects - One-sample MR on logOR scale ***************
capture drop education_MR
regress eduyears_scaled ea_weighted $max_covar $PCs
predict education_MR, xb

foreach out in  hyperten CVD_inc {
	
		local x=4
	
	foreach med in zldl  {

		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

		logit `out' education_MR $PCs $max_covar
	
	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[education_MR]-1.96*_se[education_MR]
	local ul_2sls1 = _b[education_MR]+1.96*_se[education_MR]
	
	putexcel A`x'="MR - One-sample" B`x'="Log OR" C`x'="LDL" D`x'="`out_label'" E`x'=_b[education_MR] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}

******************* Total Effects - One-sample MR on OR scale ***************
capture drop education_MR
regress eduyears_scaled ea_weighted $max_covar $PCs
predict education_MR, xb

foreach out in  hyperten CVD_inc {
	
		local x=6
	
	foreach med in zldl  {

		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

		logistic `out' education_MR $PCs $max_covar
	
	
	matrix results = r(table)
	local out_label : var label `out'
	local med_label : var label `med'
	local beta = results[1,1]
	local ll_2sls1 = results[5,1]
	local ul_2sls1 = results[6,1]
	
	putexcel A`x'="MR - One-sample" B`x'="OR" C`x'="LDL" D`x'="`out_label'" E`x'=`beta' F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}

************* Direct & Indirect Effects - Observational on RD scale ************

foreach out in sbp hyperten CVD_inc {
	
		local x=1	
	
	foreach med in zldl {
		
		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

	regress `out' eduyears_scaled `med' $max_covar
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[`med']-1.96*_se[`med']
	local ul_med = _b[`med']+1.96*_se[`med']
	
	putexcel H`x'=_b[eduyears_scaled] I`x'=`ll_2sls1' J`x'= `ul_2sls1'
	putexcel Z`x'=_b[`med'] AA`x'=`ll_med' AB`x'= `ul_med'

	}	
}


**************** Direct & Indirect Effects - Observational on logOR scale ******

foreach out in  hyperten CVD_inc {
	
		local x=3
	
	foreach med in zldl {

		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

	logit `out' eduyears_scaled `med' $max_covar
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[`med']-1.96*_se[`med']
	local ul_med = _b[`med']+1.96*_se[`med']
	
	putexcel H`x'=_b[eduyears_scaled] I`x'=`ll_2sls1' J`x'= `ul_2sls1'
	putexcel Z`x'=_b[`med'] AA`x'=`ll_med' AB`x'= `ul_med'
	
	}	
}

**************** Direct & Indirect Effects - Observational on OR scale ******

foreach out in  hyperten CVD_inc {
	
		local x=5
	
	foreach med in zldl {

		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

	logistic `out' eduyears_scaled `med' $max_covar

	matrix results = r(table)
	local out_label : var label `out'
	local med_label : var label `med'
	local direct = results[1,1]
	local ll_2sls1 = results[5,1]
	local ul_2sls1 = results[6,1]
	
	local indirect = results[1,2]
	local ll_med = results[5,2]
	local ul_med = results[6,2]
	
	putexcel H`x'=_b[eduyears_scaled] I`x'=`ll_2sls1' J`x'= `ul_2sls1'
	putexcel Z`x'=_b[`med'] AA`x'=`ll_med' AB`x'= `ul_med'
	
	}	
}

*************** Direct & Indirect Effects - One-sample MR on RD scale **********
							
										*LDL*
foreach out in sbp hyperten CVD_inc {
	
		local x=2

		local x=`x'+1
		
	putexcel set 202002_results_ldl, sheet(`out') modify

	ivreg2 `out' (eduyears_scaled zldl = ea_weighted ldl_weighted) $PCs $max_covar
	
	
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[zldl]-1.96*_se[zldl]
	local ul_med = _b[zldl]+1.96*_se[zldl]
	
	putexcel H`x'=_b[eduyears_scaled] I`x'=`ll_2sls1' J`x'= `ul_2sls1'
	putexcel Z`x'=_b[zldl] AA`x'=`ll_med' AB`x'= `ul_med'

	}

*********** Direct & Indirect Effects - One-sample MR on log OR scale **********

							*LDL*

capture drop dir_ealdl
capture drop dir_ldlea
regress eduyears_scaled ea_weighted ldl_weighted $max_covar $PCs
predict dir_ealdl, xb
regress zldl ldl_weighted ea_weighted $max_covar $PCs
predict dir_ldlea, xb


foreach out in  hyperten  CVD_inc {

	putexcel set 202002_results_ldl, sheet(`out') modify
	
	local x=4
	
	logit `out' dir_ealdl dir_ldlea $max_covar $PCs, vce(robust) // Natural Direct effect of Education
	
	local out_label : var label `out'
	local ll_2sls1 = _b[dir_ealdl]-1.96*_se[dir_ealdl]
	local ul_2sls1 = _b[dir_ealdl]+1.96*_se[dir_ealdl]
	local ll_med = _b[dir_ldlea]-1.96*_se[dir_ldlea]
	local ul_med = _b[dir_ldlea]+1.96*_se[dir_ldlea]
	
	local x=`x'+1	
		
	putexcel H`x'=_b[dir_ealdl] I`x'=`ll_2sls1' J`x'= `ul_2sls1'
	putexcel Z`x'=_b[dir_ldlea] AA`x'=`ll_med' AB`x'= `ul_med'
	
	
}	


*********** Direct & Indirect Effects - One-sample MR on OR scale **********

							*LDL*

capture drop dir_ealdl
capture drop dir_ldlea
regress eduyears_scaled ea_weighted ldl_weighted $max_covar $PCs
predict dir_ealdl, xb
regress zldl ldl_weighted ea_weighted $max_covar $PCs
predict dir_ldlea, xb


foreach out in  hyperten  CVD_inc {

	local x=6

	putexcel set 202002_results_ldl, sheet(`out') modify
	
	logistic `out' dir_ealdl dir_ldlea $max_covar $PCs, vce(robust) // Natural Direct effect of Education
	
	local x=`x'+1
	
	matrix results = r(table)
	local direct = results[1,1]
	local ll_2sls1 = results[5,1]
	local ul_2sls1 = results[6,1]
	
	local indirect = results[1,2]
	local ll_med = results[5,2]
	local ul_med = results[6,2]
	
	putexcel H`x'=`direct' I`x'=`ll_2sls1' J`x'= `ul_2sls1'
	putexcel Z`x'=`indirect' AA`x'=`ll_med' AB`x'= `ul_med'
	
	
}	

	
********* Total Effects Exposure-Mediator - Observational on RD scale **********

foreach out in sbp hyperten CVD_inc {
	
		local x=1	
	
	foreach med in zldl {
		
		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

		regress `med' eduyears_scaled $max_covar
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel W`x'=_b[eduyears_scaled] X`x'=`ll_2sls1' Y`x'=`ul_2sls1'
	
	}	
}

******** Total Effects Exposure-Mediator - Observational on log OR scale *******

/* All mediators are continuous variables and as such the exposure-mediator
	pathway cannot be estimated on the log OR scale */

foreach out in  hyperten CVD_inc {
	
		local x=3	
	
	foreach med in zldl {
		
		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

		regress `med' eduyears_scaled $max_covar
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel W`x'=_b[eduyears_scaled] X`x'=`ll_2sls1' Y`x'=`ul_2sls1'
	
	}	
}	
	
******** Total Effects Exposure-Mediator - Observational on  OR scale *******

/* All mediators are continuous variables and as such the exposure-mediator
	pathway cannot be estimated on the log OR scale */

foreach out in  hyperten CVD_inc {
	
		local x=5	
	
	foreach med in zldl {
		
		local x=`x'+1
	
	putexcel set 202002_results_ldl, sheet(`out') modify

		regress `med' eduyears_scaled $max_covar
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel W`x'=_b[eduyears_scaled] X`x'=`ll_2sls1' Y`x'=`ul_2sls1'
	
	}	
}	
	
********** Total Effects Exposure-Mediator - One-sample MR on RD scale *********

foreach out in sbp hyperten CVD_inc {
	
		local x=2
		
	foreach med in zldl {

		local x=`x'+1
		
	putexcel set 202002_results_ldl, sheet(`out') modify

		ivreg2 `med' (eduyears_scaled = ea_weighted) $PCs $max_covar
	
		
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel W`x'=_b[eduyears_scaled] X`x'=`ll_2sls1' Y`x'=`ul_2sls1'
	
	}	
}

******** Total Effects Exposure-Mediator - One-sample MR on log OR scale *******

/* All mediators are continuous variables and as such the exposure-mediator
	pathway cannot be estimated on the log OR scale */

foreach out in  hyperten CVD_inc {
	
		local x=4
		
	foreach med in zldl {

		local x=`x'+1
		
	putexcel set 202002_results_ldl, sheet(`out') modify

		ivreg2 `med' (eduyears_scaled = ea_weighted) $PCs $max_covar
	
		
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel W`x'=_b[eduyears_scaled] X`x'=`ll_2sls1' Y`x'=`ul_2sls1'
	
	}	
}	

******** Total Effects Exposure-Mediator - One-sample MR on OR scale *******

/* All mediators are continuous variables and as such the exposure-mediator
	pathway cannot be estimated on the log OR scale */

foreach out in  hyperten CVD_inc {
	
		local x=6
		
	foreach med in zldl {

		local x=`x'+1
		
	putexcel set 202002_results_ldl, sheet(`out') modify

		ivreg2 `med' (eduyears_scaled = ea_weighted) $PCs $max_covar
	
		
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel W`x'=_b[eduyears_scaled] X`x'=`ll_2sls1' Y`x'=`ul_2sls1'
	
	}	
}	


global PCs PC*
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth
global med zldl 
rename eduyears_scaled eduyears_3_6
egen eduyears_scaled = std(eduyears_3_6)
set more off


**************************************************
* Observational mediators on the Risk diff scale *
**************************************************

capture program drop boot
program boot, rclass
foreach med in zldl { 
        
     if (`med' == zldl) {
        regress zldl eduyears_scaled $max_covar
        scalar exp_zldl = _b[eduyears_scaled]
        }


foreach out in sbp hyperten CVD_inc {     

     
     if (`out' == sbp) {
	reg sbp eduyears_scaled $max_covar
	scalar total_sbp = _b[eduyears_scaled]
        regress sbp `med' eduyears_scaled $max_covar 
        scalar med_sbp = _b[`med']
		scalar direct_sbp = _b[eduyears_scaled]
		scalar indirect_diff_`med'_sbp = total_sbp-direct_sbp
        scalar indirect_`med'_sbp = exp_`med'*med_sbp
	scalar prop_`med'_sbp = (indirect_`med'_sbp/total_sbp)*100
	scalar diff_prop_`med'_sbp = (indirect_diff_`med'_sbp/total_sbp)*100
        }
        if (`out' == hyperten) {
	reg hyperten eduyears_scaled $max_covar
	scalar total_hyperten = _b[eduyears_scaled]
        regress hyperten `med' eduyears_scaled $max_covar 
        scalar med_hyperten = _b[`med']
		scalar direct_hyperten = _b[eduyears_scaled]
		scalar indirect_diff_`med'_hyperten = total_hyperten-direct_hyperten
        scalar indirect_`med'_hyperten = exp_`med'*med_hyperten
	scalar prop_`med'_hyperten = (indirect_`med'_hyperten/total_hyperten)*100
	scalar diff_prop_`med'_hyperten = (indirect_diff_`med'_hyperten/total_hyperten)*100
        }
		if (`out' == CVD_inc) {
	reg CVD_inc eduyears_scaled $max_covar
	scalar total_CVD = _b[eduyears_scaled]
        reg CVD_inc `med' eduyears_scaled $max_covar 
        scalar med_CVD = _b[`med']
		scalar direct_CVD = _b[eduyears_scaled]
		scalar indirect_diff_`med'_CVD_inc = total_CVD-direct_CVD
        scalar indirect_`med'_CVD_inc = exp_`med'*med_CVD
	scalar prop_`med'_CVD_inc = (indirect_`med'_CVD_inc/total_CVD)*100
	scalar diff_prop_`med'_CVD_inc = (indirect_diff_`med'_CVD_inc/total_CVD)*100
        }
 }
}
 end     

 foreach out in sbp hyperten  CVD_inc {
 local x=1
 local y=1
 local z=1
 local a=1
 foreach med in zldl{ 
     putexcel set 202002_results_ldl, sheet(`out') modify
       reg `med'  eduyears_scaled $max_covar
         reg `out' `med' eduyears_scaled $max_covar 
 bootstrap indirect = indirect_`med'_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]
         local x=`x'+1
		 putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'
		 
		reg `out'  eduyears_scaled $max_covar
		reg `out'  eduyears_scaled `med' $max_covar
		 bootstrap indirect_diff = indirect_diff_`med'_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect_diff = results[1,1]
         local indirect_LCI_diff = results[5,1]
         local indirect_UCI_diff = results[6,1]
         local a=`a'+1
     	 putexcel K`a'=`indirect_diff' L`a'=`indirect_LCI_diff' M`a'= `indirect_UCI_diff'

  	 
	 reg `out' eduyears_scaled $max_covar
         reg `med'  eduyears_scaled $max_covar
         reg `out' `med' eduyears_scaled $max_covar 

bootstrap prop_med = prop_`med'_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
         local y=`y'+1

      	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI' 
		
	reg `out'  eduyears_scaled $max_covar
	reg `out'  eduyears_scaled `med' $max_covar
bootstrap prop_med_diff = diff_prop_`med'_`out', reps(100) nodrop : boot 
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local z=`z'+1
	 
   		putexcel N`z'=`proportion' O`z'=`proportion_LCI' P`z'=`proportion_UCI'
 }
}


***********************************************
* Observational mediators on the log OR scale *
***********************************************

capture program drop boot
program boot, rclass
foreach med in zldl { 
        
     if (`med' == zldl) {
        regress zldl eduyears_scaled $max_covar
        scalar exp_zldl = _b[eduyears_scaled]
        }

foreach out in  hyperten CVD_inc {     
        
		if (`out' == hyperten) {
	logit hyperten eduyears_scaled $max_covar
	scalar total_hyperten = _b[eduyears_scaled]
        logit hyperten `med' eduyears_scaled $max_covar 
        scalar med_hyperten = _b[`med']
		scalar direct_hyperten = _b[eduyears_scaled]
		scalar indirect_diff_`med'_hyperten = total_hyperten-direct_hyperten
        scalar indirect_`med'_hyperten = exp_`med'*med_hyperten
	scalar prop_`med'_hyperten = (indirect_`med'_hyperten/total_hyperten)*100
	scalar diff_prop_`med'_hyperten = (indirect_diff_`med'_hyperten/total_hyperten)*100
  }
     if (`out' == CVD_inc) {
	logit CVD_inc eduyears_scaled $max_covar
	scalar total_CVD = _b[eduyears_scaled]
        logit CVD_inc `med' eduyears_scaled $max_covar 
        scalar med_CVD = _b[`med']
		scalar direct_CVD = _b[eduyears_scaled]
		scalar indirect_diff_`med'_CVD_inc = total_CVD-direct_CVD
        scalar indirect_`med'_CVD_inc = exp_`med'*med_CVD
	scalar prop_`med'_CVD_inc = (indirect_`med'_CVD_inc/total_CVD)*100
	scalar diff_prop_`med'_CVD_inc = (indirect_diff_`med'_CVD_inc/total_CVD)*100
        }
  }
}
 end     

 foreach out in hyperten CVD_inc {
 local x=3
 local y=3
 local z=3
 local a=3
 foreach med in zldl { 
 
     putexcel set 202002_results_ldl, sheet(`out') modify
		reg `med'  eduyears_scaled $max_covar
         logit `out' `med' eduyears_scaled $max_covar 
 bootstrap indirect = indirect_`med'_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]
 local x=`x'+1

		putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI' 
  
  		logit `out'  eduyears_scaled $max_covar
		logit `out'  eduyears_scaled `med' $max_covar
		 bootstrap indirect_diff = indirect_diff_`med'_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect_diff = results[1,1]
         local indirect_LCI_diff = results[5,1]
         local indirect_UCI_diff = results[6,1]
         local a=`a'+1
     	 putexcel K`a'=`indirect_diff' L`a'=`indirect_LCI_diff' M`a'= `indirect_UCI_diff'
  
		 logit `out' eduyears_scaled $max_covar
         reg `med'  eduyears_scaled $max_covar
         logit `out' `med' eduyears_scaled $max_covar  	
bootstrap prop_med = prop_`med'_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local y = `y'+1 
 
      	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI'

	logit `out'  eduyears_scaled $max_covar
	logit `out'  eduyears_scaled `med' $max_covar
bootstrap diff_prop_`med'_`out', reps(100) nodrop : boot 
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local z = `z'+1
	 
   		putexcel N`z'=`proportion' O`z'=`proportion_LCI' P`z'=`proportion_UCI'


 }
}


***********************************************
* Observational mediators on the OR scale *
***********************************************

capture program drop boot
program boot, rclass
foreach med in zldl { 
        
     if (`med' == zldl) {
        regress zldl eduyears_scaled $max_covar
        scalar exp_zldl = _b[eduyears_scaled]
        }

foreach out in  hyperten CVD_inc {     
        
		if (`out' == hyperten) {
	logistic hyperten eduyears_scaled $max_covar
	scalar total_hyperten = exp(_b[eduyears_scaled])
        logistic hyperten `med' eduyears_scaled $max_covar 
        scalar med_hyperten = exp(_b[`med'])
		scalar direct_hyperten = exp(_b[eduyears_scaled])
		scalar indirect_diff_`med'_hyperten = total_hyperten-direct_hyperten
        scalar indirect_`med'_hyperten = exp_`med'*med_hyperten
	scalar prop_`med'_hyperten = (indirect_`med'_hyperten/total_hyperten)*100
	scalar diff_prop_`med'_hyperten = (indirect_diff_`med'_hyperten/total_hyperten)*100
  }
     if (`out' == CVD_inc) {
	logistic CVD_inc eduyears_scaled $max_covar
	scalar total_CVD = exp(_b[eduyears_scaled])
        logistic CVD_inc `med' eduyears_scaled $max_covar 
        scalar med_CVD = exp(_b[`med'])
		scalar direct_CVD = exp(_b[eduyears_scaled])
		scalar indirect_diff_`med'_CVD_inc = total_CVD-direct_CVD
        scalar indirect_`med'_CVD_inc = exp_`med'*med_CVD
	scalar prop_`med'_CVD_inc = (indirect_`med'_CVD_inc/total_CVD)*100
	scalar diff_prop_`med'_CVD_inc = (indirect_diff_`med'_CVD_inc/total_CVD)*100
        }
  }
}
 end     

 foreach out in hyperten CVD_inc {
 local x=5
 local y=5
 local z=5
 local a=5
 foreach med in zldl { 
 
     putexcel set 202002_results_ldl, sheet(`out') modify
		reg `med'  eduyears_scaled $max_covar
         logistic `out' `med' eduyears_scaled $max_covar 
 bootstrap indirect = indirect_`med'_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]
 local x=`x'+1

		putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI' 
		
		logistic `out'  eduyears_scaled $max_covar
		logistic `out'  eduyears_scaled `med' $max_covar
		 bootstrap indirect_diff = indirect_diff_`med'_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect_diff = results[1,1]
         local indirect_LCI_diff = results[5,1]
         local indirect_UCI_diff = results[6,1]
         local a=`a'+1
     	 putexcel K`a'=`indirect_diff' L`a'=`indirect_LCI_diff' M`a'= `indirect_UCI_diff'
  
		 logistic `out' eduyears_scaled $max_covar
         reg `med'  eduyears_scaled $max_covar
         logistic `out' `med' eduyears_scaled $max_covar  	
bootstrap prop_med = prop_`med'_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local y = `y'+1 
 
      	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI'

	logistic `out'  eduyears_scaled $max_covar
	logistic `out'  eduyears_scaled `med' $max_covar
bootstrap diff_prop_`med'_`out', reps(100) nodrop : boot 
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local z = `z'+1
	 
   		putexcel N`z'=`proportion' O`z'=`proportion_LCI' P`z'=`proportion_UCI'


 }
}


*******************************************
* ldl MR mediation on the Risk diff scale *
*******************************************

capture program drop boot
program boot, rclass
	ivreg2 zldl (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in sbp hyperten CVD_inc {
		
	if (`out' == sbp)  {
	ivreg2 sbp (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_sbp = _b[eduyears_scaled]
	ivreg2 sbp (eduyears_scaled zldl = ea_weighted ldl_weighted) $PCs $max_covar
	scalar med_sbp = _b[zldl]
	scalar direct_sbp = _b[eduyears_scaled]
	scalar indirect_diff_sbp = total_sbp-direct_sbp
	scalar indirect_sbp = exp_med*med_sbp
	scalar prop_ldl_sbp = (indirect_sbp/total_sbp)*100 
	scalar diff_prop_ldl_sbp = (indirect_diff_sbp/total_sbp)*100 
	}
	if (`out' == hyperten) {
	ivreg2 hyperten (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_hyperten = _b[eduyears_scaled]
	ivreg2 hyperten (eduyears_scaled zldl = ea_weighted ldl_weighted) $PCs $max_covar
	scalar med_hyperten = _b[zldl]
	scalar direct_hyperten = _b[eduyears_scaled]
	scalar indirect_diff_hyperten = total_hyperten-direct_hyperten
	scalar indirect_hyperten = exp_med*med_hyperten
	scalar prop_ldl_hyperten = (indirect_hyperten/total_hyperten)*100 
	scalar diff_prop_ldl_hyperten = (indirect_diff_hyperten/total_hyperten)*100 
	}
	if (`out' == CVD_inc) {
	ivreg2 CVD_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_CVD_inc = _b[eduyears_scaled]
	ivreg2 CVD_inc (eduyears_scaled zldl = ea_weighted ldl_weighted) $PCs $max_covar
	scalar med_CVD_inc = _b[zldl]
	scalar direct_CVD_inc = _b[eduyears_scaled]
	scalar indirect_diff_CVD_inc = total_CVD_inc-direct_CVD_inc
	scalar indirect_CVD_inc = exp_med*med_CVD_inc
	scalar prop_ldl_CVD_inc = (indirect_CVD_inc/total_CVD_inc)*100 
	scalar diff_prop_ldl_CVD_inc = (indirect_diff_CVD_inc/total_CVD_inc)*100 
	}
}
end	
foreach out in sbp hyperten CVD_inc {
local x=2
local y=2
local z=2
local a=2
	putexcel set 202002_results_ldl, sheet(`out') modify
	ivreg2 zldl (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zldl = ea_weighted ldl_weighted) $PCs $max_covar

bootstrap indirect = indirect_`out', reps(100) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]
	local x=`x'+1
       putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	ivreg2 `out' (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zldl = ea_weighted ldl_weighted) $PCs $max_covar
bootstrap indirect_diff = indirect_diff_`out', reps(100) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
     
	local a=`a'+1 
	 
	putexcel K`a'=`proportion' L`a'=`proportion_LCI' M`a'=`proportion_UCI'	
			
	ivreg2 `out' (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 zldl (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zldl = ea_weighted ldl_weighted) $PCs $max_covar

bootstrap prop_med = prop_ldl_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]

	local y=`y'+1
		 
	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI'
	
	ivreg2 `out' (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zldl = ea_weighted ldl_weighted) $PCs $max_covar
bootstrap diff_prop_med = diff_prop_ldl_`out', reps(100) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
     
	local z=`z'+1 
	 
	putexcel N`z'=`proportion' O`z'=`proportion_LCI' P`z'=`proportion_UCI'
}


****************************************
* ldl MR mediation on the log OR scale *
****************************************

capture drop dir_ealdl
capture drop dir_ldlea
capture drop education
regress eduyears_scaled ea_weighted  $max_covar $PCs
predict education, xb
regress eduyears_scaled ea_weighted ldl_weighted $max_covar $PCs
predict dir_ealdl, xb
regress zldl ldl_weighted ea_weighted $max_covar $PCs
predict dir_ldlea, xb


capture program drop boot
program boot, rclass
	ivreg2 zldl (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in hyperten CVD_inc {

	if (`out' == hyperten) {
	logit hyperten education $max_covar $PCs 
	scalar total_hyperten = _b[education]
	logit hyperten dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 
	scalar med_hyperten = _b[dir_ldlea]
	scalar direct_hyperten = _b[dir_ealdl]
	scalar indirect_diff_hyperten = total_hyperten-direct_hyperten
	scalar indirect_hyperten = exp_med*med_hyperten
	scalar prop_ldl_hyperten = (indirect_hyperten/total_hyperten)*100 
	scalar diff_prop_ldl_hyperten = (indirect_diff_hyperten/total_hyperten)*100
}
	if (`out' == CVD_inc) {
	logit CVD_inc education $max_covar $PCs 
	scalar total_CVD_inc = _b[education]
	logit CVD_inc dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 
	scalar med_CVD_inc = _b[dir_ldlea]
	scalar direct_CVD_inc = _b[dir_ealdl]
	scalar indirect_diff_CVD_inc = total_CVD_inc-direct_CVD_inc
	scalar indirect_CVD_inc = exp_med*med_CVD_inc
	scalar prop_ldl_CVD_inc = (indirect_CVD_inc/total_CVD_inc)*100 
	scalar diff_prop_ldl_CVD_inc = (indirect_diff_CVD_inc/total_CVD_inc)*100
}
}
end	
foreach out in hyperten CVD_inc {
local x=4
local y=4
local z=4
local a=4
	putexcel set 202002_results_ldl, sheet(`out') modify
	ivreg2 zldl (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	logit `out' dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 
bootstrap indirect = indirect_`out', reps(100) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	local x=`x'+1
	
	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	logit `out' education $PCs $max_covar
	logit `out' dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 
bootstrap indirect_diff = indirect_diff_`out', reps(100) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
     
	local a=`a'+1 
	 
	putexcel K`a'=`proportion' L`a'=`proportion_LCI' M`a'=`proportion_UCI'	
	
	logit `out' education $max_covar $PCs 
	ivreg2 zldl (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	logit `out' dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 

bootstrap prop_med = prop_ldl_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
    
	local  y = `y'+1
	  
	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI'

	
	logit `out' education $max_covar $PCs 
	logit `out' dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 

bootstrap diff_prop_med = diff_prop_ldl_`out', reps(100) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	local z = `z'+1	
		
	putexcel N`z'=`proportion' O`z'=`proportion_LCI' P`z'=`proportion_UCI'
}

****************************************
* ldl MR mediation on the OR scale *
****************************************

capture drop dir_ealdl
capture drop dir_ldlea
capture drop education
regress eduyears_scaled ea_weighted  $max_covar $PCs
predict education, xb
regress eduyears_scaled ea_weighted ldl_weighted $max_covar $PCs
predict dir_ealdl, xb
regress zldl ldl_weighted ea_weighted $max_covar $PCs
predict dir_ldlea, xb


capture program drop boot
program boot, rclass
	ivreg2 zldl (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in hyperten CVD_inc {

	if (`out' == hyperten) {
	logistic hyperten education $max_covar $PCs 
	scalar total_hyperten = exp(_b[education])
	logistic hyperten dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 
	scalar med_hyperten = exp(_b[dir_ldlea])
	scalar direct_hyperten = exp(_b[dir_ealdl])
	scalar indirect_diff_hyperten = total_hyperten-direct_hyperten
	scalar indirect_hyperten = exp_med*med_hyperten
	scalar prop_ldl_hyperten = (indirect_hyperten/total_hyperten)*100 
	scalar diff_prop_ldl_hyperten = (indirect_diff_hyperten/total_hyperten)*100
}
	if (`out' == CVD_inc) {
	logistic CVD_inc education $max_covar $PCs 
	scalar total_CVD_inc = exp(_b[education])
	logistic CVD_inc dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 
	scalar med_CVD_inc = exp(_b[dir_ldlea])
	scalar direct_CVD_inc = exp(_b[dir_ealdl])
	scalar indirect_diff_CVD_inc = total_CVD_inc-direct_CVD_inc
	scalar indirect_CVD_inc = exp_med*med_CVD_inc
	scalar prop_ldl_CVD_inc = (indirect_CVD_inc/total_CVD_inc)*100 
	scalar diff_prop_ldl_CVD_inc = (indirect_diff_CVD_inc/total_CVD_inc)*100
}
}
end	
foreach out in hyperten CVD_inc {
local x=6
local y=6
local z=6
local a=6
	putexcel set 202002_results_ldl, sheet(`out') modify
	ivreg2 zldl (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	logistic `out' dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 
bootstrap indirect = indirect_`out', reps(100) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	local x=`x'+1
	
	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	logistic `out' education $PCs $max_covar
	logistic `out' dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 
bootstrap indirect_diff = indirect_diff_`out', reps(100) nodrop : boot
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]
	
	putexcel K`x'=`indirect' L`x'=`indirect_LCI' M`x'= `indirect_UCI'
	
     
	local a=`a'+1 
	
	logistic `out' education $max_covar $PCs 
	ivreg2 zldl (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	logistic `out' dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 

bootstrap prop_med = prop_ldl_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
    
	local  y = `y'+1
	  
	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI'

	
	logistic `out' education $max_covar $PCs 
	logistic `out' dir_ldlea dir_ealdl  $max_covar $PCs, vce(robust) 

bootstrap diff_prop_med = diff_prop_ldl_`out', reps(100) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	local z = `z'+1	
		
	putexcel N`z'=`proportion' O`z'=`proportion_LCI' P`z'=`proportion_UCI'
}
