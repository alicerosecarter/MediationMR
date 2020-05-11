
global PCs PC*
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth


********************************************************************************

foreach out in sbp hyperten CVD_inc {

	putexcel set 202002_results_combined, sheet(`out') modify
	putexcel A1="Method" B1="Scale" C1="Mediator" D1="Outcome" E1="Total Effect" F1="Total - LCI" G1="Total - UCI" H1="Direct Effect" I1="Direct - LCI" J1="Direct - UCI" K1="Indirect - difference" ///
		L1="LCI" M1="UCI" N1="% mediated - difference" O1="LCI" P1="UCI"
}
******************** Total Effects - Observational on RD scale *****************

foreach out in sbp hyperten CVD_inc {
	
		local x=2	

	
	putexcel set 202002_results_combined, sheet(`out') modify

		regress `out' eduyears_scaled $max_covar
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="Observational" B`x'="Risk Difference" C`x'="BMI+LDL" D`x'="`out_label'" E`x'=_b[eduyears_scaled] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	


******************* Total Effects - Observational on logOR scale ***************

foreach out in hyperten CVD_inc {
	
		local x=3
	
	foreach med in zbmi  {


	putexcel set 202002_results_combined, sheet(`out') modify

		logit `out' eduyears_scaled $max_covar

	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="Observational" B`x'="Log OR" C`x'="BMI+LDL" D`x'="`out_label'" E`x'=_b[eduyears_scaled] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}

******************** Total Effects - One-sample MR on RD scale *****************

foreach out in sbp hyperten CVD_inc {
	
		local x=4
		
	foreach med in zbmi  {

		
	putexcel set 202002_results_combined, sheet(`out') modify

		ivreg2 `out' (eduyears_scaled = ea_weighted) $PCs $max_covar
	
	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="MR - One-sample" B`x'="Risk Difference" C`x'="BMI+LDL" D`x'="`out_label'" E`x'=_b[eduyears_scaled] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}


******************* Total Effects - One-sample MR on logOR scale ***************
capture drop education_MR
regress eduyears_scaled ea_weighted $max_covar $PCs
predict education_MR, xb

foreach out in  hyperten CVD_inc {
	
		local x=5
	
	foreach med in zbmi  {

	
	putexcel set 202002_results_combined, sheet(`out') modify

		logit `out' education_MR $PCs $max_covar
	
	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[education_MR]-1.96*_se[education_MR]
	local ul_2sls1 = _b[education_MR]+1.96*_se[education_MR]
	
	putexcel A`x'="MR - One-sample" B`x'="Log OR" C`x'="BMI+LDL" D`x'="`out_label'" E`x'=_b[education_MR] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}

******* Direct Effects for Multiple Mediators - Observational on RD scale ******

foreach out in  sbp hyperten CVD_inc {
	
		local x=2	
	
	putexcel set 202002_results_combined, sheet(`out') modify

	regress `out' eduyears_scaled zbmi zldl  $max_covar
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel H`x'=_b[eduyears_scaled] I`x'=`ll_2sls1' J`x'= `ul_2sls1'

	}	
	
******* Direct Effects for Multiple Mediators - Observational on log OR scale ******


foreach out in  hyperten CVD_inc {
	
	local x=3	
	
	putexcel set 202002_results_combined, sheet(`out') modify

	logit `out' eduyears_scaled zbmi zldl $max_covar
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel H`x'=_b[eduyears_scaled] I`x'=`ll_2sls1' J`x'= `ul_2sls1'

	}	

*************** Direct - One-sample MR on RD scale **********	
	
	foreach out in sbp hyperten CVD_inc {
	
		local x=4
	
	putexcel set 202002_results_combined, sheet(`out') modify

	ivreg2 `out' (eduyears_scaled zbmi zldl = ea_weighted bmi_weighted ldl_weighted) $PCs $max_covar
	
	
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel H`x'=_b[eduyears_scaled] I`x'=`ll_2sls1' J`x'= `ul_2sls1'

	}

	
*********** Direct & Indirect Effects - One-sample MR on log OR scale **********

							

capture drop dir_ea
capture drop dir_bmi
capture drop dir_ldl
regress eduyears_scaled ea_weighted bmi_weighted ldl_weighted $max_covar $PCs
predict dir_ea, xb
regress zbmi bmi_weighted ea_weighted ldl_weighted $max_covar $PCs
predict dir_bmi, xb
regress zldl ldl_weighted bmi_weighted ea_weighted $max_covar $PCs
predict dir_ldl, xb


foreach out in  hyperten  CVD_inc {

	putexcel set 202002_results_combined, sheet(`out') modify
	
	local x=5
	
	logit `out' dir_ea dir_bmi dir_ldl $max_covar $PCs, vce(robust) // Natural Direct effect of Education
	
	local out_label : var label `out'
	local ll_2sls1 = _b[dir_ea]-1.96*_se[dir_ea]
	local ul_2sls1 = _b[dir_ea]+1.96*_se[dir_ea]
		
	putexcel H`x'=_b[dir_ea] I`x'=`ll_2sls1' J`x'= `ul_2sls1'
	
	
}	


global PCs PC*
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth

rename eduyears_scaled eduyears_3_6
egen eduyears_scaled = std(eduyears_3_6)

set more off
global med zbmi zldl
capture program drop boot
program boot, rclass

foreach out in sbp hyperten CVD_inc {     

     if (`out' == sbp) {
	reg sbp eduyears_scaled $max_covar
	scalar total_sbp = _b[eduyears_scaled]
        regress sbp $med eduyears_scaled $max_covar 
		scalar direct_sbp = _b[eduyears_scaled]
		scalar indirect_diff_sbp = total_sbp-direct_sbp
	scalar diff_prop_sbp = (indirect_diff_sbp/total_sbp)*100
        }
     if (`out' == hyperten) {
	reg hyperten eduyears_scaled $max_covar
	scalar total_hyperten = _b[eduyears_scaled]
        regress hyperten $med eduyears_scaled $max_covar 
		scalar direct_hyperten = _b[eduyears_scaled]
		scalar indirect_diff_hyperten = total_hyperten-direct_hyperten
	scalar diff_prop_hyperten = (indirect_diff_hyperten/total_hyperten)*100
        }
		     if (`out' == CVD_inc) {
	reg CVD_inc eduyears_scaled $max_covar
	scalar total_CVD = _b[eduyears_scaled]
        regress CVD_inc $med eduyears_scaled $max_covar 
		scalar direct_CVD = _b[eduyears_scaled]
		scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar diff_prop_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
        }
}

 end     

foreach out in sbp hyperten CVD_inc {     
 local x=2
 local y=2

     	 putexcel set 202002_results_combined, sheet(`out') modify
         reg `out' eduyears_scaled $max_covar
         reg `out' $med eduyears_scaled $max_covar 
 
 bootstrap indirect = indirect_diff_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]

     	
		 putexcel K`x'=`indirect' L`x'=`indirect_LCI' M`x'= `indirect_UCI'
  	 
		
	reg `out'  eduyears_scaled $max_covar
	reg `out'  eduyears_scaled $med $max_covar
bootstrap diff_prop_`out', reps(100) nodrop : boot 
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 
   		putexcel N`y'=`proportion' O`y'=`proportion_LCI' P`y'= `proportion_UCI'
 }

*******************************************************
* Observational multiple mediators on the LogOR scale *
*******************************************************

global med zbmi zldl
capture program drop boot
program boot, rclass

foreach out in hyperten CVD_inc {     

     
     if (`out' == hyperten) {
	logit hyperten eduyears_scaled $max_covar
	scalar total_hyperten = _b[eduyears_scaled]
        logit hyperten $med eduyears_scaled $max_covar 
		scalar direct_hyperten = _b[eduyears_scaled]
		scalar indirect_diff_hyperten = total_hyperten-direct_hyperten
	scalar diff_prop_hyperten = (indirect_diff_hyperten/total_hyperten)*100
        }
		     if (`out' == CVD_inc) {
	logit CVD_inc eduyears_scaled $max_covar
	scalar total_CVD = _b[eduyears_scaled]
        logit CVD_inc $med eduyears_scaled $max_covar 
		scalar direct_CVD = _b[eduyears_scaled]
		scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar diff_prop_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
        }
}

 end     

foreach out in hyperten CVD_inc {     
 local x=3
 local y=3

     	 putexcel set 202002_results_combined, sheet(`out') modify
	 logit `out'  eduyears_scaled $max_covar
	logit `out'  eduyears_scaled $med $max_covar
 
 bootstrap indirect = indirect_diff_`out', reps(100) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]

     	
		 putexcel K`x'=`indirect' L`x'=`indirect_LCI' M`x'= `indirect_UCI'
  	 
		
	reg `out'  eduyears_scaled $max_covar
	reg `out'  eduyears_scaled $med $max_covar
bootstrap diff_prop_`out', reps(100) nodrop : boot 
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 
   		putexcel N`y'=`proportion' O`y'=`proportion_LCI' P`y'= `proportion_UCI'
 }

 
*******************************************
* MR mediation on the Risk diff scale *
*******************************************

capture program drop boot
program boot, rclass

foreach out in sbp hyperten CVD_inc {
		
	if (`out' == sbp)  {
	ivreg2 sbp (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_sbp = _b[eduyears_scaled]
	ivreg2 sbp (eduyears_scaled zbmi zldl = ea_weighted bmi_weighted ldl_weighted) $PCs $max_covar
	scalar direct_sbp = _b[eduyears_scaled]
	scalar indirect_diff_sbp = total_sbp-direct_sbp
	scalar diff_prop_bmi_`out' = (indirect_diff_sbp/total_sbp)*100 
	}
	if (`out' == hyperten) {
	ivreg2 hyperten (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_hyperten = _b[eduyears_scaled]
	ivreg2 hyperten (eduyears_scaled zbmi zldl = ea_weighted bmi_weighted ldl_weighted) $PCs $max_covar
	scalar direct_hyperten = _b[eduyears_scaled]
	scalar indirect_diff_hyperten = total_hyperten-direct_hyperten
	scalar diff_prop_bmi_`out' = (indirect_diff_hyperten/total_hyperten)*100 
	}
	if (`out' == CVD_inc) {
	ivreg2 CVD_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_CVD_inc = _b[eduyears_scaled]
	ivreg2 CVD_inc (eduyears_scaled zbmi zldl = ea_weighted bmi_weighted ldl_weighted) $PCs $max_covar
	scalar direct_CVD_inc = _b[eduyears_scaled]
	scalar indirect_diff_CVD_inc = total_CVD_inc-direct_CVD_inc
	scalar diff_prop_bmi_`out' = (indirect_diff_CVD_inc/total_CVD_inc)*100 
	}
}
end	
foreach out in sbp hyperten CVD_inc {
local x=4
local y=4

	putexcel set 202002_results_combined, sheet(`out') modify

	ivreg2 `out' (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zbmi zldl = ea_weighted bmi_weighted ldl_weighted) $PCs $max_covar
bootstrap indirect_diff = indirect_diff_`out', reps(100) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
     
	 
	putexcel K`x'=`proportion' L`x'=`proportion_LCI' M`x'=`proportion_UCI'	
			
	
	ivreg2 `out' (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zbmi zldl = ea_weighted bmi_weighted ldl_weighted) $PCs $max_covar
bootstrap diff_prop_med = diff_prop_bmi_`out', reps(100) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
     
	putexcel N`y'=`proportion' O`y'=`proportion_LCI' P`y'=`proportion_UCI'
}


*************************************
*  MR mediation on the log OR scale *
*************************************

capture drop dir_ea
capture drop dir_bmi
capture drop dir_ldl
regress eduyears_scaled ea_weighted $max_covar $PCs
predict education, xb
regress eduyears_scaled ea_weighted bmi_weighted ldl_weighted $max_covar $PCs
predict dir_ea, xb
regress zbmi bmi_weighted ea_weighted ldl_weighted $max_covar $PCs
predict dir_bmi, xb
regress zldl ldl_weighted bmi_weighted ea_weighted $max_covar $PCs
predict dir_ldl, xb


capture program drop boot
program boot, rclass
 
foreach out in hyperten CVD_inc {

	if (`out' == hyperten) {
	logit hyperten education $max_covar $PCs 
	scalar total_hyperten = _b[education]
	logit `out' dir_ea dir_bmi dir_ldl   $max_covar $PCs, vce(robust) 
	scalar direct_hyperten = _b[dir_ea]
	scalar indirect_diff_hyperten = total_hyperten-direct_hyperten
	scalar diff_prop_bmi_hyperten = (indirect_diff_hyperten/total_hyperten)*100
}
	if (`out' == CVD_inc) {
	logit CVD_inc education $max_covar $PCs 
	scalar total_CVD_inc = _b[education]
	logit `out' dir_ea dir_bmi dir_ldl   $max_covar $PCs, vce(robust) 
	scalar direct_CVD_inc = _b[dir_ea]
	scalar indirect_diff_CVD_inc = total_CVD_inc-direct_CVD_inc
	scalar diff_prop_bmi_CVD_inc = (indirect_diff_CVD_inc/total_CVD_inc)*100
}
}
end	
foreach out in hyperten CVD_inc {
local x=5
local y=5

	putexcel set 202002_results_combined, sheet(`out') modify

	logit `out' education $PCs $max_covar
	logit `out' dir_ea dir_bmi dir_ldl   $max_covar $PCs, vce(robust) 
bootstrap indirect_diff = indirect_diff_`out', reps(100) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 
	putexcel K`x'=`proportion' L`x'=`proportion_LCI' M`x'=`proportion_UCI'	

	
	logit `out' education $max_covar $PCs 
	logit `out' dir_ea dir_bmi dir_ldl   $max_covar $PCs, vce(robust) 

bootstrap diff_prop_med = diff_prop_bmi_`out', reps(100) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]

		
	putexcel N`y'=`proportion' O`y'=`proportion_LCI' P`y'=`proportion_UCI'
}
