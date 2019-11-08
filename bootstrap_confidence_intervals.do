global PCs PC01 PC02 PC03 PC04 PC05 PC06 PC07 PC08 PC09 PC10
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth
global med zbmi 
set more off


**************************************************
* Observational mediators on the Risk diff scale *
**************************************************

capture program drop boot
program boot, rclass
foreach med in zbmi { 
        
     if (`med' == zbmi) {
        regress zbmi eduyears_scaled $max_covar
        scalar exp_zbmi = _b[eduyears_scaled]
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
 foreach med in zbmi{ 
     putexcel set 20190516_results, sheet(`out') modify
       reg `med'  eduyears_scaled $max_covar
         reg `out' `med' eduyears_scaled $max_covar 
 bootstrap indirect = indirect_`med'_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]
         local x=`x'+1
		 putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'
		 
		reg `out'  eduyears_scaled $max_covar
		reg `out'  eduyears_scaled `med' $max_covar
		 bootstrap indirect_diff = indirect_diff_`med'_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect_diff = results[1,1]
         local indirect_LCI_diff = results[5,1]
         local indirect_UCI_diff = results[6,1]
         local a=`a'+1
     	 putexcel K`a'=`indirect_diff' L`a'=`indirect_LCI_diff' M`a'= `indirect_UCI_diff'

  	 
	 reg `out' eduyears_scaled $max_covar
         reg `med'  eduyears_scaled $max_covar
         reg `out' `med' eduyears_scaled $max_covar 

bootstrap prop_med = prop_`med'_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
         local y=`y'+1

      	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI' 
		
	reg `out'  eduyears_scaled $max_covar
	reg `out'  eduyears_scaled `med' $max_covar
bootstrap prop_med_diff = diff_prop_`med'_`out', reps(1000) nodrop : boot 
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
foreach med in zbmi { 
        
     if (`med' == zbmi) {
        regress zbmi eduyears_scaled $max_covar
        scalar exp_zbmi = _b[eduyears_scaled]
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
 foreach med in zbmi { 
 
     putexcel set 20190516_results, sheet(`out') modify
		reg `med'  eduyears_scaled $max_covar
         logit `out' `med' eduyears_scaled $max_covar 
 bootstrap indirect = indirect_`med'_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]
 local x=`x'+1

		putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI' 
  
  		logit `out'  eduyears_scaled $max_covar
		logit `out'  eduyears_scaled `med' $max_covar
		 bootstrap indirect_diff = indirect_diff_`med'_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect_diff = results[1,1]
         local indirect_LCI_diff = results[5,1]
         local indirect_UCI_diff = results[6,1]
         local a=`a'+1
     	 putexcel K`a'=`indirect_diff' L`a'=`indirect_LCI_diff' M`a'= `indirect_UCI_diff'
  
		 logit `out' eduyears_scaled $max_covar
         reg `med'  eduyears_scaled $max_covar
         logit `out' `med' eduyears_scaled $max_covar  	
bootstrap prop_med = prop_`med'_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local y = `y'+1 
 
      	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI'

	logit `out'  eduyears_scaled $max_covar
	logit `out'  eduyears_scaled `med' $max_covar
bootstrap diff_prop_`med'_`out', reps(1000) nodrop : boot 
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
foreach med in zbmi { 
        
     if (`med' == zbmi) {
        regress zbmi eduyears_scaled $max_covar
        scalar exp_zbmi = _b[eduyears_scaled]
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
 foreach med in zbmi { 
 
     putexcel set 20190516_results, sheet(`out') modify
		reg `med'  eduyears_scaled $max_covar
         logistic `out' `med' eduyears_scaled $max_covar 
 bootstrap indirect = indirect_`med'_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]
 local x=`x'+1

		putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI' 
		
		logistic `out'  eduyears_scaled $max_covar
		logistic `out'  eduyears_scaled `med' $max_covar
		 bootstrap indirect_diff = indirect_diff_`med'_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect_diff = results[1,1]
         local indirect_LCI_diff = results[5,1]
         local indirect_UCI_diff = results[6,1]
         local a=`a'+1
     	 putexcel K`a'=`indirect_diff' L`a'=`indirect_LCI_diff' M`a'= `indirect_UCI_diff'
  
		 logistic `out' eduyears_scaled $max_covar
         reg `med'  eduyears_scaled $max_covar
         logistic `out' `med' eduyears_scaled $max_covar  	
bootstrap prop_med = prop_`med'_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local y = `y'+1 
 
      	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI'

	logistic `out'  eduyears_scaled $max_covar
	logistic `out'  eduyears_scaled `med' $max_covar
bootstrap diff_prop_`med'_`out', reps(1000) nodrop : boot 
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local z = `z'+1
	 
   		putexcel N`z'=`proportion' O`z'=`proportion_LCI' P`z'=`proportion_UCI'


 }
}


*******************************************
* BMI MR mediation on the Risk diff scale *
*******************************************

capture program drop boot
program boot, rclass
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in sbp hyperten CVD_inc {
		
	if (`out' == sbp)  {
	ivreg2 sbp (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_sbp = _b[eduyears_scaled]
	ivreg2 sbp (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
	scalar med_sbp = _b[zbmi]
	scalar direct_sbp = _b[eduyears_scaled]
	scalar indirect_diff_sbp = total_sbp-direct_sbp
	scalar indirect_sbp = exp_med*med_sbp
	scalar prop_bmi_sbp = (indirect_sbp/total_sbp)*100 
	scalar diff_prop_bmi_`out' = (indirect_diff_sbp/total_sbp)*100 
	}
	if (`out' == hyperten) {
	ivreg2 hyperten (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_hyperten = _b[eduyears_scaled]
	ivreg2 hyperten (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
	scalar med_hyperten = _b[zbmi]
	scalar direct_hyperten = _b[eduyears_scaled]
	scalar indirect_diff_hyperten = total_hyperten-direct_hyperten
	scalar indirect_hyperten = exp_med*med_hyperten
	scalar prop_bmi_hyperten = (indirect_hyperten/total_hyperten)*100 
	scalar diff_prop_bmi_`out' = (indirect_diff_hyperten/total_hyperten)*100 
	}
	if (`out' == CVD_inc) {
	ivreg2 CVD_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_CVD_inc = _b[eduyears_scaled]
	ivreg2 CVD_inc (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
	scalar med_CVD_inc = _b[zbmi]
	scalar direct_CVD_inc = _b[eduyears_scaled]
	scalar indirect_diff_CVD_inc = total_CVD_inc-direct_CVD_inc
	scalar indirect_CVD_inc = exp_med*med_CVD_inc
	scalar prop_bmi_CVD_inc = (indirect_CVD_inc/total_CVD_inc)*100 
	scalar diff_prop_bmi_`out' = (indirect_diff_CVD_inc/total_CVD_inc)*100 
	}
}
end	
foreach out in sbp hyperten CVD_inc {
local x=2
local y=2
local z=2
local a=2
	putexcel set 20190516_results, sheet(`out') modify
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar

bootstrap indirect = indirect_`out', reps(1000) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]
	local x=`x'+1
       putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	ivreg2 `out' (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
bootstrap indirect_diff = indirect_diff_`out', reps(1000) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
     
	local a=`a'+1 
	 
	putexcel K`a'=`proportion' L`a'=`proportion_LCI' M`a'=`proportion_UCI'	
			
	ivreg2 `out' (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar

bootstrap prop_med = prop_bmi_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]

	local y=`y'+1
		 
	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI'
	
	ivreg2 `out' (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
bootstrap diff_prop_med = diff_prop_bmi_`out', reps(1000) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
     
	local z=`z'+1 
	 
	putexcel N`z'=`proportion' O`z'=`proportion_LCI' P`z'=`proportion_UCI'
}


****************************************
* BMI MR mediation on the log OR scale *
****************************************

capture drop dir_eabmi
capture drop dir_bmiea
capture drop education
regress eduyears_scaled ea_weighted  $max_covar $PCs
predict education, xb
regress eduyears_scaled ea_weighted zbmi_weighted $max_covar $PCs
predict dir_eabmi, xb
regress zbmi zbmi_weighted ea_weighted $max_covar $PCs
predict dir_bmiea, xb


capture program drop boot
program boot, rclass
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in hyperten CVD_inc {

	if (`out' == hyperten) {
	logit hyperten education $max_covar $PCs 
	scalar total_hyperten = _b[education]
	logit hyperten dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
	scalar med_hyperten = _b[dir_bmiea]
	scalar direct_hyperten = _b[dir_eabmi]
	scalar indirect_diff_hyperten = total_hyperten-direct_hyperten
	scalar indirect_hyperten = exp_med*med_hyperten
	scalar prop_bmi_hyperten = (indirect_hyperten/total_hyperten)*100 
	scalar diff_prop_bmi_hyperten = (indirect_diff_hyperten/total_hyperten)*100
}
	if (`out' == CVD_inc) {
	logit CVD_inc education $max_covar $PCs 
	scalar total_CVD_inc = _b[education]
	logit CVD_inc dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
	scalar med_CVD_inc = _b[dir_bmiea]
	scalar direct_CVD_inc = _b[dir_eabmi]
	scalar indirect_diff_CVD_inc = total_CVD_inc-direct_CVD_inc
	scalar indirect_CVD_inc = exp_med*med_CVD_inc
	scalar prop_bmi_CVD_inc = (indirect_CVD_inc/total_CVD_inc)*100 
	scalar diff_prop_bmi_CVD_inc = (indirect_diff_CVD_inc/total_CVD_inc)*100
}
}
end	
foreach out in hyperten CVD_inc {
local x=4
local y=4
local z=4
local a=4
	putexcel set 20190516_results, sheet(`out') modify
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	logit `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
bootstrap indirect = indirect_`out', reps(1000) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	local x=`x'+1
	
	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	logit `out' education $PCs $max_covar
	logit `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
bootstrap indirect_diff = indirect_diff_`out', reps(1000) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
     
	local a=`a'+1 
	 
	putexcel K`a'=`proportion' L`a'=`proportion_LCI' M`a'=`proportion_UCI'	
	
	logit `out' education $max_covar $PCs 
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	logit `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 

bootstrap prop_med = prop_bmi_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
    
	local  y = `y'+1
	  
	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI'

	
	logit `out' education $max_covar $PCs 
	logit `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 

bootstrap diff_prop_med = diff_prop_bmi_`out', reps(1000) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	local z = `z'+1	
		
	putexcel N`z'=`proportion' O`z'=`proportion_LCI' P`z'=`proportion_UCI'
}

****************************************
* BMI MR mediation on the OR scale *
****************************************

capture drop dir_eabmi
capture drop dir_bmiea
capture drop education
regress eduyears_scaled ea_weighted  $max_covar $PCs
predict education, xb
regress eduyears_scaled ea_weighted zbmi_weighted $max_covar $PCs
predict dir_eabmi, xb
regress zbmi zbmi_weighted ea_weighted $max_covar $PCs
predict dir_bmiea, xb


capture program drop boot
program boot, rclass
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in hyperten CVD_inc {

	if (`out' == hyperten) {
	logistic hyperten education $max_covar $PCs 
	scalar total_hyperten = exp(_b[education])
	logistic hyperten dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
	scalar med_hyperten = exp(_b[dir_bmiea])
	scalar direct_hyperten = exp(_b[dir_eabmi])
	scalar indirect_diff_hyperten = total_hyperten-direct_hyperten
	scalar indirect_hyperten = exp_med*med_hyperten
	scalar prop_bmi_hyperten = (indirect_hyperten/total_hyperten)*100 
	scalar diff_prop_bmi_hyperten = (indirect_diff_hyperten/total_hyperten)*100
}
	if (`out' == CVD_inc) {
	logistic CVD_inc education $max_covar $PCs 
	scalar total_CVD_inc = exp(_b[education])
	logistic CVD_inc dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
	scalar med_CVD_inc = exp(_b[dir_bmiea])
	scalar direct_CVD_inc = exp(_b[dir_eabmi])
	scalar indirect_diff_CVD_inc = total_CVD_inc-direct_CVD_inc
	scalar indirect_CVD_inc = exp_med*med_CVD_inc
	scalar prop_bmi_CVD_inc = (indirect_CVD_inc/total_CVD_inc)*100 
	scalar diff_prop_bmi_CVD_inc = (indirect_diff_CVD_inc/total_CVD_inc)*100
}
}
end	
foreach out in hyperten CVD_inc {
local x=6
local y=6
local z=6
local a=6
	putexcel set 20190516_results, sheet(`out') modify
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	logistic `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
bootstrap indirect = indirect_`out', reps(1000) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	local x=`x'+1
	
	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	logistic `out' education $PCs $max_covar
	logistic `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
bootstrap indirect_diff = indirect_diff_`out', reps(1000) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
     
	local a=`a'+1 
	
	logistic `out' education $max_covar $PCs 
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	logistic `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 

bootstrap prop_med = prop_bmi_`out', reps(1000) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
    
	local  y = `y'+1
	  
	putexcel T`y'=`proportion' U`y'=`proportion_LCI' V`y'= `proportion_UCI'

	
	logistic `out' education $max_covar $PCs 
	logistic `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 

bootstrap diff_prop_med = diff_prop_bmi_`out', reps(1000) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	local z = `z'+1	
		
	putexcel N`z'=`proportion' O`z'=`proportion_LCI' P`z'=`proportion_UCI'
}
