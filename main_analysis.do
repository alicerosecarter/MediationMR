global PCs PC01 PC02 PC03 PC04 PC05 PC06 PC07 PC08 PC09 PC10
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth

********************************************************************************

foreach out in sbp hyperten CVD_inc {

	putexcel set 20190516_results, sheet(`out') modify
	putexcel A1="Method" B1="Scale" C1="Mediator" D1="Outcome" E1="Total Effect" F1="Total - LCI" G1="Total - UCI" H1="Direct Effect" I1="Direct - LCI" J1="Direct - UCI" K1="Indirect - difference" ///
		L1="LCI" M1="UCI" N1="% mediated - difference" O1="LCI" P1="UCI" Q1="Indirect - product" R1="LCI" S1="UCI" T1="% mediated - product" U1="LCI" V1="UCI" ///
		W1="Exposure - Mediator" X1="Exp-Med - LCI" Y1="Exp-Med - UCI" Z1="Mediator - Outcome" AA1="Med-Out - LCI" AB1="Med-Out - UCI" 
}

******************** Total Effects - Observational on RD scale *****************

foreach out in sbp hyperten CVD_inc {
	
		local x=1	
	
		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

		regress `out' eduyears_scaled $max_covar
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="Observational" B`x'="Risk Difference" C`x'="BMI" D`x'="`out_label'" E`x'=_b[eduyears_scaled] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	


******************* Total Effects - Observational on logOR scale ***************

foreach out in hyperten CVD_inc {
	
		local x=3
	
	foreach med in zbmi  {

		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

		logit `out' eduyears_scaled $max_covar

	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="Observational" B`x'="Log OR" C`x'="BMI" D`x'="`out_label'" E`x'=_b[eduyears_scaled] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}

******************* Total Effects - Observational on OR scale ***************

foreach out in hyperten CVD_inc {
	
		local x=5
	
	foreach med in zbmi  {

		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

		logit `out' eduyears_scaled $max_covar

	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local beta = exp(_b[eduyears_scaled])
	local ll_2sls1 = exp(_b[eduyears_scaled])-1.96*(_se[eduyears_scaled])
	local ul_2sls1 = exp(_b[eduyears_scaled])+1.96*(_se[eduyears_scaled])
	
	putexcel A`x'="Observational" B`x'="Log OR" C`x'="BMI" D`x'="`out_label'" E`x'=`beta' F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}

******************** Total Effects - One-sample MR on RD scale *****************

foreach out in sbp hyperten CVD_inc {
	
		local x=2
		
	foreach med in zbmi  {

		local x=`x'+1
		
	putexcel set 20190516_results, sheet(`out') modify

		ivreg2 `out' (eduyears_scaled = ea_weighted) $PCs $max_covar
	
	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="MR - One-sample" B`x'="Risk Difference" C`x'="BMI" D`x'="`out_label'" E`x'=_b[eduyears_scaled] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}


******************* Total Effects - One-sample MR on logOR scale ***************
capture drop education_MR
regress eduyears_scaled ea_weighted $max_covar $PCs
predict education_MR, xb

foreach out in  hyperten CVD_inc {
	
		local x=4
	
	foreach med in zbmi  {

		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

		logit `out' education_MR $PCs $max_covar
	
	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[education_MR]-1.96*_se[education_MR]
	local ul_2sls1 = _b[education_MR]+1.96*_se[education_MR]
	
	putexcel A`x'="MR - One-sample" B`x'="Log OR" C`x'="BMI" D`x'="`out_label'" E`x'=_b[education_MR] F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}

******************* Total Effects - One-sample MR on OR scale ***************
capture drop education_MR
regress eduyears_scaled ea_weighted $max_covar $PCs
predict education_MR, xb

foreach out in  hyperten CVD_inc {
	
		local x=6
	
	foreach med in zbmi  {

		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

		logistic `out' education_MR $PCs $max_covar
	
	
	matrix results = r(table)
	local out_label : var label `out'
	local med_label : var label `med'
	local beta = results[1,1]
	local ll_2sls1 = results[5,1]
	local ul_2sls1 = results[6,1]
	
	putexcel A`x'="MR - One-sample" B`x'="OR" C`x'="BMI" D`x'="`out_label'" E`x'=`beta' F`x'=`ll_2sls1' G`x'=`ul_2sls1'
	
	}	
}

************* Direct & Indirect Effects - Observational on RD scale ************

foreach out in sbp hyperten CVD_inc {
	
		local x=1	
	
	foreach med in zbmi {
		
		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

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
	
	foreach med in zbmi{

		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

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
	
	foreach med in zbmi{

		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

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
							
										*BMI*
foreach out in sbp hyperten CVD_inc {
	
		local x=2

		local x=`x'+1
		
	putexcel set 20190516_results, sheet(`out') modify

	ivreg2 `out' (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
	
	
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[zbmi]-1.96*_se[zbmi]
	local ul_med = _b[zbmi]+1.96*_se[zbmi]
	
	putexcel H`x'=_b[eduyears_scaled] I`x'=`ll_2sls1' J`x'= `ul_2sls1'
	putexcel Z`x'=_b[zbmi] AA`x'=`ll_med' AB`x'= `ul_med'

	}

*********** Direct & Indirect Effects - One-sample MR on log OR scale **********

							*BMI*

capture drop dir_eabmi
capture drop dir_bmiea
regress eduyears_scaled ea_weighted zbmi_weighted $max_covar $PCs
predict dir_eabmi, xb
regress zbmi zbmi_weighted ea_weighted $max_covar $PCs
predict dir_bmiea, xb


foreach out in  hyperten  CVD_inc {

	putexcel set 20190516_results, sheet(`out') modify
	
	local x=4
	
	logit `out' dir_eabmi dir_bmiea $max_covar $PCs, vce(robust) // Natural Direct effect of Education
	
	local out_label : var label `out'
	local ll_2sls1 = _b[dir_eabmi]-1.96*_se[dir_eabmi]
	local ul_2sls1 = _b[dir_eabmi]+1.96*_se[dir_eabmi]
	local ll_med = _b[dir_bmiea]-1.96*_se[dir_bmiea]
	local ul_med = _b[dir_bmiea]+1.96*_se[dir_bmiea]
	
	local x=`x'+1	
		
	putexcel H`x'=_b[dir_eabmi] I`x'=`ll_2sls1' J`x'= `ul_2sls1'
	putexcel Z`x'=_b[dir_bmiea] AA`x'=`ll_med' AB`x'= `ul_med'
	
	
}	


*********** Direct & Indirect Effects - One-sample MR on OR scale **********

							*BMI*

capture drop dir_eabmi
capture drop dir_bmiea
regress eduyears_scaled ea_weighted zbmi_weighted $max_covar $PCs
predict dir_eabmi, xb
regress zbmi zbmi_weighted ea_weighted $max_covar $PCs
predict dir_bmiea, xb


foreach out in  hyperten  CVD_inc {

	local x=6

	putexcel set 20190516_results, sheet(`out') modify
	
	logistic `out' dir_eabmi dir_bmiea $max_covar $PCs, vce(robust) // Natural Direct effect of Education
	
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
	
	foreach med in zbmi {
		
		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

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
	
	foreach med in zbmi {
		
		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

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
	
	foreach med in zbmi {
		
		local x=`x'+1
	
	putexcel set 20190516_results, sheet(`out') modify

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
		
	foreach med in zbmi {

		local x=`x'+1
		
	putexcel set 20190516_results, sheet(`out') modify

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
		
	foreach med in zbmi {

		local x=`x'+1
		
	putexcel set 20190516_results, sheet(`out') modify

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
		
	foreach med in zbmi {

		local x=`x'+1
		
	putexcel set 20190516_results, sheet(`out') modify

		ivreg2 `med' (eduyears_scaled = ea_weighted) $PCs $max_covar
	
		
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel W`x'=_b[eduyears_scaled] X`x'=`ll_2sls1' Y`x'=`ul_2sls1'
	
	}	
}	
