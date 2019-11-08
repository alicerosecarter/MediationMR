/* Note on all bootstrapping: Bootstrapping should be used to estimate confidence
intervals. In this code we use 5 replications for efficiency, but this will 
usually be much more (>1000 reps) */

********************************************************************************
*							Setting the data								   *
********************************************************************************

global exposure edu
global mediator bmi
global outcome_cont sbp
global outcome_binary hyp
global covariates /* NOTE: Covariates should include all X-M AND all M-Y confounders */
global exposure_iv pgrs_edu
global mediator_iv pgrs_bmi
global PCs 

set more off
capture program drop boot

********************************************************************************
*							Observational methods							   *
********************************************************************************

********************************************************************************
							*Risk difference scale*

/* Note: This code only shows a continous outcome, but binary could also be used */							

************************		Total Effects		****************************

* Standard multivariable regression model
regress $outcome_cont $exposure $covariates


************************		Direct Effects		****************************

* Simply include the mediator(s) in the multivariable model with ALL covariates
regress $outcome_cont $exposure $covariates $mediator


**************		Indirect Effects & Proportion mediated		****************

* Specify the bootstrap model
********************************************************************************
program boot, rclass
/* 1 - Regress the exposure on the outcome(total effect) and store the scalar */
regress $outcome_cont $exposure $covariates
scalar total = _b[$exposure]
	
/* 2 - Regress the exposure on the mediator and store the scalar */
regress $mediator $exposure $covariates
scalar exposure_mediator = _b[$exposure]
       
/* 3 - Regress the mediator on the outcome ALSO adjusting for the mediator and store 
	the scalar of 1) the mediator and 2) the exposure
	Note: This model also gives you an estimate of the direct effect of the exposure (as above)
*/	
regress $outcome_cont $mediator $exposure $covariates 
scalar mediator_outcome = _b[$mediator]
scalar direct = _b[$exposure]
		
/*4 -  Using the stored scalars estimate the indirect effect using the difference method and store */	
scalar indirect_difference = total-direct

/* 5 - Using the stored scalars estimate the indirect effect using the product of coefficents method and store */	
scalar indirect_product = exposure_mediator*mediator_outcome
	
/*4 -  Using the stored scalars estimate the proportion mediated using the difference method and store */	
scalar prop_med_difference = (indirect_difference/total)*100
	
/* 5 - Using the stored scalars estimate the proportion mediated using the product of coefficents method and store */	
scalar prop_med_prod = (indirect_product/total)*100
end 
	
* Run the bootstrap models
********************************************************************************

/* 1 - Indirect effect using the difference method */
regress $outcome_cont $exposure $covariates
regress $outcome_cont $exposure $mediator $covariates
bootstrap indirect_difference = indirect_difference, reps(5) nodrop : boot

/* 2 - Indirect effect using the product of coefficients method */
regress $outcome_cont $exposure $covariates
regress $mediator $exposure $covariates
regress $outcome_cont $mediator $exposure $covariates 
bootstrap indirect_product = indirect_product, reps(5) nodrop : boot 
  
/* 3 - Proportion mediated using the difference method */ 
regress $outcome_cont $exposure $covariates
regress $outcome_cont $exposure $mediator $covariates
bootstrap prop_med_diff = prop_med_difference, reps(5) nodrop : boot 

/* 4 - Proportion mediated using the product of coefficients method */
regress $outcome_cont $exposure $covariates
regress $mediator $exposure $covariates
regress $outcome_cont $mediator $exposure $covariates 
bootstrap prop_med_prod = prop_med_prod, reps(5) nodrop : boot 		

capture program drop boot

********************************************************************************
						*Risk difference scale using sureg*
* More efficient code, where both models for the product method are included *

sureg ($mediator $exposure $max_covar) /// exposure-mediator model
	  ($outcome_cont $mediator $exposure $max_covar) /* mediator-outcome model/direct effect*/
	  
*Bootstrapped confidence intervals *

capture program drop bootmm

program boot, rclass
syntax [if] [in]
sureg (bmi edu) (sbp bmi edu ) `if' `in'
return scalar indirect = [bmi]_b[edu]*[sbp]_b[bmi]
end

bootstrap r(indirect), bca reps(5): boot

********************************************************************************
							*Log OR scale*
												

************************		Total Effects		****************************

* Standard multivariable regression model
logit $outcome_binary $exposure $covariates


************************		Direct Effects		****************************

* Simply include the mediator(s) in the multivariable model with ALL covariates
logit $outcome_binary $exposure $covariates $mediator


**************		Indirect Effects & Proportion mediated		****************

* Specify the bootstrap model
********************************************************************************
program boot, rclass
/* 1 - Regress the exposure on the outcome(total effect) and store the scalar */
logit $outcome_binary $exposure $covariates
scalar total = _b[$exposure]
	
/* 2 - Regress the exposure on the mediator and store the scalar */
regress $mediator $exposure $covariates
scalar exposure_mediator = _b[$exposure]
       
/* 3 - Regress the mediator on the outcome ALSO adjusting for the mediator and store 
	the scalar of 1) the mediator and 2) the exposure
	Note: This model also gives you an estimate of the direct effect of the exposure (as above)
*/	
logit $outcome_binary $mediator $exposure $covariates 
scalar mediator_outcome = _b[$mediator]
scalar direct = _b[$exposure]
		
/*4 -  Using the stored scalars estimate the indirect effect using the difference method and store */	
scalar indirect_difference = total-direct

/* 5 - Using the stored scalars estimate the indirect effect using the product of coefficents method and store */	
scalar indirect_product = exposure_mediator*mediator_outcome
	
/*4 -  Using the stored scalars estimate the proportion mediated using the difference method and store */	
scalar prop_med_difference = (indirect_difference/total)*100
	
/* 5 - Using the stored scalars estimate the proportion mediated using the product of coefficents method and store */	
scalar prop_med_prod = (indirect_product/total)*100
end
	
* Run the bootstrap models
********************************************************************************

/* 1 - Indirect effect using the difference method */
logit $outcome_binary $exposure $covariates
logit $outcome_binary $exposure $mediator $covariates
bootstrap indirect_difference = indirect_difference, reps(5) nodrop : boot

/* 2 - Indirect effect using the product of coefficients method */
logit $outcome_binary $exposure $covariates
regress $mediator $exposure $covariates
logit $outcome_binary $mediator $exposure $covariates 
bootstrap indirect_product = indirect_product, reps(5) nodrop : boot 
  
/* 3 - Proportion mediated using the difference method */ 
logit $outcome_binary $exposure $covariates
logit $outcome_binary $exposure $mediator $covariates
bootstrap prop_med_diff = prop_med_difference, reps(5) nodrop : boot 

/* 4 - Proportion mediated using the product of coefficients method */
logit $outcome_binary $exposure $covariates
regress $mediator $exposure $covariates
logit $outcome_binary $mediator $exposure $covariates 
bootstrap prop_med_prod = prop_med_prod, reps(5) nodrop : boot 		



capture program drop boot

********************************************************************************
*									MR methods								   *
********************************************************************************

/* NOTE: Not all MR analyses will need adjudting for covariates, especially given 
	the assumptions of MR. However, when using an exposure such as education 
	which we know is at least associated with population structure we should
	consider controlling for covariates including principal components */

********************************************************************************
							*Risk difference scale*

/* Note: This code only shows a continous outcome, but binary could also be used */									
						
************************		Total Effects		****************************

*Standard univariate MR
ivreg2 $outcome_cont ($exposure = $exposure_iv) $PCs $max_covar


************************		Direct Effects		****************************

*Multivariable MR model
ivreg2 $outcome_cont ($exposure $mediator = $exposure_iv $mediator_iv) $PCs $covariates


**************		Indirect Effects & Proportion mediated		****************

* Specify the bootstrap model
********************************************************************************
program boot, rclass
/* 1 - Regress the exposure on the outcome(total effect) and store the scalar */
ivreg2 $outcome_cont ($exposure = $exposure_iv) $PCs $max_covar
scalar total = _b[$exposure]
	
/* 2 - Regress the exposure on the mediator and store the scalar */
ivreg2 $mediator ($exposure = $exposure_iv) $PCs $max_covar
scalar exposure_mediator = _b[$exposure]
       
/* 3 - Regress the mediator on the outcome ALSO adjusting for the mediator and store 
	the scalar of 1) the mediator and 2) the exposure
	Note: This model also gives you an estimate of the direct effect of the exposure (as above)
*/	
ivreg2 $outcome_cont ($mediator $exposure = $mediator_iv $exposure_iv) $PCs $covariates
scalar mediator_outcome = _b[$mediator]
scalar direct = _b[$exposure]
		
/*4 -  Using the stored scalars estimate the indirect effect using the difference method and store */	
scalar indirect_difference = total-direct

/* 5 - Using the stored scalars estimate the indirect effect using the product of coefficents method and store */	
scalar indirect_product = exposure_mediator*mediator_outcome
	
/*4 -  Using the stored scalars estimate the proportion mediated using the difference method and store */	
scalar prop_med_difference = (indirect_difference/total)*100
	
/* 5 - Using the stored scalars estimate the proportion mediated using the product of coefficents method and store */	
scalar prop_med_prod = (indirect_product/total)*100
end
	
* Run the bootstrap models
********************************************************************************

/* 1 - Indirect effect using the difference method */
ivreg2 $outcome_cont ($exposure = $exposure_iv) $PCs $max_covar
ivreg2 $outcome_cont ($mediator $exposure = $mediator_iv $exposure_iv) $PCs $covariates
bootstrap indirect_difference = indirect_difference, reps(5) nodrop : boot

/* 2 - Indirect effect using the product of coefficients method */
ivreg2 $outcome_cont ($exposure = $exposure_iv) $PCs $max_covar
ivreg2 $mediator ($exposure = $exposure_iv) $PCs $max_covar
ivreg2 $outcome_cont ($mediator $exposure = $mediator_iv $exposure_iv) $PCs $covariates
bootstrap indirect_product = indirect_product, reps(5) nodrop : boot 
  
/* 3 - Proportion mediated using the difference method */ 
ivreg2 $outcome_cont ($exposure = $exposure_iv) $PCs $max_covar
ivreg2 $outcome_cont ($mediator $exposure = $mediator_iv $exposure_iv) $PCs $covariates
bootstrap prop_med_diff = prop_med_difference, reps(5) nodrop : boot 

/* 4 - Proportion mediated using the product of coefficients method */
ivreg2 $outcome_cont ($exposure = $exposure_iv) $PCs $max_covar
ivreg2 $mediator ($exposure = $exposure_iv) $PCs $max_covar
ivreg2 $outcome_cont ($mediator $exposure = $mediator_iv $exposure_iv) $PCs $covariates
bootstrap prop_med_prod = prop_med_prod, reps(5) nodrop : boot 		



capture program drop boot

********************************************************************************
							*Log OR scale*

************************		Total Effects		****************************

/* Estimate the gene-exposure association and store the values */
capture drop exposure
regress $exposure $exposure_iv $covariates $PCs
predict exposure, xb	

/* Use the predicted values of G-X in a regression for the outcome with robust standard errors */
logit $outcome_binary exposure $PCs $covariates, vce(robust)

************************		Direct Effects		****************************

/* Estimate the gene-exposure associations for both the exposure and mediator
	controlling for each other and store the values for both regressions */
capture drop gene_exposure
capture drop gene_mediator

* Exposure
regress $exposure $exposure_iv $mediator_iv $covariates $PCs
predict gene_exposure, xb

* Mediator
regress $mediator $mediator_iv $exposure_iv $covariates $PCs
predict gene_mediator, xb			

* Use the predicted values in the regression model for the outcome
logit $outcome_binary gene_exposure gene_mediator $covariates $PCs, vce(robust)	



**************		Indirect Effects & Proportion mediated		****************

* Estimate the first stage regression models (gene-exposure models)

capture drop gene_exposure
capture drop gene_mediator
capture drop exposure

regress $exposure $exposure_iv $covariates $PCs
predict exposure, xb

regress $exposure $exposure_iv $mediator_iv $covariates $PCs
predict gene_exposure, xb

regress $mediator $mediator_iv $exposure_iv $covariates $PCs
predict gene_mediator, xb


* Specify the bootstrap model
********************************************************************************
program boot, rclass
/* 1 - Regress the exposure on the outcome(total effect) and store the scalar */
logit $outcome_binary exposure $covariates $PCs 
scalar total = _b[exposure]
	
/* 2 - Regress the exposure on the mediator and store the scalar */
ivreg2 $mediator ($exposure = $exposure_iv) $PCs $max_covar
scalar exposure_mediator = _b[$exposure]
       
/* 3 - Regress the mediator on the outcome ALSO adjusting for the mediator and store 
	the scalar of 1) the mediator and 2) the exposure
	Note: This model also gives you an estimate of the direct effect of the exposure (as above)
*/	
logit $outcome_binary gene_mediator gene_exposure $covariates $PCs, vce(robust) 
scalar mediator_outcome = _b[gene_mediator]
scalar direct = _b[gene_exposure]
		
/*4 -  Using the stored scalars estimate the indirect effect using the difference method and store */	
scalar indirect_difference = total-direct

/* 5 - Using the stored scalars estimate the indirect effect using the product of coefficents method and store */	
scalar indirect_product = exposure_mediator*mediator_outcome
	
/*4 -  Using the stored scalars estimate the proportion mediated using the difference method and store */	
scalar prop_med_difference = (indirect_difference/total)*100
	
/* 5 - Using the stored scalars estimate the proportion mediated using the product of coefficents method and store */	
scalar prop_med_prod = (indirect_product/total)*100
end
	
* Run the bootstrap models
********************************************************************************

/* 1 - Indirect effect using the difference method */
logit $outcome_binary exposure $covariates $PCs 
logit $outcome_binary gene_mediator gene_exposure $covariates $PCs, vce(robust) 
bootstrap indirect_difference = indirect_difference, reps(5) nodrop : boot

/* 2 - Indirect effect using the product of coefficients method */
logit $outcome_binary exposure $covariates $PCs 
ivreg2 $mediator ($exposure = $exposure_iv) $PCs $max_covar
logit $outcome_binary gene_mediator gene_exposure $covariates $PCs, vce(robust) 
bootstrap indirect_product = indirect_product, reps(5) nodrop : boot 
  
/* 3 - Proportion mediated using the difference method */ 
logit $outcome_binary exposure $covariates $PCs 
logit $outcome_binary gene_mediator gene_exposure $covariates $PCs, vce(robust) 
bootstrap prop_med_diff = prop_med_difference, reps(5) nodrop : boot 

/* 4 - Proportion mediated using the product of coefficients method */
logit $outcome_binary exposure $covariates $PCs 
ivreg2 $mediator ($exposure = $exposure_iv) $PCs $max_covar
logit $outcome_binary gene_mediator gene_exposure $covariates $PCs, vce(robust) 
bootstrap prop_med_prod = prop_med_prod, reps(5) nodrop : boot 		


capture program drop boot
