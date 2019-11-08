#install.packages(foreign)
library(foreign)
sample_mediation <- read.dta("simulated_dataset_stata12.dta")

#install.packages("AER")
library("AER")

#install.packages("boot")
library ("boot")

## 1 - Code for two-step MR - Product of coefficients method
# Total effect of X (education) on continuous Y (SBP)
total_sbp <- ivreg(sbp ~ edu | pgrs_edu, data = sample_mediation)

# Effect of X on M (BMI)
edu_bmi <- ivreg(bmi ~ edu | pgrs_edu, data = sample_mediation)

# Effect of M on Y
bmi_sbp <- ivreg(sbp ~ bmi + edu | pgrs_bmi + pgrs_edu, data = sample_mediation)

# Indirect effect (including bootstrapped confidence intervals)
indirect_effect_product <- edu_bmi$coef[2]*bmi_sbp$coef[2]

set.seed(1234)
indirect_product <- function(data, indices) { 
  sample <- data[indices,]
  edu_bmi <- ivreg(bmi ~ edu | pgrs_edu, data = sample)
  bmi_sbp <- ivreg(sbp ~ bmi + edu | pgrs_bmi + pgrs_edu, data = sample)
  return(edu_bmi$coef[2]*bmi_sbp$coef[2])
  }

boot_indirect_product <- boot(sample_mediation, indirect_product,R=500)
boot.ci(boot.out= boot_indirect_product, type = c("norm"))

# Proportion mediated (including bootstrapped confidence intervals)
proportion_effect_product <- indirect_effect_product/total_sbp$coef[2]

set.seed(1234)
proportion_product <- function(data, indices) { 
  sample <- data[indices,]
  total_sbp <- ivreg(sbp ~ edu | pgrs_edu, data = sample)
  edu_bmi <- ivreg(bmi ~ edu | pgrs_edu, data = sample)
  bmi_sbp <- ivreg(sbp ~ bmi + edu | pgrs_bmi + pgrs_edu, data = sample)
  return((edu_bmi$coef[2]*bmi_sbp$coef[2])/total_sbp$coef[2])
}

boot_proportion_product <- boot(sample_mediation, proportion_product,R=500)
boot.ci(boot.out= boot_proportion_product, type = c("norm"))

## 2 - Code for two-step MR - Difference in coefficients method
# Total effect of X (education) on continuous Y (SBP)
total_sbp <- ivreg(sbp ~ edu | pgrs_edu, data = sample_mediation)

# Direct effect of X on Y controlling for M
direct <- ivreg(sbp ~ edu + bmi | pgrs_edu + pgrs_bmi, data = sample_mediation)

# Indirect effect (including bootstrapped confidence intervals)
indirect_effect_difference <- total_sbp$coef[2]-direct$coef[2]

set.seed(1234)
indirect_difference <- function(data, indices) { 
  sample <- data[indices,]
  total_sbp <- ivreg(sbp ~ edu | pgrs_edu, data = sample)
  direct <- ivreg(sbp ~ edu + bmi | pgrs_edu + pgrs_bmi, data = sample)
  return(total_sbp$coef[2]-direct$coef[2])
}

boot_indirect_difference <- boot(sample_mediation, indirect_difference,R=500)
boot.ci(boot.out= boot_indirect_product, type = c("norm"))

# Proportion mediated (including bootstrapped confidence intervals)
proportion_effect_difference <- indirect_effect_difference/total_sbp$coef[2]

set.seed(1234)
proportion_difference <- function(data, indices) { 
  sample <- data[indices,]
  total_sbp <- ivreg(sbp ~ edu | pgrs_edu, data = sample)
  direct <- ivreg(sbp ~ edu + bmi | pgrs_edu + pgrs_bmi, data = sample)
  return((total_sbp$coef[2]-direct$coef[2])/total_sbp$coef[2])
}

boot_proportion_difference <- boot(sample_mediation, proportion_difference, R=500)
boot.ci(boot.out= boot_proportion_difference, type = c("norm"))
