---
title: "Effect of misclassification on effect size"
output: 
  html_document:
    toc: true
    keep_md: true
date: '2022-12-10'
---



## What we do

We look into how misclassifications may affect the estimated effect size of a genetic variant on OSA. Misclassification is assumed only in one direction: an individuals with OSA may not have any indication as having OSA, and is therefore misclassified as a control. For now, we ignore potential age/gender effects. 




```r
library(data.table)
expit <- function(x){exp(x)/(1 + exp(x))}
n <- 100000 # large n for accurate prevalence

g <- rbinom(n, size = 2, prob = 0.3)

# compute OSA prevalence for a few potential simulation settings: 
prev <- data.frame(intercept = c(-2, -1.5, -1, -0.5), prevalence = NA)

for (i in 1:nrow(prev)){
  OSA_prob <- expit(prev$intercept[i] + 0.1*g)
  OSA <- rbinom(n = length(OSA_prob), size = 1, prob = OSA_prob)
  prev[i, "prevalence"] <- round(mean(OSA)*100, 1)
}

prev
```

```
##   intercept prevalence
## 1      -2.0       12.5
## 2      -1.5       19.2
## 3      -1.0       28.2
## 4      -0.5       39.2
```


Now let's perform a simulation study. We will vary the prevalence of  OSA, as well as the sensitivity (rate of misclassification of OSA cases as controls). 


```r
n <- 20000 # 
# choose intercept values that lead to realistic prevalence
OSA_prob_intercept <- c(-1.5, -1, -0.5)
# choose realistic miscalssifcation rates: 
OSA_misclassification_rate <- c(0.4, 0.6, 0.8)
n_rep <- 1000 # number of repetition of simulations for each scenario
# prepare results data frame:
res <- data.frame(OSA_prob_intercept = rep(NA, length(OSA_prob_intercept)*length(OSA_misclassification_rate)*n_rep), 
                        OSA_misclassification_rate = NA, 
                        OSA_true_est = NA,
                        pval_true_est = NA,
                        OSA_observed_est = NA,
                        pval_observed_est = NA)

row_ind <- 1
for (p in 1:length(OSA_prob_intercept)){
  for (m in 1:length(OSA_misclassification_rate)){
    for (i in 1:n_rep){
      
      res$OSA_prob_intercept[row_ind] <- OSA_prob_intercept[p]
      res$OSA_misclassification_rate[row_ind] <- OSA_misclassification_rate[m]
        
      # simulate genotype
      g <- rbinom(n, size = 2, prob = 0.3)
      
      # OSA probability for each person
      OSA_prob <- expit(OSA_prob_intercept[p] + 0.1*g)
      
      # OSA  based on probability
      OSA <- rbinom(n = length(OSA_prob), size = 1, prob = OSA_prob)
      
      # observed OSA status
      OSA_status <- ifelse(OSA== 0, 0, rbinom(sum(OSA), size = 1, prob = 1-OSA_misclassification_rate[m]))
      
      
      res$OSA_true_est[row_ind] <- summary(glm(OSA ~ g, family = "binomial"))$coef["g", "Estimate"]
      res$OSA_observed_est[row_ind] <- summary(glm(OSA_status ~ g, family = "binomial"))$coef["g", "Estimate"]
      
      res$pval_true_est[row_ind] <- summary(glm(OSA ~ g, family = "binomial"))$coef["g", "Pr(>|z|)"]
      res$pval_observed_est[row_ind] <- summary(glm(OSA_status ~ g, family = "binomial"))$coef["g", "Pr(>|z|)"]

      row_ind <- row_ind + 1
    }
  }
}
```



Summarize the results:

```r
res <- data.table(res)
summarized_res <- res[, .(mean_estimate_true = mean(OSA_true_est), 
                          mean_estimate_observed = mean(OSA_observed_est), 
                          mean_bias_true = mean(0.1 - OSA_true_est), 
                          mean_bias_observed = mean(0.1 - OSA_observed_est), 
                          power_true_0.05 = mean(pval_true_est < 0.05), 
                          power_observed_0.05 = mean(pval_observed_est < 0.05)), 
                           by = c("OSA_prob_intercept", 
                             "OSA_misclassification_rate")]


summarized_res
```

```
##    OSA_prob_intercept OSA_misclassification_rate mean_estimate_true
## 1:               -1.5                        0.4         0.09980589
## 2:               -1.5                        0.6         0.10048059
## 3:               -1.5                        0.8         0.09903014
## 4:               -1.0                        0.4         0.09995585
## 5:               -1.0                        0.6         0.09913497
## 6:               -1.0                        0.8         0.10125527
## 7:               -0.5                        0.4         0.09938142
## 8:               -0.5                        0.6         0.09958494
## 9:               -0.5                        0.8         0.09921285
##    mean_estimate_observed mean_bias_true mean_bias_observed power_true_0.05
## 1:             0.08963131   1.941062e-04         0.01036869           0.956
## 2:             0.08611771  -4.805897e-04         0.01388229           0.959
## 3:             0.08078696   9.698572e-04         0.01921304           0.951
## 4:             0.08602495   4.415443e-05         0.01397505           0.985
## 5:             0.07974199   8.650326e-04         0.02025801           0.980
## 6:             0.07719226  -1.255271e-03         0.02280774           0.984
## 7:             0.07896281   6.185813e-04         0.02103719           0.993
## 8:             0.07057452   4.150569e-04         0.02942548           0.991
## 9:             0.06517548   7.871530e-04         0.03482452           0.995
##    power_observed_0.05
## 1:               0.762
## 2:               0.568
## 3:               0.292
## 4:               0.843
## 5:               0.654
## 6:               0.380
## 7:               0.863
## 8:               0.640
## 9:               0.367
```

```r
write.csv(res, file = "estimates_OSA_and_misclassification.csv")
write.csv(summarized_res, file = "summarized_estimates_OSA_and_misclassification.csv")
```
