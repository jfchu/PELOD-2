AgeMonths <- function(age.days) {
  return (age.days / 365 * 12)
}

#Global vector of age cutoffs measured in months
ages <<- c(144, 60, 24, 12, 1)

FindAgeCutoffs <- function(X, age) {
  if (age >= ages[X]) 
    return (X)
  if (age < ages[5])
    return (6)
}

#Logical values for score functions -> if clinical value should be less than a certain threshold

#Glasgow Coma Score
#True
GComa <- function(id) {
  cutoffs <- list(c(4, 4), c(10, 1))
  min.val <- FindExtremeValueMin(id, pelod2.datalist[["GCS"]])
  FindCutoffs <- function(X, score) {
    if (is.na(score)) {
      return (NA)
    } else if (score <= cutoffs[[X]][1]) {
      return (cutoffs[[X]][2])
    }
    return (0)
  }
  return (max(sapply(seq_len(length(cutoffs)), FindCutoffs, score = min.val)))
}

#Pupil reactivity should be entered as follows: 1 if true or 0 if false
#True
Pupillary <- function(id) {
  cutoffs <- list(c(0, 5))
  min.val1 <- FindExtremeValuePupil(id, pelod2.datalist[["pup.left"]])
  min.val2 <- FindExtremeValuePupil(id, pelod2.datalist[["pup.right"]])
  if (is.na(min.val1) || is.na(min.val2))
    return (NA)
  else if (min.val1 <= cutoffs[[1]][1] && min.val2 <= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0)
}

#Molarity measured in mmol/L
#False
#lac.wholeblood contains data with ">" character 
Lactatemia <- function(id) {
  cutoffs <- list(c(11.0, 4), c(5.0, 1))
  max.val1 <- FindExtremeValueMax(id, pelod2.datalist[["lac"]])
  if (!is.numeric(as.numeric(FindValues(id, pelod2.datalist[["lac.wholeblood"]])))) {
    return (4)
  } else {
    max.val2 <- FindExtremeValueMax(id, pelod2.datalist[["lac.wholeblood"]])
  }
  FindCutoffs <- function(X, molarity) {
    if (is.na(molarity))
      return (NA)
    else if (molarity >= cutoffs[[X]][1])
      return (cutoffs[[X]][2])
    return (0)
  }
  if (is.na(max.val1) && is.na(max.val2))
    return (NA)
  else if (is.na(max.val2) || (!is.na(max.val1) && max.val1 >= max.val2)) 
    return (max(sapply(seq_len(length(cutoffs)), FindCutoffs, molarity = max.val1)))
  else
    return (max(sapply(seq_len(length(cutoffs)), FindCutoffs, molarity = max.val2)))
}

#Age measured in months
#Mean arterial pressure measure in mm Hg
#True
MAP <- function(id, age) {
  cutoffs <- list(list(c(37, 6), c(51, 3), c(66, 2)), list(c(35, 6), c(48, 3), c(64, 2)), list(c(31, 6), c(44, 3), c(61, 2)), list(c(30, 6), c(43, 3), c(59, 2)), list(c(24, 6), c(38, 3), c(54, 2)), list(c(16, 6), c(30, 3), c(45, 2)))
  min.val <- FindExtremeValueMin(id, pelod2.datalist[["map"]])
  age.index <- min(unlist(sapply(seq_len(length(ages)), FindAgeCutoffs, age = age), recursive = F))
  FindCutoffs <- function(X, pressure) {
    if (is.na(pressure))
      return (NA)
    else if (pressure <= cutoffs[[age.index]][[X]][1])
      return (cutoffs[[age.index]][[X]][2])
    return (0)
  }
  return (max(sapply(seq_len(length(cutoffs[[1]])), FindCutoffs, pressure = min.val)))
}

#Age measured in months
#Molarity measured in μmol/L
#False
#Conversion factor from raw data to units used in score calculation is 88.44
#cr contains data with "<" character
Creatinine <- function(id, age) {
  cutoffs <- list(list(c(93, 2)), list(c(59, 2)), list(c(51, 2)), list(c(35, 2)), list(c(23, 2)), list(c(70, 2)))
  max.val <- FindExtremeValueMax(id, pelod2.datalist[["cr"]]) * 88.44
  age.index <- min(unlist(sapply(seq_len(length(ages)), FindAgeCutoffs, age = age), recursive = F))
  FindCutoffs <- function(X, cr.molarity) {
    if (is.na(cr.molarity))
      return (NA)
    else if (cr.molarity >= cutoffs[[age.index]][[X]][1])
      return (cutoffs[[age.index]][[X]][2])
    return (0)
  }
  return (max(sapply(seq_len(length(cutoffs[[1]])), FindCutoffs, cr.molarity = max.val)))
}

#Carrico index = PaO2(measured in mm Hg)/FIO2
#PaO2 -> partial pressure arterial oxygen
#FiO2 -> fraction of inspired oxygen
#True
#Finds min value of pao2 and fio2 value with closest timestamp (before or after)
#FiO2 data measured in percent in file
Carrico <- function(id) {
  cutoffs <- list(c(60, 2))
  pao2.frame <- pelod2.datalist[["pao2"]][which(pelod2.datalist[["pao2"]]$encid == id), ]
  fio2.frame <- pelod2.datalist[["fio2"]][which(pelod2.datalist[["fio2"]]$encid == id), ]
  min.val <- FindExtremeValueMin(id, pao2.frame)
  min.val.time <- as.Date(pao2.frame[match(min.val, pao2.frame[, "Clinical.Event.Result"]), "Clinical.Events.Verified.Date/Time"], origin="1970-01-01")
  FindClosestTime <- function(t) {
    if(identical(nrow(fio2.frame), 0L))
      return (NA)
    vals <- as.Date(fio2.frame[, "Clinical.Events.Verified.Date/Time"], origin="1970-01-01")
    closest.time.index <- which.closest(vals, t)
    return (closest.time.index)
  }
  max.val <- as.numeric(FindValues(id, fio2.frame))[FindClosestTime(min.val.time)] / 100
  if (is.na(min.val) || is.na(max.val))
    return (NA)
  else if (min.val / max.val <= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0) 
}

#PaCO2(measured in mm HG) -> partial pressure arterial carbon dioxide
#True
#paco2 contains data with ">" character
PaCO2 <- function(id) {
  cutoffs <- list(c(95, 3), c(59, 1))
  if (!is.numeric(as.numeric(FindValues(id, pelod2.datalist[["paco2"]]))))
    return (3)
  max.val <- FindExtremeValueMax(id, pelod2.datalist[["paco2"]])
  FindCutoffs <- function(X, partial.pressure) {
    if (is.na(partial.pressure))
      return (NA)
    else if (partial.pressure >= cutoffs[[X]][1])
      return (cutoffs[[X]][2])
    return (0)
  }
  return (max(sapply(seq_len(length(cutoffs)), FindCutoffs, partial.pressure = max.val)))
}

#Presence of invasive ventilation should be entered as follows: 1 if true or 0 if false
#False
Ventilation <- function(id) {
  cutoffs <- list(c(1, 3))
  max.val <- FindExtremeValueVent(id, pelod2.datalist[["vent"]])
  if (is.na(max.val))
    return (NA)
  else if (max.val >= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0)
}

#wbc.count * 10^9/L
#True
#wbc contains data with "<" character 
WBC <- function(id) {
  cutoffs <- list(c(2, 2))
  if (!is.numeric(as.numeric(FindValues(id, pelod2.datalist[["wbc"]])))) {
    return (2)
  }
  min.val = FindExtremeValueMin(id, pelod2.datalist[["wbc"]])
  if (is.na(min.val))
    return (NA)
  else if (min.val <= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0)
}

#plate.count * 10^9/L
#True
Platelet <- function(id) {
  cutoffs <- list(c(76, 2), c(141, 1))
  min.val = FindExtremeValueMin(id, pelod2.datalist[["plate"]])
  FindCutoffs <- function(X, plate.count) {
    if (is.na(plate.count))
      return (NA)
    else if (plate.count <= cutoffs[[X]][1])
      return (cutoffs[[X]][2])
    return (0)
  }
  return (max(sapply(seq_len(length(cutoffs)), FindCutoffs, plate.count = min.val)))
}

ProbMortality <- function(pelod2.scores) {
  return (1 / (1 + exp(6.61 - 0.47 * pelod2.scores)))
}

ProbMortalityNew <- function(pelod2) {
  return (1 / (1 + exp((-6.76204 + 0.3904402 * pelod2$pelod2.gcs + 0.5149909 * pelod2$pelod2.pup + 0.1381793 * pelod2$pelod2.map - 0.06871416 * pelod2$pelod2.cr + 0.2626874 * pelod2$pelod2.carrico + 0.8777295 * pelod2$pelod2.paco2 + 0.1311851 * pelod2$pelod2.vent + 0.7447805 * pelod2$pelod2.wbc + 0.2743563 * pelod2$pelod2.plate) * -1)))
}

##Logistic regression results on all subscores with full data
#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-1.9376  -0.0873  -0.0469  -0.0320   3.8929  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  -7.57695    0.38917 -19.469  < 2e-16 ***
#  neuroscores   0.45135    0.04675   9.653  < 2e-16 ***
#  cardioscores  0.27382    0.05749   4.763 1.91e-06 ***
#  renalscores   0.01417    0.14363   0.099    0.921    
#respscores    0.38281    0.08386   4.565 5.00e-06 ***
#  hemascores    0.49196    0.09600   5.125 2.98e-07 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 1054.78  on 5117  degrees of freedom
#Residual deviance:  502.55  on 5112  degrees of freedom
#AIC: 514.55

#Number of Fisher Scoring iterations: 8
ProbMortalityNewSubscores <- function(pelod2) {
  return (1 / (1 + exp((-7.576947 + 0.4513461 * pelod2$neuroscores + 0.2738175 * pelod2$cardioscores + 0.01417186 * pelod2$renalscores + 0.38281 * pelod2$respscores + 0.4919619 * pelod2$hemascores) * -1)))
}

##Hosmer-Lemeshow results on full dataset after downsampling
#$Table_HLtest
#total meanpred meanobs predicted observed
#[0.0178,0.0268)    58    0.025   0.000      1.43        0
#[0.0268,0.0275)    16    0.027   0.000      0.43        0
#[0.0275,0.0331)    46    0.032   0.000      1.45        0
#0.0331             31    0.033   0.000      1.03        0
#[0.0338,0.0392)    19    0.037   0.000      0.71        0
#0.0392            270    0.039   0.000     10.58        0
#0.0397            337    0.040   0.000     13.39        0
#[0.0406,0.0430)    50    0.042   0.000      2.08        0
#[0.0430,0.0456)    18    0.044   0.000      0.80        0
#0.0456             43    0.046   0.000      1.96        0
#[0.0473,0.0473)     6    0.047   0.167      0.28        1
#[0.0473,0.0494)    24    0.048   0.000      1.16        0
#[0.0494,0.0503)   110    0.050   0.000      5.51        0
#[0.0503,0.0567)    12    0.053   0.000      0.63        0
#[0.0567,0.0577)    26    0.057   0.000      1.48        0
#[0.0577,0.0584)    65    0.058   0.000      3.78        0
#[0.0584,0.0603)  1658    0.060   0.001     99.40        1
#[0.0603,0.0678)    14    0.064   0.000      0.89        0
#0.0678             78    0.068   0.000      5.29        0
#[0.0681,0.0695)    25    0.069   0.000      1.71        0
#0.0695            417    0.069   0.005     28.97        2
#[0.0698,0.0724)    18    0.072   0.000      1.29        0
#[0.0724,0.0759)    18    0.073   0.000      1.32        0
#[0.0759,0.0801)    28    0.078   0.000      2.19        0
#[0.0801,0.0849)    17    0.083   0.000      1.41        0
#[0.0849,0.0872)    38    0.086   0.000      3.29        0
#[0.0872,0.1034)   268    0.103   0.000     27.59        0
#[0.1034,0.1079)     8    0.104   0.000      0.83        0
#[0.1079,0.1171)    41    0.115   0.024      4.72        1
#[0.1171,0.1234)    10    0.121   0.000      1.21        0
#[0.1234,0.1752)    32    0.152   0.000      4.87        0
#[0.1752,0.1882)    14    0.179   0.000      2.51        0
#[0.1882,0.2090)    24    0.192   0.000      4.62        0
#[0.2090,0.2200)    82    0.218   0.000     17.87        0
#[0.2200,0.2382)    21    0.228   0.048      4.79        1
#[0.2382,0.2544)    18    0.243   0.000      4.37        0
#[0.2544,0.2670)    55    0.263   0.018     14.49        1
#[0.2670,0.2934)    26    0.288   0.038      7.48        1
#[0.2934,0.3060)   275    0.305   0.011     83.81        3
#0.3060              8    0.306   0.000      2.45        0
#[0.3071,0.3417)    58    0.332   0.017     19.23        1
#[0.3417,0.3491)    16    0.346   0.000      5.54        0
#[0.3491,0.3936)    20    0.371   0.050      7.43        1
#[0.3936,0.4037)    57    0.403   0.035     22.95        2
#[0.4037,0.4378)    17    0.422   0.000      7.18        0
#[0.4378,0.4689)    24    0.452   0.042     10.86        1
#[0.4689,0.5015)    24    0.490   0.000     11.75        0
#[0.5015,0.5472)    23    0.522   0.043     12.00        1
#[0.5472,0.5652)    26    0.555   0.000     14.43        0
#[0.5652,0.5986)    47    0.591   0.000     27.79        0
#[0.5986,0.6510)    25    0.618   0.080     15.46        2
#[0.6510,0.6675)    23    0.652   0.000     15.00        0
#[0.6675,0.7063)    27    0.684   0.074     18.47        2
#[0.7063,0.7806)    21    0.734   0.095     15.41        2
#[0.7806,0.7925)    25    0.789   0.000     19.74        0
#[0.7925,0.8331)    26    0.816   0.038     21.20        1
#[0.8331,0.8591)    50    0.852   0.020     42.59        1
#[0.8591,0.8784)    28    0.872   0.036     24.40        1
#[0.8784,0.8934)    17    0.888   0.000     15.09        0
#[0.8934,0.9188)    24    0.903   0.167     21.67        4
#[0.9188,0.9333)    23    0.927   0.087     21.33        2
#[0.9333,0.9562)    25    0.944   0.200     23.61        5
#[0.9562,0.9663)    25    0.961   0.440     24.04       11
#[0.9663,0.9761)    27    0.972   0.222     26.25        6
#[0.9761,0.9817)    20    0.979   0.350     19.59        7
#[0.9817,0.9889)    27    0.985   0.296     26.60        8
#[0.9889,0.9949)    21    0.992   0.286     20.82        6
#[0.9949,0.9994)    24    0.998   0.542     23.95       13
#[0.9994,1.0000]    24    1.000   0.875     23.99       21

#$Chi_square
#[1] Inf

#$df
#[1] 210

#$p_value
#[1] 0
ProbMortalityDownsample <- function(pelod2) {
  return (1 / (1 + exp((-2.75234 + 0.5901375 * pelod2$pelod2.gcs + 0.6751759 * pelod2$pelod2.pup + 0.2286607 * pelod2$pelod2.lac + 0.0653454 * pelod2$pelod2.map - 0.223401 * pelod2$pelod2.cr + 0.5362901 * pelod2$pelod2.carrico + 1.081174 * pelod2$pelod2.paco2 - 0.1441723 * pelod2$pelod2.vent + 2.450229 * pelod2$pelod2.wbc - 0.1907878 * pelod2$pelod2.plate) * -1)))
}

##Hosmer-Lemeshow results of subscores after downsampling
#$Table_HLtest
#total meanpred meanobs predicted observed
#[0.0258,0.0307)    32    0.028   0.031      0.89        1
#[0.0307,0.0617)    20    0.052   0.000      1.04        0
#[0.0617,0.0967)    22    0.090   0.136      1.98        3
#[0.0967,0.3014)    15    0.172   0.200      2.57        3
#[0.3014,0.4527)    20    0.348   0.400      6.95        8
#[0.4527,0.7436)    22    0.622   0.500     13.69       11
#[0.7436,0.9106)    22    0.855   0.818     18.81       18
#[0.9106,0.9655)    22    0.945   1.000     20.80       22
#[0.9655,0.9868)    22    0.974   1.000     21.43       22
#[0.9868,0.9983]    21    0.993   1.000     20.84       21

#$Chi_square
#[1] 5.675

#$df
#[1] 8

#$p_value
#[1] 0.6836

##Hosmer-Lemeshow results on full dataset of subscores after downsampling
#$Table_HLtest
#total meanpred meanobs predicted observed
#0.0258            270    0.026   0.000      6.95        0
#[0.0279,0.0301)  1656    0.028   0.001     46.93        1
#0.0301              9    0.030   0.000      0.27        0
#[0.0307,0.0359)    86    0.033   0.000      2.83        0
#0.0359             13    0.036   0.000      0.47        0
#[0.0360,0.0453)    29    0.042   0.000      1.21        0
#[0.0453,0.0482)   268    0.046   0.000     12.42        0
#0.0482             32    0.048   0.000      1.54        0
#[0.0493,0.0521)     3    0.050   0.000      0.15        0
#[0.0521,0.0529)    46    0.053   0.000      2.42        0
#0.0529            109    0.053   0.000      5.77        0
#[0.0532,0.0562)    41    0.054   0.000      2.21        0
#[0.0562,0.0584)   343    0.058   0.000     19.79        0
#[0.0584,0.0622)    22    0.061   0.000      1.34        0
#[0.0622,0.0767)    21    0.069   0.048      1.44        1
#[0.0767,0.0852)    56    0.083   0.000      4.67        0
#0.0852             24    0.085   0.000      2.04        0
#[0.0857,0.0906)    18    0.088   0.000      1.59        0
#[0.0906,0.0952)   421    0.093   0.005     39.00        2
#[0.0952,0.0967)    10    0.096   0.000      0.96        0
#0.0967             55    0.097   0.000      5.32        0
#[0.0972,0.1034)    17    0.099   0.000      1.68        0
#[0.1034,0.1060)    35    0.105   0.000      3.67        0
#[0.1060,0.1105)    33    0.107   0.000      3.53        0
#[0.1105,0.1152)     8    0.111   0.000      0.89        0
#[0.1152,0.1388)    23    0.120   0.000      2.76        0
#[0.1388,0.1567)    24    0.148   0.000      3.56        0
#[0.1567,0.1693)    74    0.163   0.014     12.04        1
#[0.1693,0.1754)    22    0.171   0.000      3.77        0
#[0.1754,0.1854)    67    0.183   0.030     12.29        2
#[0.1854,0.1916)     9    0.187   0.000      1.68        0
#[0.1916,0.2148)    21    0.204   0.000      4.29        0
#[0.2148,0.2726)    28    0.255   0.000      7.13        0
#0.2726             19    0.273   0.000      5.18        0
#[0.2730,0.2855)    63    0.279   0.000     17.59        0
#[0.2855,0.3010)    78    0.299   0.000     23.36        0
#[0.3010,0.3023)     4    0.301   0.250      1.20        1
#[0.3023,0.3234)   279    0.321   0.011     89.49        3
#[0.3234,0.3391)    19    0.333   0.000      6.32        0
#[0.3391,0.3608)    64    0.357   0.016     22.82        1
#[0.3608,0.3850)    24    0.375   0.000      9.01        0
#[0.3850,0.4179)    49    0.407   0.020     19.96        1
#[0.4179,0.4409)    24    0.430   0.083     10.33        2
#[0.4409,0.4525)    25    0.449   0.000     11.23        0
#[0.4525,0.4771)    50    0.473   0.020     23.67        1
#[0.4771,0.5116)    22    0.487   0.045     10.72        1
#[0.5116,0.5167)    33    0.515   0.030     16.99        1
#[0.5167,0.5372)    17    0.534   0.000      9.08        0
#[0.5372,0.5549)    22    0.543   0.000     11.95        0
#[0.5549,0.5720)    27    0.566   0.074     15.29        2
#[0.5720,0.5795)    22    0.576   0.045     12.68        1
#[0.5795,0.6129)    28    0.602   0.036     16.85        1
#[0.6129,0.6361)    23    0.630   0.043     14.50        1
#[0.6361,0.6714)    26    0.658   0.038     17.12        1
#[0.6714,0.7083)    20    0.693   0.050     13.86        1
#[0.7083,0.7385)    27    0.724   0.037     19.55        1
#[0.7385,0.7688)    19    0.748   0.053     14.20        1
#[0.7688,0.8215)    25    0.794   0.080     19.86        2
#[0.8215,0.8602)    38    0.853   0.105     32.43        4
#[0.8602,0.8650)    10    0.864   0.100      8.64        1
#[0.8650,0.8841)    25    0.877   0.120     21.93        3
#[0.8841,0.9002)    23    0.891   0.261     20.50        6
#[0.9002,0.9310)    25    0.915   0.240     22.87        6
#[0.9310,0.9572)    23    0.945   0.565     21.74       13
#[0.9572,0.9729)    24    0.965   0.625     23.17       15
#[0.9729,0.9887)    25    0.981   0.520     24.52       13
#[0.9887,0.9983]    23    0.993   0.826     22.84       19

#$Chi_square
#[1] 2179.25

#$df
#[1] 210

#$p_value
#[1] 0
ProbMortalityDownsampleSubscores <- function(pelod2) {
  return (1 / (1 + exp((-3.534816 + 0.5110541 * pelod2$neuroscores + 0.0811069 * pelod2$cardioscores - 0.04917815 * pelod2$renalscores + 0.2473746 * pelod2$respscores + 0.6501 * pelod2$hemascores) * -1)))
}

##Logistic regression with variables selected from all subset selection with glmulti package
#Coefficients:
#  (Intercept)    pelod2.gcs    pelod2.pup    pelod2.lac  
#-6.5677        0.5135        0.5271        0.6153  
#pelod2.paco2    pelod2.wbc  
#1.0866        0.8862  

#Degrees of Freedom: 5117 Total (i.e. Null);  5112 Residual
#Null Deviance:	    1055 
#Residual Deviance: 494.9 	AIC: 506.9


##Hosmer-Lemeshow results on full dataset after variable selection with glmulti
#$Table_HLtest
#total meanpred meanobs predicted observed
#0.00140            2722    0.001   0.000      3.82        1
#0.00234            1005    0.002   0.002      2.35        2
#0.00259              23    0.003   0.000      0.06        0
#[0.00415,0.00691)    33    0.004   0.061      0.14        2
#0.00691              20    0.007   0.000      0.14        0
#[0.00765,0.01084)   105    0.008   0.010      0.86        1
#0.01084             665    0.011   0.015      7.21       10
#[0.01271,0.01507)    19    0.014   0.053      0.26        1
#[0.01507,0.01987)    39    0.018   0.000      0.71        0
#0.01987              21    0.020   0.000      0.42        0
#[0.02493,0.03171)   129    0.031   0.031      4.02        4
#0.03171              41    0.032   0.049      1.30        2
#[0.03499,0.05713)    35    0.055   0.057      1.92        2
#[0.05713,0.08848)    32    0.060   0.031      1.94        1
#[0.08848,0.11382)    14    0.104   0.071      1.46        1
#[0.11382,0.16048)    59    0.129   0.102      7.63        6
#0.16048              14    0.160   0.000      2.25        0
#[0.22046,0.26128)    26    0.222   0.231      5.76        6
#[0.26128,0.34556)    36    0.294   0.389     10.58       14
#[0.34556,0.47354)    17    0.426   0.471      7.23        8
#[0.47354,0.69141)    17    0.587   0.647      9.97       11
#[0.69141,0.88045)    30    0.786   0.733     23.57       22
#[0.88045,0.99637]    16    0.962   0.938     15.39       15

#$Chi_square
#[1] 40.7

#$df
#[1] 210

#$p_value
#[1] 1

##Cross validation results for all subsets regression on individual scores
##Results for 2558 rows
#Confusion Matrix and Statistics

#Reference
#Prediction    0    1
#0 2498   30
#1    6   24

#Accuracy : 0.9859          
#95% CI : (0.9806, 0.9901)
#No Information Rate : 0.9789          
#P-Value [Acc > NIR] : 0.0056725       

#Kappa : 0.5649          
#Mcnemar's Test P-Value : 0.0001264       

#Sensitivity : 0.9976          
#Specificity : 0.4444          
#Pos Pred Value : 0.9881          
#Neg Pred Value : 0.8000          
#Prevalence : 0.9789          
#Detection Rate : 0.9765          
#Detection Prevalence : 0.9883          
#Balanced Accuracy : 0.7210          

#'Positive' Class : 0

##Results for 2560 rows
#Confusion Matrix and Statistics

#Reference
#Prediction    0    1
#0 2499   33
#1    6   22

#Accuracy : 0.9848          
#95% CI : (0.9792, 0.9891)
#No Information Rate : 0.9785          
#P-Value [Acc > NIR] : 0.01389         

#Kappa : 0.5232          
#Mcnemar's Test P-Value : 3.136e-05       

#Sensitivity : 0.9976          
#Specificity : 0.4000          
#Pos Pred Value : 0.9870          
#Neg Pred Value : 0.7857          
#Prevalence : 0.9785          
#Detection Rate : 0.9762          
#Detection Prevalence : 0.9891          
#Balanced Accuracy : 0.6988          

#'Positive' Class : 0
ProbMortalityGlmulti <- function(pelod2) {
  return (1 / (1 + exp((-6.567728 + 0.5135225 * pelod2$pelod2.gcs + 0.5270587 * pelod2$pelod2.pup + 0.6153427 * pelod2$pelod2.lac + 1.086556 * pelod2$pelod2.paco2 + 0.8862088 * pelod2$pelod2.wbc) * -1)))
}

##Logistic regression with subscores selected with all subset regression with glmulti package
#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-1.9365  -0.0866  -0.0470  -0.0320   3.8929  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)         -7.57671    0.38920 -19.467  < 2e-16 ***
#  neuroscores   0.45210    0.04613   9.800  < 2e-16 ***
#  cardioscores  0.27501    0.05624   4.890 1.01e-06 ***
#  respscores    0.38367    0.08345   4.598 4.27e-06 ***
#  hemascores    0.49463    0.09209   5.371 7.83e-08 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 1054.78  on 5117  degrees of freedom
#Residual deviance:  502.56  on 5113  degrees of freedom
#AIC: 512.56

#Number of Fisher Scoring iterations: 8


##Hosmer-Lemeshow results after all subset regression with glmulti on subscores
#$Table_HLtest
#total meanpred meanobs predicted observed
#0.000512             1925    0.001   0.001      0.99        1
#0.000674                8    0.001   0.000      0.01        0
#[0.000751,0.000839)   285    0.001   0.000      0.23        0
#0.000839              141    0.001   0.000      0.12        0
#0.000887               87    0.001   0.000      0.08        0
#[0.001059,0.001231)    24    0.001   0.000      0.03        0
#[0.001231,0.001376)    35    0.001   0.000      0.05        0
#0.001376               72    0.001   0.000      0.10        0
#0.001393               38    0.001   0.000      0.05        0
#0.001454               19    0.001   0.000      0.03        0
#0.001617              380    0.002   0.000      0.61        0
#[0.001713,0.001810)     4    0.002   0.000      0.01        0
#[0.001810,0.002254)    31    0.002   0.000      0.06        0
#[0.002254,0.002371)    19    0.002   0.000      0.04        0
#[0.002371,0.002649)   470    0.003   0.006      1.19        3
#0.002649               40    0.003   0.000      0.11        0
#[0.002660,0.003115)    21    0.003   0.000      0.06        0
#0.003115               55    0.003   0.036      0.17        2
#[0.003134,0.003691)    16    0.003   0.000      0.05        0
#0.003691               58    0.004   0.000      0.21        0
#[0.003721,0.004112)    14    0.004   0.000      0.05        0
#[0.004112,0.004174)    72    0.004   0.000      0.30        0
#0.004174                3    0.004   0.000      0.01        0
#[0.004336,0.004582)    42    0.004   0.000      0.18        0
#0.004582                5    0.005   0.000      0.02        0
#[0.004854,0.004916)    25    0.005   0.000      0.12        0
#[0.004916,0.005788)    26    0.005   0.077      0.14        2
#[0.005788,0.007091)    34    0.007   0.000      0.23        0
#[0.007091,0.007492)    13    0.007   0.000      0.09        0
#[0.007492,0.009783)    23    0.008   0.043      0.19        1
#0.009783              342    0.010   0.009      3.35        3
#[0.009840,0.012355)    22    0.011   0.045      0.25        1
#[0.012355,0.015289)    22    0.014   0.000      0.30        0
#[0.015289,0.016835)    66    0.016   0.015      1.05        1
#0.016835               67    0.017   0.015      1.13        1
#[0.017091,0.019083)    11    0.018   0.091      0.20        1
#[0.019083,0.022358)    25    0.021   0.000      0.54        0
#[0.022358,0.024516)    35    0.024   0.029      0.83        1
#0.024516               40    0.025   0.000      0.98        0
#0.025881               28    0.026   0.036      0.72        1
#0.027314               28    0.027   0.036      0.76        1
#[0.028824,0.030916)    13    0.030   0.154      0.39        2
#[0.030916,0.035652)    26    0.034   0.038      0.88        1
#[0.035652,0.039582)    24    0.036   0.000      0.87        0
#[0.039582,0.044818)    34    0.042   0.059      1.42        2
#[0.044818,0.051466)    20    0.050   0.000      0.99        0
#[0.051466,0.057161)    17    0.053   0.000      0.90        0
#[0.057161,0.066671)    25    0.062   0.000      1.54        0
#[0.066671,0.077798)    23    0.068   0.043      1.57        1
#[0.077798,0.090241)    41    0.084   0.073      3.46        3
#[0.090241,0.099775)     8    0.092   0.000      0.74        0
#[0.099775,0.122052)    27    0.109   0.037      2.95        1
#[0.122052,0.154073)    26    0.139   0.000      3.62        0
#[0.154073,0.185653)    21    0.169   0.190      3.54        4
#[0.185653,0.239561)    23    0.210   0.217      4.83        5
#[0.239561,0.342203)    24    0.289   0.333      6.94        8
#[0.342203,0.460368)    23    0.392   0.435      9.02       10
#[0.460368,0.607608)    24    0.521   0.667     12.51       16
#[0.607608,0.779419)    24    0.676   0.708     16.23       17
#[0.779419,0.986482]    24    0.873   0.833     20.96       20

#$Chi_square
#[1] 84.925

#$df
#[1] 210

#$p_value
#[1] 1
ProbMortalityGlmultiSubscores <- function(pelod2) {
  return (1 / (1 + exp((-7.57671 + 0.452099 * pelod2$neuroscores + 0.2750051 * pelod2$cardioscores + 0.3836698 * pelod2$respscores + 0.4946317 * pelod2$hemascores) * -1)))
}

#poor model generated with linear regression
ProbMortalityGlmulti2 <- function(pelod2) {
  return (1 / (1 + exp((-0.005095304 + 0.0339569 * pelod2$pelod2.pup + 0.08460044 * pelod2$pelod2.lac + 0.03506476 * pelod2$pelod2.carrico + 0.09284607 * pelod2$pelod2.paco2 + 0.03479897 * pelod2$pelod2.wbc) * -1)))
}

AllSubsetReg <- function() {
  return (glmulti(deceased~pelod2.gcs+pelod2.pup+pelod2.lac+pelod2.map+pelod2.cr+pelod2.carrico+pelod2.paco2+pelod2.vent+pelod2.wbc+pelod2.plate, intercept = T, level = 1, data = PELOD2Scores(pelod2.datalist), crit = "bic", fitfunction = "glm", family = binomial, confsetsize = 10))
}

AllSubsetRegSubscores <- function() {
  return (glmulti(deceased~neuroscores+cardioscores+renalscores+respscores+hemascores, intercept = T, level = 1, data = PELOD2Scores(pelod2.datalist), crit = "bic", fitfunction = "glm", confsetsize = 10, family = binomial))
}

FindValues <- function(id, frame) {
  vals <- frame[which(frame$encid == id), "Clinical.Event.Result"]
  if(identical(nrow(vals), 0L))
    return(NA)
  if (all(is.na(vals)))
    return (NA)
  else
    return (vals[complete.cases(vals)])
}

FindExtremeValueMax <- function(id, frame) {
  vals <- as.numeric(frame[which(frame$encid == id), "Clinical.Event.Result"])
  if (all(is.na(vals)))
    return (NA)
  else
    return (max(vals, na.rm = T))
}

FindExtremeValueMin <- function(id, frame) {
  vals <- as.numeric(frame[which(frame$encid == id), "Clinical.Event.Result"])
  if (all(is.na(vals)))
    return (NA)
  else
    return (min(vals, na.rm = T))
}

FindExtremeValuePupil <- function(id, frame) {
  vals <- frame[which(frame$encid == id), "Clinical.Event.Result"]
  FilterNonreactive <- function(X) {if (X == "Nonreactive") return (0) else if (X == "Unable to assess") return (NA) else return (1)}
  if (all(is.na(sapply(vals, FilterNonreactive))))
    return (NA)
  else
    return (min(unlist(sapply(vals, FilterNonreactive), recursive = F), na.rm = T))
}

FindExtremeValueVent <- function(id, frame) {
  vals <- frame[which(frame$encid == id), "Clinical.Event.Result"]
  FilterCharacter <- function(X) {if (is.character(X)) return (1) else return (0)}
  if (all(is.na(sapply(vals, FilterCharacter))))
    return (NA)
  else
    return (max(unlist(sapply(vals, FilterCharacter), recursive = F)))
}

FindDateVals <- function(id, date, frame) {
  start.time <- frame[min(which(frame$encid == id)), "Clinical.Event.Performed.Date/Time"]
  vals <- frame[which(frame$encid == id), ]
  vals <- vals[, "Clinical.Event.Result"][which(vals$"Clinical.Event.Performed.Date/Time" >= start.time + date - 1 & vals$"Clinical.Event.Performed.Date/Time" <= start.time + date)]
  return (vals)
}

FindExtremeValueMaxDate <- function(id, date, frame) {
  vals <- FindDateVals(id, date, frame)
  vals <- as.numeric(vals)
  if (all(is.na(vals)))
    return (NA)
  else
    return (max(vals, na.rm = T))
}

FindExtremeValueMinDate <- function(id, date, frame) {
  vals <- FindDateVals(id, date, frame)
  vals <- as.numeric(vals)
  if (all(is.na(vals)))
    return (NA)
  else
    return (min(vals, na.rm = T))
}

FindExtremeValuePupDate <- function(id, date, frame) {
  vals <- FindDateVals(id, date, frame) 
  FilterNonreactive <- function(X) {if (X == "Nonreactive") return (0) else if (X == "Unable to assess") return (NA) else return (1)}
  if (all(is.na(sapply(vals, FilterNonreactive))))
    return (NA)
  else
    return (min(unlist(sapply(vals, FilterNonreactive), recursive = F), na.rm = T))
}

FindExtremeValueVentDate <- function(id, date, frame) {
  vals <- FindDateVals(id, date, frame)
  FilterCharacter <- function(X) {if (is.character(X)) return (1) else return (0)}
  if (all(is.na(sapply(vals, FilterCharacter))))
    return (NA)
  else
    return (max(unlist(sapply(vals, FilterCharacter), recursive = F)))
}

PELOD2Scores <- function(frame.list) {
  pelod2.frame <- data.frame(
    pelod2.id = read.xlsx("PELOD2_calculator_deid_unencrypt.xlsx")$encid,
    pelod2.age = AgeMonths(read.xlsx("admit1picu_deid_unencrypt.xlsx")$AGEDAYS))
    pelod2.frame$pelod2.gcs = sapply(pelod2.frame$pelod2.id, GComa)
    pelod2.frame$pelod2.pup = sapply(pelod2.frame$pelod2.id, Pupillary)
    pelod2.frame$pelod2.lac = sapply(pelod2.frame$pelod2.id, Lactatemia)
    pelod2.frame$pelod2.map = mapply(MAP, pelod2.frame$pelod2.id, pelod2.frame$pelod2.age)
    pelod2.frame$pelod2.cr = mapply(Creatinine, pelod2.frame$pelod2.id, pelod2.frame$pelod2.age)
    pelod2.frame$pelod2.carrico = sapply(pelod2.frame$pelod2.id, Carrico)
    pelod2.frame$pelod2.paco2 = sapply(pelod2.frame$pelod2.id, PaCO2)
    pelod2.frame$pelod2.vent = sapply(pelod2.frame$pelod2.id, Ventilation)
    pelod2.frame$pelod2.wbc = sapply(pelod2.frame$pelod2.id, WBC)
    pelod2.frame$pelod2.plate = sapply(pelod2.frame$pelod2.id, Platelet)
    pelod2.frame[is.na(pelod2.frame)] <- 0
    pelod2.frame$neuroscores = pelod2.frame$pelod2.gcs + pelod2.frame$pelod2.pup
    pelod2.frame$cardioscores = pelod2.frame$pelod2.lac + pelod2.frame$pelod2.map
    pelod2.frame$renalscores = pelod2.frame$pelod2.cr
    pelod2.frame$respscores = pelod2.frame$pelod2.carrico + pelod2.frame$pelod2.paco2 + pelod2.frame$pelod2.vent
    pelod2.frame$hemascores = pelod2.frame$pelod2.wbc + pelod2.frame$pelod2.plate
    pelod2.frame$pelod2.scores = pelod2.frame$neuroscores + pelod2.frame$cardioscores + pelod2.frame$renalscores + pelod2.frame$respscores + pelod2.frame$hemascores
#consider whether to calculate total score by adding neuro, cardio, ... or adding gcs, pup, ...
#see whether NA values will influence above decision
#    pelod2.frame$pelod2.scores = pelod2.frame$pelod2.gcs + pelod2.frame$pelod2.pup + pelod2.frame$pelod2.lac + pelod2.frame$pelod2.map + pelod2.frame$pelod2.cr + pelod2.frame$pelod2.carrico + pelod2.frame$pelod2.paco2 + pelod2.frame$pelod2.vent + pelod2.frame$pelod2.wbc + pelod2.frame$pelod2.plate
    #pelod2.frame$pred = ProbMortalityNew(pelod2.frame)
    pelod2.frame$deceased = FindDeceased(read.xlsx("admit1picu_deid_unencrypt.xlsx"))
    return (pelod2.frame)
}

RemoveNA <- function(df) {
  return(df[complete.cases(df), ])
}

FindDeceased <- function(frame) {
  IsDeceased <- function(id, indices) {
    if (is.element(id, indices))
      return (1)
    else 
      return (0)
  }
  return (sapply(read.xlsx("PELOD2_calculator_deid_unencrypt.xlsx")$encid, IsDeceased, indices = which(frame$Discharge.Disposition == "Deceased")))
}

FindDifferences <- function(frame) {
  master <- read.xlsx("PELOD2_calculator_deid_unencrypt.xlsx")
  differences <- data.frame(
    gcs = frame$pelod2.gcs - master$GCS_Points,
    pup = frame$pelod2.pup - master$Pupils_Points,
    lac = frame$pelod2.lac - master$Lactate_Points,
    map = frame$pelod2.map - master$MAP_Points,
    cr = frame$pelod2.cr - master$Cr_Points,
    carrico = frame$pelod2.carrico - master$PaO2FiO2_Points_Near,
    paco2 = frame$pelod2.paco2 - master$PaCO2_Points,
    vent = frame$pelod2.vent - master$MV_Points,
    wbc = frame$pelod2.wbc - master$WBC_Points,
    plate = frame$pelod2.plate - master$Plt_Points
  )
  return (differences)
}

#Adapted code from https://github.com/joyofdata/joyofdata-articles/blob/master/roc-auc/plot_pred_type_distribution.R 
plot_pred_type_distribution <- function(df, threshold) {
  v <- rep(NA, nrow(df))
  v <- ifelse(df$pred >= threshold & df$deceased == 1, "TP", v)
  v <- ifelse(df$pred >= threshold & df$deceased == 0, "FP", v)
  v <- ifelse(df$pred < threshold & df$deceased == 1, "FN", v)
  v <- ifelse(df$pred < threshold & df$deceased == 0, "TN", v)
  
  df$pred_type <- v
  
  ggplot(data=df, aes(x=deceased, y=pred)) + 
    geom_violin(fill=rgb(1,1,1,alpha=0.6), color=NA) + 
    geom_jitter(aes(color=pred_type), alpha=0.6) +
    geom_hline(yintercept=threshold, color="red", alpha=0.6) +
    scale_color_discrete(name = "type") +
    labs(title=sprintf("Threshold at %.2f", threshold))
}

#Adapted code from https://github.com/joyofdata/joyofdata-articles/blob/master/roc-auc/calculate_roc.R
calculate_roc <- function(df, cost_of_fp, cost_of_fn, n) {
  tpr <- function(df, threshold) {
    sum(df$pred >= threshold & df$deceased == 1) / sum(df$deceased == 1)
  }
  
  fpr <- function(df, threshold) {
    sum(df$pred >= threshold & df$deceased == 0) / sum(df$deceased == 0)
  }
  
  cost <- function(df, threshold, cost_of_fp, cost_of_fn) {
    sum(df$pred >= threshold & df$deceased == 0) * cost_of_fp + 
      sum(df$pred < threshold & df$deceased == 1) * cost_of_fn
  }
  
  roc <- data.frame(threshold = seq(0,1,length.out=n), tpr=NA, fpr=NA)
  roc$tpr <- sapply(roc$threshold, function(th) tpr(df, th))
  roc$fpr <- sapply(roc$threshold, function(th) fpr(df, th))
  roc$cost <- sapply(roc$threshold, function(th) cost(df, th, cost_of_fp, cost_of_fn))
  
  return(roc)
}

#File names:
#admit1picu_deid_unencrypt.xlsx
#PELOD2_2015-2016_deid_unencrypt.xlsx
pelod2.datalist <- list(
  GCS = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "GCS")),
  pup.left = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Pupilary Reaction Left")),
  pup.right = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Pupilary Reaction Right")),
  lac = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Lactate")),
  lac.wholeblood = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Lactate Whole Blood")),
  map = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Mean Arterial Pressure")),
  cr = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Creatine")),
  pao2 = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "PaO2 Verified Date")),
  fio2 = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "FiO2 Verified Date")),
  paco2 = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "PaCo2")),
  vent = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Mechanical Vent")),
  wbc = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "WBC")),
  plate = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Platelets"))
)

id = read.xlsx("PELOD2_calculator_deid_unencrypt.xlsx")$encid
pelod2.raw <- data.frame(
  age = AgeMonths(read.xlsx("admit1picu_deid_unencrypt.xlsx")$AGEDAYS),
  GCS = sapply(id, FindExtremeValueMin, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "GCS"))),
  pup.left = sapply(id, FindExtremeValuePupil, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Pupilary Reaction Left"))),
  pup.right = sapply(id, FindExtremeValuePupil, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Pupilary Reaction Right"))),
  map = sapply(id, FindExtremeValueMin, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Mean Arterial Pressure"))),
  cr = sapply(id, FindExtremeValueMax, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Creatine"))),
  pao2 = sapply(id, FindExtremeValueMin, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "PaO2"))),
  fio2 = sapply(id, FindExtremeValueMax, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "FiO2"))),
  paco2 = sapply(id, FindExtremeValueMax, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "PaCo2"))),
  vent = sapply(id, FindExtremeValueVent, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Mechanical Vent"))),
  wbc = sapply(id, FindExtremeValueMin, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "WBC"))),
  plate = sapply(id, FindExtremeValueMin, frame = RemoveNA(read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", sheet = "Platelets")))
)
pelod2.raw[is.na(pelod2.raw)] <- 0

#LOOCV code borrowed from https://www.r-bloggers.com/predicting-creditability-using-logistic-regression-in-r-cross-validating-the-classifier-part-2-2/
LOOCV <- function(data, cv.model) {
acc <- NULL
for(i in 1:nrow(data))
{
  # Train-test splitting
  # 499 samples -> fitting
  # 1 sample -> testing
  train <- data[-i,]
  test <- data[i,]
  
  # Predict results
  results_prob <- predict(cv.model,subset(test,select=c(2:9)),type='response')
  
  # If prob > 0.5 then 1, else 0
  results <- ifelse(results_prob > 0.5,1,0)
  
  # Actual answers
  answers <- test$deceased
  
  # Calculate accuracy
  misClassificError <- mean(answers != results)
  
  # Collecting results
  acc[i] <- 1-misClassificError
}

  # Average accuracy of the model
  cat("mean accuracy: ", mean(acc))

  # Histogram of the model accuracy
  hist(acc,xlab='Accuracy',ylab='Freq',main='Accuracy LOOCV',
     col='cyan',border='blue',density=30)
}

LOOCV2 <- function(df) {
  for (k in 1:n_train) {
    train.data <- df[-k, ]
    test.data <- df[k, ]
    x <- train.data$x
    y <- train.data$y
    fitted_models <- apply(t(df), 2, function(degf) lm(y ~ ns(x, df = degf)))
    pred <- mapply(function(obj, degf) predict(obj, data.frame(x = test_xy$x)), fitted_models, df)
    loocv_tmp[k, ] <- (test_xy$y - pred)^2
  }
  loocv <- colMeans(loocv_tmp)
  
  plot(df, mse, type = "l", lwd = 2, col = gray(.4), ylab = "Prediction error",
       xlab = "Flexibilty (spline's degrees of freedom [log scaled])",
       main = "Leave-One-Out Cross-Validation", ylim = c(.1, .8), log = "x")
  lines(df, cv, lwd = 2, col = "steelblue2", lty = 2)
  lines(df, loocv, lwd = 2, col = "darkorange")
  legend(x = "topright", legend = c("Training error", "10-fold CV error", "LOOCV error"),
         lty = c(1, 2, 1), lwd = rep(2, 3), col = c(gray(.4), "steelblue2", "darkorange"),
         text.width = .3, cex = .85)
}

##Results of Hosmer-Lemeshow GOF test (original model)
#$Table_HLtest
#total meanpred meanobs predicted observed
#0.00135            1655    0.001   0.001      2.23        1
#0.00215             377    0.002   0.000      0.81        0
#0.00344             434    0.003   0.000      1.49        0
#0.00549             481    0.005   0.000      2.64        0
#0.00875             605    0.009   0.008      5.29        5
#0.01393             186    0.014   0.005      2.59        1
#0.02210             156    0.022   0.000      3.45        0
#0.03489             334    0.035   0.009     11.65        3
#0.05468              94    0.055   0.032      5.14        3
#0.08471             195    0.085   0.026     16.52        5
#0.12898             108    0.129   0.028     13.93        3
#0.19155              90    0.192   0.022     17.24        2
#0.27488              82    0.275   0.024     22.54        2
#0.37754              56    0.378   0.036     21.14        2
#0.49250              53    0.493   0.019     26.10        1
#0.60826              32    0.608   0.031     19.46        1
#0.71300              28    0.713   0.143     19.96        4
#0.79899              28    0.799   0.250     22.37        7
#0.86413              25    0.864   0.200     21.60        5
#0.91052              15    0.911   0.333     13.66        5
#0.94213              19    0.942   0.579     17.90       11
#[0.96303,0.98523)    26    0.970   0.615     25.23       16
#[0.98523,0.99635)    23    0.990   0.739     22.78       17
#[0.99635,0.99986]    16    0.998   0.938     15.97       15

#$Chi_square
#[1] 786.73

#$df
#[1] 210

#$p_value
#[1] 0


##Results of logistic regression on Children's PICU data (updated model)
#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-2.4143  -0.0936  -0.0584  -0.0481   3.6778  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)    -6.76204    0.42801 -15.799  < 2e-16 ***
#  pelod2.gcs      0.39044    0.13843   2.821  0.00479 ** 
#  pelod2.pup      0.51499    0.05739   8.973  < 2e-16 ***
#  pelod2.lac      0.50890    0.10605   4.799 1.60e-06 ***
#  pelod2.map      0.13818    0.09022   1.532  0.12565    
#pelod2.cr      -0.06871    0.15618  -0.440  0.65996    
#pelod2.carrico  0.26269    0.17666   1.487  0.13703    
#pelod2.paco2    0.87773    0.16773   5.233 1.67e-07 ***
#  pelod2.vent     0.13119    0.19437   0.675  0.49973    
#pelod2.wbc      0.74478    0.17899   4.161 3.17e-05 ***
#  pelod2.plate    0.27436    0.19434   1.412  0.15803    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 1054.78  on 5117  degrees of freedom
#Residual deviance:  484.51  on 5107  degrees of freedom
#AIC: 506.51

#Number of Fisher Scoring iterations: 9


##Results of logistic regression on selected variables
#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-2.6357  -0.0721  -0.0530  -0.0530   3.6247  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  -6.56773    0.35852 -18.319  < 2e-16 ***
#  pelod2.gcs    0.51352    0.10660   4.817 1.46e-06 ***
#  pelod2.pup    0.52706    0.05651   9.327  < 2e-16 ***
#  pelod2.lac    0.61534    0.09730   6.324 2.54e-10 ***
#  pelod2.paco2  1.08656    0.14723   7.380 1.58e-13 ***
#  pelod2.wbc    0.88621    0.15731   5.634 1.76e-08 ***
#  ---
#  Signif. codes:  
#  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 1054.78  on 5117  degrees of freedom
#Residual deviance:  494.91  on 5112  degrees of freedom
#AIC: 506.91

#Number of Fisher Scoring iterations: 8


##Distributions of variables (raw values)

#GCS:
#10     11     12     13     14     15      3      4      5      6 
#37678  70534  11436   9081  24092 167148  70337   3618   3476  14682 
#7      8      9 
#12659  15067  43408 

#pup.left:
#              Brisk     Brisk, Brisk      Nonreactive         Sluggish      Unable to assess 
#364           372806                1            10358            41853             8394 

#pup.right:
#               Brisk     Brisk, Brisk      Nonreactive 
#365           372819                1            10525 
#Sluggish     Unable to assess 
#43189             6732 

#lac:
#0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0 10.2 10.3 10.9  1.1 11.0 11.6 
#1    6    8   18   28   54   51   48    1    1    1   41    1    1 
#1.2 12.2 12.6  1.3  1.4 14.6  1.5 15.0  1.6 16.5  1.7  1.8  1.9 19.0 
#47    1    1   34   49    1   40    1   32    1   38   27   25    1 
#2.0  2.1 21.0  2.2  2.3  2.4  2.5  2.6  2.7  2.8  2.9  3.0  3.1  3.2 
#16   17    1   20   17   15   12   11    7   15    7   13    4    6 
#3.3  3.4  3.5  3.6  3.7  3.8  3.9  4.0  4.1  4.2  4.3  4.4  4.5  4.6 
#8    8    2    6    5    2    3    4    2    8    6    7    3    3 
#4.7  4.8  4.9  5.0  5.1  5.2  5.4  5.5  5.6  5.7  5.8  5.9  6.0  6.1 
#5    5    1    6    1    2    3    1    3    1    2    1    1    2 
#6.2  6.3  6.4  6.5  6.6  6.7  6.8  6.9  7.0  7.1  7.2  7.4  7.5  8.2 
#2    2    1    1    1    1    1    1    3    1    2    1    1    1 
#8.5  9.2  9.4  9.5  9.6 
#1    2    1    1    1 

#map:
#0   10  100  101  102  103  104  105  106  107  108  109  110  111 
#1   11 2227   21 1874   18 1526   17 1293   15 1047   14  844    3 
#112  113  114  115  116  117  118  119   12  120  121  122  124  125 
#668    4  519    5  455    1  338    1    9  286    3  225  180    2 
#126  127  128  130  132  134  135  136  138   14  140  142  144  146 
#148    2  119  100   65   62    1   38   32   14   29   25   16   21 
#148  150  152  153  154  156  158   16  160  162  164  166  168  170 
#23    7   11    1   10   10    8    7    3    5    5    9    7    2 
#172  174  176  178   18  180  182  184  186  188  190  192  194  196 
#3    4    7    2    9    3    6    3    2    4    3    4    2    2 
#198   20  200  202  204  206  208  210  212  214  216  218   22  220 
#2   10    3    4    5    3    2    1    2    2    5    1   11    2 
#222  224  226  228  230  232  234  236  238   24  240  242  244  246 
#1    1    1    1    1    2    1    4    1    7    6    3    6    3 
#250  252  256  258   26  260  262  264  268  270  272  274  276  278 
#4    2    2    1   17    3    3    1    2    2    7    1    1    2 
#28  280  282  284  286  288   29  290  292  294  296  298   30  300 
#27    3    3    4    7    8    2    5    5    4    1    1   62    1 
#32   34   35   36   37   38   39   40   41   42   43   44   45   46 
#73  150    2  208    3  270    5  403    8  500    6  799    5 1057 
#47   48   49   50   51   52   53   54   55   56   57   58   59   60 
#12 1484   19 2033   18 2768   26 3435   62 4260   44 5222   48 5873 
#61   62   63   64   65   66   67   68   69   70   71   72   73   74 
#50 6617   65 7407   76 8114   79 8508   80 8690   79 8535   83 8348 
#75   76   77   78   79   80   81   82   83   84   85   86   87   88 
#73 8135   85 7577   75 7009   71 6669   65 6072   61 5566   56 4905 
#89   90   91   92   93   94   95   96   97   98   99 
#55 4357   41 4027   32 3335   29 3041   23 2612   27

#cr:
#<0.1      0.1    <0.15     0.15     0.16     0.17     0.18     0.19 
#142     2797     2432      120      314      237      289      203 
#0.2     0.20     0.21     0.22     0.23     0.24     0.25     0.26 
#9517      261      235      267      177      266      196      227 
#0.27     0.28     0.29      0.3     0.30     0.31     0.32     0.33 
#203      223      179     7580      206      167      207      163 
#0.34     0.35     0.36     0.37     0.38     0.39      0.4     0.40 
#183      162      185      167      170      151     5591      197 
#0.41     0.42     0.43     0.44     0.45     0.46     0.47     0.48 
#137      179      132      164      129      157      127      158 
#0.49      0.5     0.50     0.51     0.52     0.53     0.54     0.55 
#117     3973      150      130      141      118      154      117 
#0.56     0.57     0.58     0.59      0.6     0.60     0.61     0.62 
#152      132      133      104     2817      148       93      129 
#0.63     0.64     0.65     0.66     0.67     0.68     0.69      0.7 
#97      133      116       95       98      102       66     1946 
#0.70     0.71     0.72     0.73     0.74     0.75     0.76     0.77 
#115       93       92       63       88       66       85       67 
#0.78     0.79      0.8     0.80     0.81     0.82     0.83     0.84 
#78       70     1254       56       64       74       55       64 
#0.85     0.86     0.87     0.88     0.89      0.9     0.90     0.91 
#51       63       43       62       45      797       48       42 
#0.92     0.93     0.94     0.95     0.96     0.97     0.98     0.99 
#59       48       42       33       50       41       44       32 
#1.0     1.00     1.01    10.10     1.02     1.03     10.3     1.04 
#579       42       34        1       33       45        4       44 
#10.4     1.05     10.5     1.06     10.6     1.07     10.7     1.08 
#2       41        1       33        3       22        2       21 
#10.8     1.09     10.9      1.1     1.10     11.0     1.11     11.1 
#1       25        4      427       26        1       20        3 
#1.12     11.2     1.13     11.3     1.14     11.4     1.15     11.5 
#19        1       35        2       25        2       19        2 
#11.50     1.16     1.17     11.7     1.18     11.8     1.19     11.9 
#1       35       26        1       18        2       20        1 
#1.2     1.20     12.0    12.00     1.21     12.1    12.10     1.22 
#375       14        3        1       22        4        1       23 
#12.2     1.23     12.3     1.24     12.4     1.25     12.5     1.26 
#4       28        3       17        2       20        1       19 
#12.6     1.27     12.7     1.28     12.8     1.29     12.9      1.3 
#2       25        1       20        3       17        1      290 
#1.30     13.0     1.31     13.1     1.32     13.2     1.33     13.3 
#17        3       11        3       13        2        8        1 
#1.34     13.4     1.35     1.36     1.37     1.38     13.8     1.39 
#12        2       12       11        9       12        1       16 
#13.9      1.4     1.40     14.0     1.41     14.1     1.42     14.2 
#1      208       14        1        6        1        5        1 
#14.20     1.43     1.44     1.45     14.5    14.50     1.46     14.6 
#1       13       13        6        3        1       14        2 
#1.47     14.7    14.70     1.48     14.8     1.49     14.9      1.5 
#9        1        1        9        1        4        2      207 
#1.50     1.51     15.1     1.52     15.2     1.53    15.30     1.54 
#11       10        1       11        2       10        1        6 
#1.55     1.56     1.57     1.58     1.59      1.6     1.60     1.61 
#7        5        8       10        7      206        8        2 
#1.62     1.63     1.64     16.4     1.65     1.66     1.67     1.68 
#4        7        4        1       13        7        8        8 
#16.8     1.69     16.9      1.7     1.70    17.00     1.71     1.72 
#1        6        1      138        7        1        2        3 
#1.73     1.74     1.75     1.76     1.77     1.78     1.79      1.8 
#4        7        7        5        2        5        3      119 
#1.80    18.00     1.81     1.82     1.83     1.84     1.85     1.86 
#3        1        6        7        5        5        3        2 
#1.87     1.88     1.89      1.9     1.90     1.91     1.92     1.93 
#3        3        3      106        2       11        4        9 
#1.95     1.96     1.97     1.98     1.99      2.0     2.00     2.01 
#5        3        1        4        5       77        6        9 
#2.02     2.03     2.04     2.05     2.06     2.07     2.08     2.09 
#1        9        5        4        5        1        4        3 
#2.1     2.10     2.11     2.12     2.13     2.14     2.15     2.16 
#87        2        4        2        5        9        4        3 
#2.17     2.18     2.19      2.2     2.20     2.21     2.22     2.23 
#5        1        5       77        4        2        2        2 
#2.24     2.25     2.26     2.27     2.28     2.29      2.3     2.30 
#5        3        3        4        6        2       59        7 
#2.31     2.32     2.33     2.34     2.35     2.36     2.37     2.38 
#2        3        4        4        2        2        6        3 
#2.39      2.4     2.40     2.41     2.42     2.43     2.45     2.46 
#3       47        2        4        2        3        3        5 
#2.47     2.49      2.5     2.50     2.51     2.52     2.54     2.55 
#5        2       55        2        2        1        2        2 
#2.56     2.57     2.58     2.59      2.6     2.60     2.61     2.62 
#3        2        3        3       36        1        1        3 
#2.63     2.64     2.66     2.67     2.69      2.7     2.70     2.71 
#1        2        4        2        2       37        1        1 
#2.72     2.73     2.77     2.78     2.79      2.8     2.80     2.81 
#2        1        1        4        1       35        1        2 
#2.82     2.83     2.85     2.89      2.9     2.90     2.91     2.92 
#1        2        2        1       30        1        4        2 
#2.94     2.95     2.96     2.97     2.98     2.99      3.0     3.00 
#2        2        1        1        1        1       21        2 
#3.01     3.02     3.04     3.05     3.07     3.08     3.09      3.1 
#3        1        2        3        2        1        1       22 
#3.10     3.11     3.12     3.14     3.15     3.16     3.17     3.18 
#2        1        2        2        1        3        1        1 
#3.2     3.20     3.21     3.22     3.23     3.25     3.27     3.29 
#24        4        1        1        1        3        2        1 
#3.3     3.30     3.31     3.32     3.35     3.37     3.38      3.4 
#22        1        1        1        1        1        1       18 
#3.40     3.42     3.43      3.5     3.51     3.54     3.57      3.6 
#2        1        1       10        2        2        2       21 
#3.60     3.62     3.64     3.66     3.67      3.7     3.70     3.71 
#1        2        1        2        1       13        1        1 
#3.72     3.76      3.8     3.87      3.9     3.92     3.93     3.94 
#1        1       18        1       22        1        1        1 
#3.95     3.96     3.97     3.99      4.0     4.00     4.03     4.05 
#1        1        1        1       12        3        1        1 
#4.09      4.1     4.11     4.12     4.16     4.17     4.18     4.19 
#1       20        1        2        1        2        1        1 
#4.2     4.24     4.25     4.26     4.27     4.29      4.3     4.32 
#14        1        1        1        2        1       18        2 
#4.35     4.37     4.38     4.39      4.4     4.42     4.46     4.48 
#1        2        1        1       12        1        1        1 
#4.5     4.50     4.51     4.52     4.53     4.55     4.56     4.59 
#17        1        1        1        2        2        1        1 
#4.6     4.60     4.62     4.63     4.64     4.65     4.66      4.7 
#10        1        1        2        1        3        2       12 
#4.71     4.74     4.75     4.77     4.78     4.79      4.8     4.84 
#2        1        2        1        1        2        6        1 
#4.85     4.88      4.9     4.91     4.94     4.96      5.0 5.000000 
#1        1        8        1        2        1        9        9 
#5.01     5.03     5.05     5.07      5.1     5.15     5.18      5.2 
#1        1        1        2       11        3        1        8 
#5.22     5.23     5.26     5.27      5.3     5.31     5.34     5.35 
#1        1        3        1        5        1        1        1 
#5.38      5.4     5.45      5.5      5.6     5.64      5.7     5.72 
#1       14        1        5        5        1        6        1 
#5.73     5.78      5.8     5.83     5.89      5.9     5.91      6.0 
#1        2        4        1        1        5        1        7 
#6.05     6.07      6.1      6.2     6.21      6.3      6.4     6.42 
#1        1        4        3        1        4        3        1 
#6.5     6.56      6.6     6.62      6.7     6.70     6.78      6.8 
#8        1        3        1        6        1        1        5 
#6.9     6.95      7.0     7.07      7.1     7.13     7.18      7.2 
#3        1        4        1        3        1        1        8 
#7.39      7.4      7.5      7.6      7.7     7.79     7.83      7.9 
#2        2        3        2        3        2        1        1 
#7.95      8.0      8.1     8.11      8.2     8.24     8.28      8.3 
#1        3        5        1        5        1        1        4 
#8.34     8.37     8.38      8.5     8.55     8.61     8.69      8.7 
#1        1        1        1        1        1        1        2 
#8.78      8.8     8.87     8.89      8.9      9.0      9.1     9.15 
#1        2        1        1        3        5        4        1 
#9.17     9.21      9.3     9.30     9.37     9.38     9.45      9.5 
#2        1        2        1        1        1        1        2 
#9.52     9.57      9.6      9.7     9.73     9.76      9.8     9.85 
#1        1        2        2        1        1        2        1 
#9.9     9.97 
#1        1 

#pao2:
#100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 
#304 326 326 297 284 323 288 318 276 255 264 281 259 268 265 224 265 258 
#118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 
#249 269 243 244 236 238 232 222 229 218 196 206 192 235 197 216 206 176 
#136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 
#196 212 175 184 174 161 155 174 169 167 159 150 165 146 147 154 134 128 
#154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169  17 170 
#130 147 119 132 136 122 102 136 133 110 129  96 106 111 105 103   2 102 
#171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 
#87  85  96  89  85  79  77  76  70  91  75  72  58  68  68  67  54  67 
#189 190 191 192 193 194 195 196 197 198 199  20 200 201 202 203 204 205 
#82  63  65  60  51  66  64  48  62  55  66   1  65  45  49  47  47  56 
#206 207 208 209  21 210 211 212 213 214 215 216 217 218 219  22 220 221 
#40  40  52  45   2  49  50  47  35  44  45  48  36  42  47   1  32  38 
#222 223 224 225 226 227 228 229  23 230 231 232 233 234 235 236 237 238 
#44  37  37  31  32  35  30  30   4  47  34  32  36  33  37  33  46  22 
#239  24 240 241 242 243 244 245 246 247 248 249  25 250 251 252 253 254 
#28   5  24  31  28  32  28  30  26  33  36  29   5  34  30  32  25  25 
#255 256 257 258 259  26 260 261 262 263 264 265 266 267 268 269  27 270 
#26  32  27  25  31   2  26  21  26  26  24  22  25  29  24  17   5  15 
#271 272 273 274 275 276 277 278 279  28 280 281 282 283 284 285 286 287 
#26  16  19  19  14  13  18  23  23  10  15  18  17  15  16  12  14  15 
#288 289  29 290 291 292 293 294 295 296 297 298 299  30 300 301 302 303 
#11  18  11  18  16  17  15  19  17  11   9  10   9  16  19  14   8   7 
#304 305 306 307 308 309  31 310 311 312 313 314 315 316 317 318 319  32 
#14  14  13   6   8  10  18  10  11   9   7  10   8  16   6   2   7  27 
#320 321 322 323 324 325 326 327 328 329  33 330 331 332 333 334 335 336 
#9   4   8  13  12   6   6  10  11   7  33   3   8   8   7  11   8   4 
#337 338 339  34 340 341 342 343 344 345 346 347 348 349  35 350 351 352 
#13   4   2  34   3   5   5   3   7   3   7   3   3   5  43   5   4   4 
#353 354 355 356 357 358 359  36 360 361 362 363 364 365 366 367 368 369 
#6   9   8   5   6   7   7  71   5   4   2   2   5   4   4   5   7   2 
#37 370 371 372 373 374 375 376 377 378 379  38 380 381 382 383 384 385 
#58   5   4   5   4   8   7   6   5   6   9  80   3   2   1   4   5   4 
#386 387 388 389  39 390 391 392 393 394 395 396 397 398 399  40 402 403 
#1   4   4   1  88   1   5   6   5   3   5   6   1   1   4  80   7   2 
#404 405 406 407 408 409  41 410 411 412 413 414 415 416 417 418  42 420 
#3   6   6   5   1   2  86   5   2   6   2   1   2   1   2   2  84   3 
#421 422 423 424 425 426 427 429  43 431 432 433 434 435 436 437 438 439 
#1   3   1   2   1   1   1   1  82   1   3   1   4   1   4   3   2   2 
#44 440 442 443 444 445 446 448 449  45 450 451 452 453 454 455 456 457 
#98   3   4   2   2   1   2   2   2 100   7   6   3   2   1   2   2   5 
#458 459  46 461 462 463 464 465 466 467 468 469  47 470 471 472 473 474 
#1   4  94   1   6   4   4   2   2   6   1   2  88   3   4   1   1   2 
#476 477 478 479  48 480 481 482 483 484 485 486 488 489  49 490 491 492 
#2   4   3   3 103   3   1   2   1   2   2   2   1   2 106   1   4   1 
#493 494 495 496 497 498 499  50 500 501 502 506 507 508 509  51 510 511 
#3   2   3   2   5   2   3 126   2   1   1   5   3   2   1 142   1   3 
#512 513 514 518 519  52 520 521 522 524 525 526 527 528 529  53 530 531 
#2   3   4   1   2 150   1   2   3   1   1   2   2   2   1 159   2   2 
#532 533 534 535 536 537 538 539  54 540 541 542 543 544 545 546 548  55 
#2   1   5   2   1   1   2   2 164   1   3   3   3   2   4   2   1 197 
#551 552 553 554 556 557 558 559  56 563 564 565 566 567  57 570 571 572 
#2   1   1   2   1   1   2   3 181   3   2   1   1   2 213   3   1   1 
#573 576 578  58 582 584 587  59 598  60 603  61 613  62 625  63  64  65 
#1   3   1 225   2   1   1 241   1 275   2 291   2 285   1 332 318 340 
#659  66 669  67  68  69  70  71  72  73  74  75  76  77  78  79   8  80 
#1 364   1 346 393 398 464 410 428 424 456 411 467 472 463 425   1 454 
#81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98 
#437 392 440 446 406 417 419 370 404 392 377 347 366 363 386 350 343 312 
#99 
#298

#fio2:
#10   100    11    12    13    14    15    16    17    18    19    20 
#31 16825     2     7     1    36    49    32    26    24    17    90 
#21    22    23    24    25    26    27    28    29    30    31    32 
#30602   343  1021  1292  7480  1527  1348  2786  1019 22890  1222  1044 
#33    34    35    36    37    38    39    40    41    42    43    44 
#723   783 13127   675   387   745   596 19477   704   361   320   537 
#45    46    47    48    49    50    51    52    53    54    55    56 
#6005   325   173   428   481  8867   396   184   138   232  2542   274 
#57    58    59    60    61    62    63    64    65    66    67    68 
#285   315   218  4528   221   196   102   129  1356   139   191   110 
#69    70    71    72    73    74    75    76    77    78    79    80 
#156  1846    93   126    89   111   880   110    80    91   119  1403 
#81    82    83    84    85    86    87    88    89    90    91    92 
#93    76    63    75   537    70    53   100   114   614   109   153 
#93    94    95    96    97    98    99 
#221   284   418   312   394   480   374 

#paco2:
#10  100  101  102  103  104  105  106  107  108  109   11  110  111 
#8   13    7   11    9    3    5    7    8    5    8    7   11    3 
#112  113  114  116  117  118  119   12  120  121  122  123  124  125 
#7    5    9    6    6    4    3   11    2    6    2    1    2    3 
#126  127  128   13  130 >130   14   15   16   17   18   19   20   21 
#8    6    1    5    3   44    7    7    8   13   23   20   33   38 
#22   23   24   25   26   27   28   29   30   31   32   33   34   35 
#66   74   79  103  153  169  305  370  574  561  961  935 1180 1132 
#36   37   38   39   40   41   42   43   44   45   46   47   48   49 
#1515 1314 1725 1479 1901 1501 1827 1484 1772 1484 1683 1293 1423 1038 
#50   51   52   53   54   55   56   57   58   59   60   61   62   63 
#1162  860  872  703  738  531  579  407  475  321  305  271  298  178 
#64   65   66   67   68   69   70   71   72   73   74   75   76   77 
#219  149  173  130  148  110  114  102   90   75   77   57   66   43 
#78   79    8   80   81   82   83   84   85   86   87   88   89    9 
#40   40    1   39   36   33   24   27   21   22   19   24   20    3 
#90   91   92   93   94   95   96   97   98   99 
#20   14   18   13   18   11   13    6   11    8 

#vent:
#
#4 
#Avea 
#7830 
#Bird VIP Gold 
#15 
#High Frequency Oscillator 3100A 
#1617 
#High Frequency Oscillator 3100B 
#837 
#Home Care Device (specify in other) 
#19 
#Infant Nasal CPAP 
#7 
#LTV 
#262 
#Other: 
#  2 
#Other: 1150 
#1 
#Other: airlife 
#26 
#Other: air life 
#14 
#Other: Airlife 
#42 
#Other: Air life 
#1 
#Other: AirLife 
#8 
#Other: Air Life 
#13 
#Other: Air Life. 
#1 
#Other: AIRLIFE 
#5 
#Other: Airlife CPAP 
#3 
#Other: Air Life CPAP 
#1 
#Other: Airlife mask 
#1 
#Other: air life nasal cpap 
#1 
#Other: airlife nasal cpap 
#2 
#Other: Air Life Nasal Cpap 
#1 
#Other: airlife nasal mask 
#5 
#Other: airlife nasal prongs 
#2 
#Other: Ait Life CPAP 
#1 
#Other: AL 
#1 
#Other: Bagging and compressions 
#1 
#Other: Bagging endotube. 
#1 
#Other: Bag mask 
#1 
#Other: BAG VALVE MASK 
#1 
#Other: Cool aerosol mist 
#1 
#Other: cycled cpap 
#1 
#Other: hand ventilation 
#2 
#Other: HME 
#2 
#Other: HME trial 
#2 
#Other: home 
#1 
#Other: home trilogy 
#2 
#Other: Home Trilogy in use 
#1 
#Other: home unit 
#1 
#Other: LTV 
#1 
#Other: LTV 1150 
#3 
#Other: Maplesom Bag 
#1 
#Other: mapleson 
#1 
#Other: mapleson bag 
#2 
#Other: Mapleson bag 
#5 
#Other: MRI 
#1 
#Other: MRI MVP 10 
#1 
#Other: mri vent 
#1 
#Other: MRI ventilator 
#1 
#Other: MVP10 
#1 
#Other: nasal cannula, pt off cpap for the day 
#1 
#Other: niv/ps 
#1 
#Other: patient off bipap for the day on 1L NC 
#1 
#Other: pt placed on HME, sats 100%, tolerating well 
#1 
#Other: ram cannula 
#1 
#Other: ram  cannula 
#1 
#Other: S/T 
#2 
#Other: tpiece 
#1 
#Other: t piece 
#1 
#Other: TPiece on trach 
#2 
#Other: T Piece trial 
#2 
#Other: T-Piece trial 
#5 
#Other: T Piece Trial 
#1 
#Other: TPIECE trial 
#2 
#Other: Trilogy 
#1 
#Other: vent check not done at this time 
#1 
#Servo-i 
#70661 
#Trilogy 
#38902 
#Vision 
#93

#wbc:
#too many values to print, see histogram

#plate:
#1    10   100 10006  1007  1009   101  1011  1015  1018   102  1020 
#1    57   146     1     1     1   134     1     1     1   141     1 
#1021  1024  1026  1027  1028   103  1035  1039   104  1042  1044  1045 
#1     1     1     1     1   135     1     1   137     1     1     1 
#1046  1047  1049   105  1051  1055   106  1062  1065   107  1073  1077 
#1     1     1   145     2     2   131     1     1   143     1     2 
#1078   108  1080  1083   109    11   110  1101  1102  1108   111  1112 
#1   140     1     1   128    81   134     1     1     1   138     1 
#112  1120  1121  1126  1129   113  1135   114  1140  1141  1146   115 
#126     1     1     2     1   145     1   124     1     2     1   136 
#1152  1154   116  1160  1161  1162  1167  1168   117  1170  1173  1177 
#1     2   110     1     2     2     1     1   145     1     1     1 
#118   119  1193    12   120  1205  1206   121  1213   122  1227   123 
#120   121     1    63   133     2     1   122     1   124     1   124 
#1232  1233   124  1245   125  1253  1255  1259   126  1267   127  127. 
#1     1   146     1   104     1     1     1   123     1   117     1 
#1273   128   129    13   130   131  1316   132  1328   133   134   135 
#1   138   136    88   101   113     1   116     1   139   128   121 
#1354   136  1362   137   138  1383   139    14   140  1402   141  1414 
#2   107     1   130   131     1    94    87   127     1   112     1 
#142   143   144   145   146   147  1470  1475   148  1484   149    15 
#102   103   111   122   126   108     1     1   139     1   124   120 
#150   151   152  1522   153   154   155   156   157   158   159    16 
#110   106   121     1   103   106   104   113   102   114   104   103 
#160   161   162   163   164   165   166   167   168   169    17   170 
#88   111   100   101    96   128   103   110    98   113   114   101 
#171   172   173   174   175   176   177   178   179    18   180   181 
#96   105   102    84    97   118    95   111    82   133   108   107 
#182   183   184   185   186   187   188   189    19   190   191   192 
#105   104    93    94    94   107    90   117   129   104   105   109 
#193   194   195   196   197   198   199     2    20   200   201   202 
#119    92    98   118    85   106   103     5    96   100    93    95 
#203   204   205   206   207   208   209    21   210   211   212   213 
#120    92    88   102   124    92   102    80   113   110    87   107 
#214   215   216   217   218   219    22   220   221   222   223   224 
#110    97   101   118    94    84   105    89    89    92   106    95 
#225   226   227   228   229    23   230   231   232   233   234   235 
#111    92    96    78   112    94    85   109    94    82    96    90 
#236   237   238   239    24   240   241   242   243   244   245   246 
#82    93    95   111   134    93   105    89    86    98   104    99 
#247   248   249    25   250   251   252   253   254   255   256   257 
#100    95   108   114    92    71    80    87   101   100    79    78 
#258   259    26   260   261   262   263   264   265   266   267   268 
#108    87    84    97    84    71    92    94    75    83    90    89 
#269    27   270   271   272   273   274   275   276   277   278   279 
#91    91    68    65    89    90    76    79    77    72    68   103 
#28   280   281   282   283   284   285   286   287   288   289    29 
#81    86    70    86    78    61    84    76    90    65    88    88 
#290   291   292   293   294   295   296   297   298   299     3    30 
#74    79    82    75    81    63    72    84    74    64    38   106 
#300   301   302   303   304   305   306   307   308   309    31   310 
#70    71    73    63    68    77    59    69    68    57   102    55 
#311   312   313   314   315   316   317   318   319    32   320   321 
#61    77    77    70    56    55    58    54    84   122    63    56 
#322   323   324   325   326   327   328   329    33   330   331   332 
#38    65    70    64    63    70    67    61   108    55    66    59 
#333   334   335   336   337   338   339    34   340   341   342   343 
#59    57    48    46    51    54    52   100    51    65    55    46 
#344   345   346   347   348   349    35   350   351   352   353   354 
#39    45    64    54    63    52   115    62    43    48    52    49 
#355   356   357   358   359    36   360   361   362   363   364   365 
#50    42    44    52    45   115    50    41    48    41    47    38 
#366   367   368   369    37   370   371   372   373   374   375   376 
#36    43    45    48   111    40    34    40    43    39    45    42 
#377   378   379    38   380   381   382   383   384   385   386   387 
#57    45    43   109    53    45    60    37    34    42    41    34 
#388   389    39   390   391   392   393   394   395   396   397   398 
#42    38   121    39    28    40    34    36    40    45    32    32 
#399     4    40   400   401   402   403   404   405   406   407   408 
#35    41   131    36    30    36    32    40    39    34    29    36 
#409    41   410   411   412   413   414   415   416   417   418   419 
#37   118    32    33    31    37    28    30    33    33    29    34 
#42   420   421   422   423   424   425   426   427   428   429    43 
#96    28    24    30    28    33    35    23    29    30    28   117 
#430   431   432   433   434   435   436   437   438   439    44   440 
#28    32    31    24    30    23    27    30    24    21   131    24 
#441   442   443   444   445   446   447   448   449    45   450   451 
#19    27    23    26    34    24    27    25    25   103    21    22 
#452   453   454   455   456   457   458   459    46   460   461   462 
#22    42    21    20    12    28    18    29   115    19    19    24 
#463   464   465   466   467   468   469    47   470   471   472   473 
#21    16    17    22    17    24    20   119    26    21    14    19 
#474   475   476   477   478   479    48   480   481   482   483   484 
#17    22    23    16    13    18   106    23    24    15    20    15 
#485   486   487   488   489    49   490   491   492   493   494   495 
#22    10    21    18    19   118    19    16    15     6    16    14 
#496   497   498   499     5    50   500   501   502   503   504   505 
#17    16     6    16    39   128    21    16    19    19    12    21 
#506   507   508   509    51   510   511   512   513   514   515   516 
#10    27    19    15   138    12     8    18    15    12    21    18 
#517   518   519    52   520   521   522   523   524   525   526   527 
#10    16    14   123    16    15    11    23    10    10    14    11 
#528   529    53   530   531   532   533   534   535   536   537   538 
#16    17   116    12    13     9    15    17     7    11    13    18 
#539    54   540   541   542   543   544   545   546   547   548   549 
#10   103    11    15     9     9     6    14    21     5    10     8 
#55   550   551   552   553   554   555   556   557   558   559    56 
#133    10    10    11     7    13     6    10    10    10    13   121 
#560   561   562   563   564   565   566   567   568   569    57   570 
#12     9     6    18    14    11     9     6    11     8   121     6 
#571   572   573   574   575   576   577   578   579    58   580   581 
#14    13    13    14    19     7     5     7     7   117    15     5 
#582   583   584   585   586   587   588   589    59   590   591   592 
#6     8     8     9    11    10     6    10   132     9     4     5 
#593   594   595   596   597   598   599     6    60   600   601   602 
#9     9     3    14     8     6     7    32   128     7     9    11 
#603   604   605   606   607   608   609    61   610   611   612   613 
#6     5    10     3     5     5     4   112    11     4     7    10 
#614   615   616   617   618   619    62   620   621   622   623   624 
#6     7     6     4     3     8   134     1     3     5     2     3 
#625   626   627   628   629    63   630   631   632   633   634   635 
#3     1     3     6     2   133     7     6     5     4     4     1 
#636   637   638   639    64   640   641   642   643   644   645   646 
#7     2     5     4   125     2     8     5     4     1     3     5 
#647   648   649    65   650   651   652   653   654   655   656   657 
#4     2     2   120     4     5     4     2     8     5     6     3 
#658   659    66   660   661   662   663   664   665   666   667   668 
#4     6   116     4     6     2     2     4     3     3     3     4 
#669    67   670   671   672   673   674   675   676   677   678   679 
#4   120     2     2     2     5     6     4     4     3     1     1 
#68   680   681   682   683   684   685   686   687   688   689    69 
#128     4     4     3     1     4     2     2     3     3     2   122 
#690   691   692   693   694   695   697   698   699     7    70   700 
#3     1     5     3     4     3     7     3     1    35   138     1 
#702   703   704   705   706   707   708   709    71   711   712   713 
#2     2     7     6     7     1     4     3   156     2     2     4 
#714   715   716   717   718   719    72   720   721   722   723   724 
#4     4     2     2     2     4   114     2     4     3     2     2 
#725   726   727   728   729    73   730   731   732   733   734   735 
#4     2     1     2     3   165     1     1     3     5     4     1 
#736   737   738   739    74   741   743   744   745   746   747   748 
#2     1     1     2   131     4     1     3     3     4     2     2 
#749    75   750   751   752   754   755   756   757   758   759    76 
#2   133     2     3     1     2     1     4     3     2     5   128 
#762   763   765   766   767   768   769    77   771   773   775   776 
#1     1     4     2     1     1     1   136     2     2     3     2 
#777   778   779    78   780   781   782   783   784   785   786   787 
#1     1     1   112     2     3     3     4     4     1     2     6 
#788   789    79   790   792   793   794   795   796   797   798   799 
#3     2   145     2     2     3     1     1     4     2     1     1 
#8    80   800   801   802   804   805   806   807   808   809    81 
#67   134     3     6     2     3     1     1     2     1     1   138 
#811   812   813   816   817   818   819    82   822   823   824   825 
#1     2     1     1     1     2     2   134     2     1     1     2 
#826   828   829    83   830   831   833   835   836   837   838    84 
#1     1     2   159     1     4     1     2     1     1     1   123 
#840   842   843   844   845   848   849    85   853   854   856   857 
#2     1     2     1     2     1     2   137     2     2     1     1 
#858   859    86   860   861   864   867   869    87   870   871   873 
#1     3   145     1     1     2     1     1   139     1     1     1 
#875   876   878   879    88   880   881   886    89   890   892   894 
#1     3     1     3   115     1     1     1   135     1     3     1 
#895   897   898   899     9    90   900   902   905   907   908   909 
#1     1     1     1    53   139     1     1     3     2     1     1 
#91   910   911   913   915   917   918   919    92   920   921   924 
#143     1     3     1     1     2     1     2   121     2     1     1 
#926   927    93   933   935   936   937    94   940   942   944   948 
#1     1   151     1     2     1     1   129     1     1     1     1 
#949    95   952   958   959    96   961   962   964   965   966   968 
#1   140     1     1     2   128     1     2     2     1     1     2 
#969    97   970    98   980   981   982   983   985   986   988    99 
#1   130     1   125     1     1     2     1     1     1     2   129 
#991   992   994 
#2     1     2 

##Distributions of variables (PELOD-2 scores)

#pelod2.gcs:
#0    1    4 
#2903 1116 1099 

#pelod2.pup:
#0    5 
#4909  209

#pelod2.lac:
#0    1    4 
#4891  164   63

#pelod2.map:
#0    2    3    6 
#4275  596  178   69

#pelod2.cr:
#0    2 
#4138  980

#pelod2.carrico:
#0    2 
#4986  132

#pelod2.paco2:
#0    1    3 
#4795  276   47

#pelod2.vent:
#0    3 
#2909 2209

#pelod2.wbc:
#0    2 
#4882  236

#pelod2.plate:
#0    1    2 
#4030  602  486

##Distribution of outcomes

#(0 indicates survival and 1 indicates death)
#0    1 
#5009  109 


##To-do:
#1. Calculate AUROC and Hosmer-Lemeshow calibration of model for PICU data (done)
#2. Simple imputation of missing data (possibly unnecessary)
#3. Variable correlation (Spearman or Pearson)
#4. Try to reduce variables used (calculate AIC and/or BIC of each configuration) (at least partially done with all subsets regression)
#5. Find new cutoffs with CHAID method and decision trees/random forest
#6. Try generalized additive models, decision trees, support vector machine, etc.
#7. Calculate AUPRC of models


##Websites to look at:
#https://cran.r-project.org/web/packages/PredictABEL/PredictABEL.pdf
#https://www.r-bloggers.com/how-to-perform-a-logistic-regression-in-r/ 
#https://cran.r-project.org/web/packages/givitiR/vignettes/givitiR.html 
#https://www.r-bloggers.com/calculating-auc-the-area-under-a-roc-curve/ 
#https://www.r-bloggers.com/illustrated-guide-to-roc-and-auc/ 
#http://thestatsgeek.com/2014/05/05/area-under-the-roc-curve-assessing-discrimination-in-logistic-regression/
#http://thestatsgeek.com/2014/02/16/the-hosmer-lemeshow-goodness-of-fit-test-for-logistic-regression/   
#http://myrcodes.blogspot.com/2013/12/area-under-curve-auc-proc-package.html 
#https://rpubs.com/Wangzf/pROC 
#https://cran.r-project.org/web/packages/pROC/pROC.pdf 
#https://r-forge.r-project.org/R/?group_id=343 
#https://rstudio-pubs-static.s3.amazonaws.com/2897_9220b21cfc0c43a396ff9abf122bb351.html
#https://cran.r-project.org/web/packages/glmulti/glmulti.pdf
#https://cran.r-project.org/web/packages/PRROC/PRROC.pdf 
#http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf 
#http://www.milanor.net/blog/cross-validation-for-predictive-analytics-using-r/ 
