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
#Fix score calculation; eg. encid 7 has incorrect score
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
#!is.numeric(unlist(FindValues(id, pelod2.datalist[["lac.wholeblood"]]), recursive = F))
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
#max never selects the specific non-numeric data in this sheet but a workaround would be ideal
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

ProbMortality <- function(pelod2) {
  return (1 / (1 + exp(6.61 - 0.47 * pelod2)))
}

ProbMortalityNew <- function(pelod2) {
  return (1 / (1 + exp(-6.76204 + 0.3904402 * pelod2$pelod2.gcs + 0.5149909 * pelod2$pelod2.pup + 0.1381793 * pelod2$pelod2.map - 0.06871416 * pelod2$pelod2.cr + 0.2626874 * pelod2$pelod2.carrico + 0.8777295 * pelod2$pelod2.paco2 + 0.1311851 * pelod2$pelod2.vent + 0.7447805 * pelod2$pelod2.wbc + 0.2743563 * pelod2$pelod2.plate)))
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
#    pelod2.frame$pelod2.neuroscores = pelod2.frame$pelod2.gcs + pelod2.frame$pelod2.pup
#    pelod2.frame$pelod2.cardioscores = pelod2.frame$pelod2.lac + pelod2.frame$pelod2.map
#    pelod2.frame$pelod2.renalscores = pelod2.frame$pelod2.cr
#    pelod2.frame$pelod2.respscores = pelod2.frame$pelod2.carrico + pelod2.frame$pelod2.paco2 + pelod2.frame$pelod2.vent
#    pelod2.frame$pelod2.hemascores = pelod2.frame$pelod2.wbc + pelod2.frame$pelod2.plate
#    pelod2.frame$pelod2.scores = pelod2.frame$pelod2.neuroscores + pelod2.frame$pelod2.cardioscores + pelod2.frame$pelod2.renalscores + pelod2.frame$pelod2.respscores + pelod2.frame$pelod2.hemascores
#consider whether to calculate total score by adding neuro, cardio, ... or adding gcs, pup, ...
#see whether NA values will influence above decision
#    pelod2.frame$pelod2.scores = pelod2.frame$pelod2.gcs + pelod2.frame$pelod2.pup + pelod2.frame$pelod2.lac + pelod2.frame$pelod2.map + pelod2.frame$pelod2.cr + pelod2.frame$pelod2.carrico + pelod2.frame$pelod2.paco2 + pelod2.frame$pelod2.vent + pelod2.frame$pelod2.wbc + pelod2.frame$pelod2.plate
    pelod2.frame$pred = ProbMortalityNew(pelod2.frame)
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
LOOCV <- function(data) {
acc <- NULL
for(i in 1:nrow(data))
{
  # Train-test splitting
  # 499 samples -> fitting
  # 1 sample -> testing
  train <- data[-i,]
  test <- data[i,]
  
  # Fitting
  model2 <- glm(deceased~.,family=binomial,data=train)
  
  # Predict results
  results_prob <- predict(model2,subset(test,select=c(2:9)),type='response')
  
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
###################################################################

##To-do:
#1. Calculate AUROC and Hosmer-Lemeshow calibration of model for PICU data (done)
#2. Simple imputation of missing data (possibly unnecessary)
#3. Variable correlation (Spearman or Pearson)
#4. Try to reduce variables used (calculate AIC and/or BIC of each configuration)
#5. Find new cutoffs with CHAID method and decision trees/random forest
#6. Try generalized additive models, decision trees, support vector machine, etc.


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
