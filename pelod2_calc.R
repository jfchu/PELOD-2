AgeMonths <- function(age.days) {
  return (age.days / 365 * 12)
}

#Global vector of age cutoffs measured in months
ages <<- c(144, 60, 24, 12, 1)

#Logical values for score functions -> if clinical value should be less than a certain threshold

#Glasgow Coma Score
#True
GComa <- function(id) {
  cutoffs <- list(c(4, 4), c(10, 1))
  min.val <- FindExtremeValueMin(id, pelod2.datalist[["GCS"]])
  FindCutoffs <- function(X, score) {
    if (is.na(score))
      return (NA)
    else if (score <= cutoffs[[1]][1])
      return (cutoffs[[X]][2])
    return (0)
  }
  sapply(seq_len(length(cutoffs)), FindCutoffs, score = min.val)
}

#Pupil reactivity should be entered as follows: 1 if true or 0 if false
#True
Pupillary <- function(id) {
  cutoffs <- list(c(0, 5))
  min.val1 <- FindExtremeValueMin(id, pelod2.datalist[["pup.left"]])
  min.val2 <- FindExtremeValueMin(id, pelod2.datalist[["pup.right"]])
  if (is.na(min.val1) || is.na(min.val2))
    return (NA)
  else if (min.val1 <= cutoffs[[1]][1] && min.val2 <= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0)
}

#Molarity measured in mmol/L
#False
#!is.numeric(unlist(FindValues(id, pelod2.datalist[["lac.wholeblood"]]), recursive = F))
Lactatemia <- function(id) {
  cutoffs <- list(c(11.0, 4), c(5.0, 1))
  max.val1 <- FindExtremeValueMax(id, pelod2.datalist[["lac"]])
  if (!is.numeric(FindValues(id, pelod2.datalist[["lac.wholeblood"]]))) 
    return (4)
  else
    max.val2 <- FindExtremeValueMax(id, pelod2.datalist[["lac.wholeblood"]])
  FindCutoffs <- function(X, molarity) {
    if (is.na(molarity))
      return (NA)
    else if (molarity >= cutoffs[[X]][1])
      return (cutoffs[[X]][2])
    return (0)
  }
  if (is.na(max.val1) || is.na(max.val2))
    return (NA)
  else if (max.val1 >= max.val2) 
    sapply(seq_len(length(cutoffs)), FindCutoffs, molarity = max.val1)
  else
    sapply(seq_len(length(cutoffs)), FindCutoffs, molarity = max.val2)
}
Lactatemia(1093)
#Age measured in months
#Mean arterial pressure measure in mm Hg
#True
MAP <- function(id, age) {
  cutoffs <- list(list(c(37, 6), c(51, 3), c(66, 2)), list(c(35, 6), c(48, 3), c(64, 2)), list(c(31, 6), c(44, 3), c(61, 2)), list(c(30, 6), c(43, 3), c(59, 2)), list(c(24, 6), c(38, 3), c(54, 2)), list(c(16, 6), c(30, 3), c(45, 2)))
  min.val <- FindExtremeValueMin(id, pelod2.datalist[["map"]])
  age.index <- 6
  FindAgeCutoffs <- function(X, age) {
    if (age >= ages[X]) 
      return (X)
  }
  age.index <- max(unlist(sapply(seq_len(length(ages)), FindAgeCutoffs, age = age), recursive = F))
  FindCutoffs <- function(X, pressure) {
    if (is.na(pressure))
      return (NA)
    else if (pressure <= cutoffs[[age.index]][[X]][1])
      print(X)
      return (cutoffs[[age.index]][[X]][2])
    return (0)
  }
  sapply(seq_len(length(cutoffs[[1]])), FindCutoffs, pressure = min.val)
}

#Age measured in months
#Molarity measured in μmol/L
#False
Creatinine <- function(id, age) {
  cutoffs <- list(list(c(93, 2)), list(c(59, 2)), list(c(51, 2)), list(c(35, 2)), list(c(23, 2)), list(c(70, 2)))
  ReplaceNaN <- function()
  if (!is.numeric(FindValues(id, pelod2.datalist[["cr"]])))
    
  max.val <- FindExtremeValueMax(id, pelod2.datalist[["cr"]])
  age.index <- 6
  FindAgeCutoffs <- function(X, age) {
    if (age > ages[X])
      return (X)
  }
  age.index <- max(unlist(sapply(seq_len(length(ages)), FindAgeCutoffs, age = age), recursive = F))
  FindCutoffs <- function(X, cr.molarity) {
    if (is.na(cr.molarity))
      return (NA)
    else if (cr.molarity >= cutoffs[[age.index]][[X]][1])
      return (cutoffs[[age.index]][[X]][2])
    return (0)
  }
  sapply(seq_len(length(cutoffs[[1]])), FindCutoffs, cr.molarity = max.val)
}

#Carrico index = PaO2(measured in mm Hg)/FIO2
#PaO2 -> partial pressure arterial oxygen
#FIO2 -> fraction of inspired oxygen
#True
Carrico <- function(id) {
  cutoffs <- list(c(60, 2))
  min.val <- FindExtremeValueMin(id, pelod2.datalist[["pao2"]])
  max.val <- FindExtremeValueMax(id, pelod2.datalist[["fio2"]])
  if (is.na(min.val) || is.na(max.val))
    return (NA)
  else if (min.val / max.val <= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0) 
}

#PaCO2(measured in mm HG) -> partial pressure arterial carbon dioxide
#True
PaCO2 <- function(id) {
  cutoffs <- list(c(95, 3), c(59, 1))
  min.val <- FindExtremeValueMin(id, pelod2.datalist[["paco2"]])
  FindCutoffs <- function(X, partial.pressure) {
    if (is.na(partial.pressure))
      return (NA)
    else if (partial.pressure <= cutoffs[[X]][1])
      return (cutoffs[[X]][2])
    return (0)
  }
  sapply(seq_len(length(cutoffs)), FindCutoffs, partial.pressure = min.val)
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
WBC <- function(id) {
  cutoffs <- list(c(2, 2))
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
  sapply(seq_len(length(cutoffs)), FindCutoffs, plate.count = min.val)
}

ProbMortality <- function(pelod2) {
  return (1 / (1 + exp(6.61 - 0.47 * pelod2)))
}

FindValues <- function(id, frame) {
  vals <- frame[which(frame$encid == id), "Clinical.Event.Result"]
  return (vals)
}
FindExtremeValueMax <- function(id, frame) {
  vals <- frame[which(frame$encid == id), "Clinical.Event.Result"]
  return (max(vals))
}

FindExtremeValueMin <- function(id, frame) {
  vals <- frame[which(frame$encid == id), "Clinical.Event.Result"]
  return (min(vals))
}

FindExtremeValuePupil <- function(id, frame) {
  vals <- frame[which(frame$encid == id), "Clinical.Event.Result"]
  FilterNonreactive <- function(X) {if (X == "Nonreactive") return (0) else return (1)}
  return (min(unlist(sapply(vals, FilterNonreactive), recursive = F)))
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
    pelod2.frame$pelod2.lac = sapply(pelod2.frame$pelod2.id, Pupillary)
    pelod2.frame$pelod2.map = mapply(MAP, pelod2.frame$pelod2.id, pelod2.frame$pelod2.age)
    pelod2.frame$pelod2.cr = mapply(Creatinine, pelod2.frame$pelod2.id, pelod2.frame$pelod2.age)
    pelod2.frame$pelod2.carrico = sapply(pelod2.frame$pelod2.id, Carrico) 
    pelod2.frame$pelod2.paco2 = sapply(pelod2.frame$pelod2.id, PaCO2)
    pelod2.frame$pelod2.vent = sapply(pelod2.frame$pelod2.id, Ventilation)
    pelod2.frame$pelod2.wbc = sapply(pelod2.frame$pelod2.id, WBC)
    pelod2.frame$pelod2.plate = sapply(pelod2.frame$pelod2.id, Platelet)
    pelod2.frame$pelod2.neuroscores = pelod2.frame$pelod2.gcs + pelod2.frame$pelod2.pup
    pelod2.frame$pelod2.cardioscores = pelod2.frame$pelod2.lac + pelod2.frame$pelod2.map
    pelod2.frame$pelod2.renalscores = pelod2.frame$pelod2.cr
    pelod2.frame$pelod2.respscores = pelod2.frame$pelod2.carrico + pelod2.frame$pelod2.paco2 + pelod2.frame$pelod2.vent
    pelod2.frame$pelod2.hemascores = pelod2.frame$pelod2.wbc + pelod2.frame$pelod2.plate
    pelod2.frame$pelod2.scores = pelod2.frame$pelod2.neuroscores + pelod2.frame$pelod2.cardioscores + pelod2.frame$pelod2.renalscores + pelod2.frame$pelod2.respscores + pelod2.frame$pelod2.hemascores
#consider whether to calculate total score by adding neuro, cardio, ... or adding gcs, pup, ...
#see whether NA values will influence above decision
#    pelod2.frame$pelod2.scores = pelod2.frame$pelod2.gcs + pelod2.frame$pelod2.pup + pelod2.frame$pelod2.lac + pelod2.frame$pelod2.map + pelod2.frame$pelod2.cr + pelod2.frame$pelod2.carrico + pelod2.frame$pelod2.paco2 + pelod2.frame$pelod2.vent + pelod2.frame$pelod2.wbc + pelod2.frame$pelod2.plate
    return (pelod2.frame)
}

RemoveNA <- function(df) {
  return(df[complete.cases(df), ])
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

View(PELOD2Scores(pelod2.datalist))

#for (i in 1:13) 
#  pelod2.datalist[[i]][["Clinical.Event.Performed.Date/Time"]] <- convertToDateTime(pelod2.datalist[[i]][["Clinical.Event.Performed.Date/Time"]])

#pelod2.datalist[[1]][["Clinical.Event.Performed.Date/Time"]] <- convertToDateTime(pelod2.datalist[[1]][["Clinical.Event.Performed.Date/Time"]])
#print (pelod2.datalist[[1]]["Clinical.Event.Performed.Date/Time"])

#Issues to consider:
#1. Handling non-numeric data
#Non-numeric values in Lactate Whole Blood(>), Creatine(<), PaCO2(>), WBC(<) 
#2. Unlisting list containing character and numeric data (coercing numeric data to character)
#3. Creatine vs. creatinine data value differences in different files
