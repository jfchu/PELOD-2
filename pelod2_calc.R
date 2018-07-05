AgeMonths <- function(age.days) {
  return (age.days / 365 * 12)
}

#Global vector of age cutoffs measured in months
ages <<- c(144, 60, 24, 12, 1)

#Logical values for score functions -> if clinical value should be less than a certain threshold

#Glasgow Coma Score
#True
GComa <- function(score) {
  cutoffs <- list(c(4, 4), c(10, 1))
  for (i in seq_len(length(cutoffs))) {
    if (score <= cutoffs[[i]][1])
      return (cutoffs[[i]][2])
  }
  return (0)
}

#Pupil reactivity should be entered as follows: 1 if true or 0 if false
#True
Pupillary <- function(left.reactive, right.reactive) {
  cutoffs <- list(c(0, 5))
  if (left.reactive <= cutoffs[[1]][1] && right.reactive <= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0)
}

#Molarity measured in mmol/L
#False
Lactatemia <- function(lac.molarity, lac.wholeblood.molarity) {
  cutoffs <- list(c(11.0, 4), c(5.0, 1))
  FindCutoffs <- function(X, molarity) {
    if (molarity >= cutoffs[[X]][1])
      return (cutoffs[[X]][2])
  }
  if (is.na(lac.molarity) || is.na(lac.wholeblood.molarity))
    return (NA)
  else if (lac.molarity >= lac.wholeblood.molarity) 
    sapply(seq_len(length(cutoffs)), FindCutoffs, molarity = lac.molarity)
  else
    sapply(seq_len(length(cutoffs)), FindCutoffs, molarity = lac.wholeblood.molarity)
  return (0)
}

#Age measured in months
#Mean arterial pressure measure in mm Hg
#True
MAP <- function(age, pressure) {
  cutoffs <- list(list(c(37, 6), c(51, 3), c(66, 2)), list(c(35, 6), c(48, 3), c(64, 2)), list(c(31, 6), c(44, 3), c(61, 2)), list(c(30, 6), c(43, 3), c(59, 2)), list(c(24, 6), c(38, 3), c(54, 2)), list(c(16, 6), c(30, 3), c(45, 2)))
  age.index <- 6
  FindAgeCutoffs <- function(X, age) {
    if (age >= ages[X]) {
      return (X)
    }
  }
  age.index <- max(sapply(seq_len(length(ages)), FindAgeCutoffs, age = age))
  FindCutoffs <- function(X, pressure) {
    if (pressure <= cutoffs[[age.index]][[X]][1])
      return (cutoffs[[age.index]][[X]][2])
  }
  sapply(seq_len(length(cutoffs)), FindCutoffs, pressure = pressure)
  return (0)
}

#Age measured in months
#Molarity measured in μmol/L
#False
Creatinine <- function(age, cr.molarity) {
  cutoffs <- list(list(c(93, 2)), list(c(59, 2)), list(c(51, 2)), list(c(35, 2)), list(c(23, 2)), list(c(70, 2)))  
  age.index <- 6
  for (i in seq_len(length(ages))) {
    if (age >= ages[i]) {
      age.index <- i
      break
    }
  }
  for (i in seq_len(length(cutoffs[[1]]))) {
    if (cr.molarity <= cutoffs[[age.index]][[i]][1])
      return (cutoffs[[age.index]][[i]][2])
  }
  return (0)
}

#Carrico index = PaO2(measured in mm Hg)/FIO2
#PaO2 -> partial pressure arterial oxygen
#FIO2 -> fraction of inspired oxygen
#True
Carrico <- function(pao2, fio2) {
  cutoffs <- list(c(60, 2))
  if (pao2 / fio2 <= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0) 
}

#PaCO2(measured in mm HG) -> partial pressure arterial carbon dioxide
#True
PaCO2 <- function(partial.pressure) {
  cutoffs <- list(c(95, 3), c(59, 1))
  for (i in seq_len(length(cutoffs))) {
    if (partial.pressure >= cutoffs[[i]][1]) 
      return (cutoffs[[i]][2])
  }
  return (0)
}

#Presence of invasive ventilation should be entered as follows: 1 if true or 0 if false
#False
Ventilation <- function(is.invasive) {
  cutoffs <- list(c(1, 3))
  if (is.invasive >= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0)
}

#wbc.count * 10^9/L
#True
WBC <- function(wbc.count) {
  cutoffs <- list(c(2, 2))
  if (wbc.count <= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0)
}

#plate.count * 10^9/L
#True
Platelet <- function(plate.count) {
  cutoffs <- list(c(76, 2), c(141, 1))
  for (i in seq_len(length(cutoffs))) {
    if (plate.count <= cutoffs[[i]][1])
      return (cutoffs[[i]][2])
  }
  return (0)
}

#11 parameters
#Age measured in months
HolisticPELOD2 <- function(age, gcoma, pup, lac, pressure, cr, carrico, co2.pressure, invasive.vent, wbc.count, plate.count) {
  return (NeuroPELOD2(gcoma, pup) + CardioPELOD2(age, lac, pressure) + RenalPELOD2(age, cr) + RespPELOD2(carrico, co2.pressure, invasive.vent) + HemaPELOD2(wbc.count, plate.count)) 
}

NeuroPELOD2 <- function(gcoma.score, pup.reactive) {
  return (GComa(gcoma.score) + Pupillary(pup.reactive))
}

CardioPELOD2 <- function(age.months, lac.molar, map) {
  return (Lactatemia(lac.molar) + MAP(age.months, map))
}

RenalPELOD2 <- function(age.months, cr.molar) {
  return (Creatinine(age.months, cr.molar))
}

RespPELOD2 <- function(carrico.index, paco2, vent) {
  return (Carrico(carrico.index) + PaCO2(paco2) + Ventilation(vent))
}

HemaPELOD2 <- function(wbc, plate) {
  return (WBC(wbc) + Platelet(plate))
}

#Extracts patient data from data frame (layout below) and calculates PELOD-2 score
#Precondition: frame contains at least 1 row
#HolisticPELOD2Frame <- function(frame) {
#  frame <- cbind(frame, AGEMONTHS = as.numeric(pelod2.data$AGEDAYS) / 365 * 12)
#  pelod2.score <- c(HolisticPELOD2(frame$AGEMONTHS[1], frame[1, 3], frame[1, 4], frame[1, 5], frame[1, 6], frame[1, 7], frame[1, 8], frame[1, 9], frame[1, 10], frame[1, 11], frame[1, 12]))
#  for (i in seq_len(nrow(frame))[-1]) {
#    pelod2.score <- c(pelod2.score, c(HolisticPELOD2(frame[i, 13], frame[i, 3], frame[i, 4], frame[i, 5], frame[i, 6], frame[i, 7], frame[i, 8], frame[i, 9], frame[i, 10], frame[i, 11], frame[i, 12])))
#  }
#  return (cbind(frame, pelod2.score))
#}

#Precondition: frame contains at least 1 row
#NeuroPELOD2Frame <- function(frame) {
#  pelod2.neuroscore <- c(NeuroPELOD2(frame[1, 3], frame[1, 4]))
#  for (i in seq_len(nrow(frame))[-1]) {
#    pelod2.neuroscore <- c(pelod2.neuroscore, c(NeuroPELOD2(frame[i, 3], frame[i, 4])))
#  }
#  return (cbind(frame, pelod2.neuroscore))
#}

#Precondition: frame contains at least 1 row
#CardioPELOD2Frame <- function(frame) {
#  frame <- cbind(frame, pelod2.months = as.numeric(pelod2.data$pelod2.age) / 365 * 12)
#  pelod2.cardioscore <- c(CardioPELOD2(frame[1, 13], frame[1, 5], frame[1, 6]))
#  for (i in seq_len(nrow(frame))[-1]) {
#    pelod2.cardioscore <- c(pelod2.cardioscore, c(CardioPELOD2(frame[i, 13], frame[i, 5], frame[i, 6])))
#  }
#  return (cbind(frame, pelod2.cardioscore))
#}

#Precondition: frame contains at least 1 row
#RenalPELOD2Frame <- function(frame) {
#  frame <- cbind(frame, pelod2.months = as.numeric(pelod2.data$pelod2.age) / 365 * 12)
#  pelod2.renalscore <- c(RenalPELOD2(as.numeric((frame$pelod2.months)[1]), frame[1, 7]))
#  for (i in seq_len(nrow(frame))[-1]) {
#    pelod2.renalscore <- c(pelod2.renalscore, c(RenalPELOD2(frame[i, 13], frame[i, 5], frame[i, 6])))
#  }
#  return (cbind(frame, pelod2.renalscore))
#}

#Precondition: frame contains at least 1 row
#RespPELOD2Frame <- function(frame) {
#  pelod2.respscore <- c(RespPELOD2(frame[1, 8], frame[1, 9], frame[1, 10]))
#  for (i in seq_len(nrow(frame))[-1]) {
#    pelod2.respscore <- c(pelod2.respscore, c(RespPELOD2(frame[i, 8], frame[i, 9], frame[i, 10])))
#  }
#  return (cbind(frame, pelod2.respscore))
#}

#Precondition: frame contains at least 1 row
#HemaPELOD2Frame <- function(frame) {
#  pelod2.hemascore <- c(HemaPELOD2(frame[1, 11], frame[1, 12]))
#  for (i in seq_len(nrow(frame))[-1]) {
#    pelod2.hemascore <- c(pelod2.hemascore, c(HemaPELOD2(frame[i, 11], frame[i, 12])))
#  }
#  return (cbind(frame, pelod2.hemascore))
#}

ProbMortality <- function(pelod2) {
  return (1 / (1 + exp(6.61 - 0.47 * pelod2)))
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
  return (min(unlist(sapply(vals, FilterNonreactive))))
}

FindExtremeValueVent <- function(id, frame) {
  vals <- frame[which(frame$encid == id), "Clinical.Event.Result"]
  FilterCharacter <- function(X) {if (is.character(X)) return (1) else return (0)}
  return (max(sapply(vals, FilterCharacter)))
}

PELOD2Scores <- function(frame.list) {
  pelod2.frame <- data.frame(
    pelod2.id = read.xlsx("PELOD2_calculator_deid_unencrypt.xlsx")$encid,
    pelod2.age = AgeMonths(read.xlsx("admit1picu_deid_unencrypt.xlsx")$AGEDAYS))
    pelod2.frame$pelod2.gcs = GComa(sapply(pelod2.frame$pelod2.id, FindExtremeValueMin, frame = frame.list[["GCS"]]))
    pelod2.frame$pelod2.pup = Pupillary(sapply(pelod2.frame$pelod2.id, FindExtremeValuePupil, frame = frame.list[["pup.left"]]), sapply(pelod2.frame$pelod2.id, FindExtremeValuePupil, frame = frame.list[["pup.right"]]))
    pelod2.frame$pelod2.lac = Lactatemia(sapply(pelod2.frame$pelod2.id, FindExtremeValueMax, frame = frame.list[["lac"]]), sapply(pelod2.frame$pelod2.id, FindExtremeValueMax, frame = frame.list[["lac.wholeblood"]]))
    pelod2.frame$pelod2.map = MAP(pelod2.frame$pelod2.age, sapply(pelod2.frame$pelod2.id, FindExtremeValueMin, frame = frame.list[["map"]])) 
    pelod2.frame$pelod2.cr = Creatinine(pelod2.frame$pelod2.age, sapply(pelod2.frame$pelod2.id, FindExtremeValueMin, frame = frame.list[["cr"]]))
    pelod2.frame$pelod2.carrico = Carrico(sapply(pelod2.frame$pelod2.id, FindExtremeValueMax, frame = frame.list[["pao2"]]), sapply(pelod2.frame$pelod2.id, FindExtremeValueMin, frame = frame.list[["fio2"]]))
    pelod2.frame$pelod2.paco2 = PaCO2(sapply(pelod2.frame$pelod2.id, FindExtremeValueMax, frame = frame.list[["paco2"]]))
    pelod2.frame$pelod2.vent = Ventilation(sapply(pelod2.frame$pelod2.id, FindExtremeValueVent, frame = frame.list[["vent"]]))
    pelod2.frame$pelod2.wbc = WBC(sapply(pelod2.frame$pelod2.id, FindExtremeValueMin, frame = frame.list[["wbc"]]))
    pelod2.frame$pelod2.plate = Platelet(sapply(pelod2.frame$pelod2.id, FindExtremeValueMin, frame = frame.list[["plate"]]))
    pelod2.frame$pelod2.scores = pelod2.frame$pelod2.gcs + pelod2.frame$pelod2.pup + pelod2.frame$pelod2.lac + pelod2.frame$pelod2.map + pelod2.frame$pelod2.cr + pelod2.frame$pelod2.carrico + pelod2.frame$pelod2.paco2 + pelod2.frame$pelod2.vent + pelod2.frame$pelod2.wbc + pelod2.frame$pelod2.plate
}

PELOD2Ages <- function(frame) {
  return (frame$AGEYRS * 365.25 + frame$AGEDAYS)
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

PELOD2Scores(pelod2.datalist)
#for (i in 1:13) 
#  pelod2.datalist[[i]][["Clinical.Event.Performed.Date/Time"]] <- convertToDateTime(pelod2.datalist[[i]][["Clinical.Event.Performed.Date/Time"]])

#pelod2.datalist[[1]][["Clinical.Event.Performed.Date/Time"]] <- convertToDateTime(pelod2.datalist[[1]][["Clinical.Event.Performed.Date/Time"]])
#print (pelod2.datalist[[1]]["Clinical.Event.Performed.Date/Time"])

#To-do:
#finish implement changes detailed in pelod2_calc_revised
#validate data in pelod2 calculator file

PELOD2LacScores <- function(frame.list) {
  pelod2.frame <- data.frame(
    pelod2.id = read.xlsx("PELOD2_calculator_deid_unencrypt.xlsx")$encid,
    pelod2.age = AgeMonths(PELOD2Ages(read.xlsx("admit1picu_deid_unencrypt.xlsx"))))
  pelod2.frame$pelod2.lac = Lactatemia(sapply(pelod2.frame$pelod2.id, FindExtremeValueMax, frame = frame.list[["lac"]]), sapply(pelod2.frame$pelod2.id, FindExtremeValueMax, frame = frame.list[["lac.wholeblood"]]))
}
