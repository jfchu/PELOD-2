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
Pupillary <- function(reactive) {
  cutoffs <- list(c(0, 5))
  if (reactive <= cutoffs[[1]][1])
    return (cutoffs[[1]][2])
  return (0)
}

#Molarity measured in mmol/L
#False
Lactatemia <- function(lac.molarity) {
  cutoffs <- list(c(11.0, 4), c(5.0, 1))
  for (i in seq_len(length(cutoffs))) {
    if (lac.molarity >= cutoffs[[i]][1])
      return (cutoffs[[i]][2])
  }
  return (0)
}

#Age measured in months
#Mean arterial pressure measure in mm Hg
#True
MAP <- function(age, pressure) {
  cutoffs <- list(list(c(37, 6), c(51, 3), c(66, 2)), list(c(35, 6), c(48, 3), c(64, 2)), list(c(31, 6), c(44, 3), c(61, 2)), list(c(30, 6), c(43, 3), c(59, 2)), list(c(24, 6), c(38, 3), c(54, 2)), list(c(16, 6), c(30, 3), c(45, 2)))
  age.index <- 6
  for (i in seq_len(length(ages))) {
    if (age >= ages[i]) {
      age.index <- i
      break
    }
  }
  for (i in seq_len(length(cutoffs[[1]]))) {
    if (pressure <= cutoffs[[age.index]][[i]][1])
      return (cutoffs[[age.index]][[i]][2])
  }
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
Carrico <- function(index) {
  cutoffs <- list(c(60, 2))
  if (index <= cutoffs[[1]][1])
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
HolisticPELOD2Frame <- function(frame) {
  frame <- cbind(frame, AGEMONTHS = as.numeric(pelod2.data$AGEDAYS) / 365 * 12)
  pelod2.score <- c(HolisticPELOD2(frame$AGEMONTHS[1], frame[1, 3], frame[1, 4], frame[1, 5], frame[1, 6], frame[1, 7], frame[1, 8], frame[1, 9], frame[1, 10], frame[1, 11], frame[1, 12]))
  for (i in seq_len(nrow(frame))[-1]) {
    pelod2.score <- c(pelod2.score, c(HolisticPELOD2(frame[i, 13], frame[i, 3], frame[i, 4], frame[i, 5], frame[i, 6], frame[i, 7], frame[i, 8], frame[i, 9], frame[i, 10], frame[i, 11], frame[i, 12])))
  }
  return (cbind(frame, pelod2.score))
}

#Precondition: frame contains at least 1 row
NeuroPELOD2Frame <- function(frame) {
  pelod2.neuroscore <- c(NeuroPELOD2(frame[1, 3], frame[1, 4]))
  for (i in seq_len(nrow(frame))[-1]) {
    pelod2.neuroscore <- c(pelod2.neuroscore, c(NeuroPELOD2(frame[i, 3], frame[i, 4])))
  }
  return (cbind(frame, pelod2.neuroscore))
}

#Precondition: frame contains at least 1 row
CardioPELOD2Frame <- function(frame) {
  frame <- cbind(frame, pelod2.months = as.numeric(pelod2.data$pelod2.age) / 365 * 12)
  pelod2.cardioscore <- c(CardioPELOD2(frame[1, 13], frame[1, 5], frame[1, 6]))
  for (i in seq_len(nrow(frame))[-1]) {
    pelod2.cardioscore <- c(pelod2.cardioscore, c(CardioPELOD2(frame[i, 13], frame[i, 5], frame[i, 6])))
  }
  return (cbind(frame, pelod2.cardioscore))
}

#Precondition: frame contains at least 1 row
RenalPELOD2Frame <- function(frame) {
  frame <- cbind(frame, pelod2.months = as.numeric(pelod2.data$pelod2.age) / 365 * 12)
  pelod2.renalscore <- c(RenalPELOD2(frame[1, 13], frame[1, 7]))
  for (i in seq_len(nrow(frame))[-1]) {
    pelod2.renalscore <- c(pelod2.renalscore, c(RenalPELOD2(frame[i, 13], frame[i, 5], frame[i, 6])))
  }
  return (cbind(frame, pelod2.renalscore))
}

#Precondition: frame contains at least 1 row
RespPELOD2Frame <- function(frame) {
  pelod2.respscore <- c(RespPELOD2(frame[1, 8], frame[1, 9], frame[1, 10]))
  for (i in seq_len(nrow(frame))[-1]) {
    pelod2.respscore <- c(pelod2.respscore, c(RespPELOD2(frame[i, 8], frame[i, 9], frame[i, 10])))
  }
  return (cbind(frame, pelod2.respscore))
}

#Precondition: frame contains at least 1 row
HemaPELOD2Frame <- function(frame) {
  pelod2.hemascore <- c(HemaPELOD2(frame[1, 11], frame[1, 12]))
  for (i in seq_len(nrow(frame))[-1]) {
    pelod2.hemascore <- c(pelod2.hemascore, c(HemaPELOD2(frame[i, 11], frame[i, 12])))
  }
  return (cbind(frame, pelod2.hemascore))
}

ProbMortality <- function(pelod2) {
  return (1 / (1 + exp(6.61 - 0.47 * pelod2)))
}

FindExtremeValue <- function(id, frame) {
  event = (frame$Clinical.Event)[1]
  found.initial = F
  initial.index = 0
  final.index = 0
  for (i in seq_len(nrow(frame))) {
    if ((frame$encid)[i] == id && !found.initial)
      initial.index = i
    else if ((frame$encid)[i] == id && found.initial)
      final.index = i
  }
    
  if (event == "Peds Coma Score" || event == "Pupil Left Reaction" || event == "Pupil Right Reaction" || event == "Arterial Mean Pressure" || event == "PaO2" || event == "WBC" || event == "Platelets") {
    extremum = (frame$ #extracting with variable name instead of column name?
    for (i in initial.index:final.index) {
      
    }
  }
    
}

#Creating data frame to test algorithm
#pelod2.age measured in days
#pelod2.data <- data.frame( 
#  pelod2.id = c(1:10),
#  pelod2.age = c(seq(from = 10, to = 2710, length.out = 10)),
#  pelod2.gcoma = c(seq(from = 3, to = 12, length.out = 10)),
#  pelod2.pup = rep(c(1, 0), 5),
#  pelod2.lac = c(seq(from = 3.0, to = 16.5, length.out = 10)),
#  pelod2.map = c(seq(from = 15, to = 105, length.out = 10)),
#  pelod2.cr = c(seq(from = 20, to = 110, length.out = 10)),
#  pelod2.carrico = c(seq(from = 30, to = 75, length.out = 10)),
#  pelod2.paco2 = c(seq(from = 40, to = 130, length.out = 10)),
#  pelod2.vent = rep(c(1, 0), 5), 
#  pelod2.wbc = c(seq(from = 1, to = 2.8, length.out = 10)),
#  pelod2.plate = c(seq(from = 60, to = 150, length.out = 10))
#) 

#set.seed(100)
#pelod2.data <- data.frame(
#  pelod2.id = c(1:5000),
#  pelod2.age = rnorm(5000, 1200, 15),
#  pelod2.gcoma = rnorm(5000, 7, 1),
#  pelod2.pup = rbinom(5000, 1, 0.8),
#  pelod2.lac = rnorm(5000, 8, 1),
#  pelod2.map = rnorm(5000, 40, 7),
#  pelod2.cr = rnorm(5000, 80, 4),
#  pelod2.carrico = rnorm(5000, 60, 6),
#  pelod2.paco2 = rnorm(5000, 76, 8),
#  pelod2.vent = rbinom(5000, 1, 0.3),
#  pelod2.wbc = rnorm(5000, 1.6, 0.3),
#  pelod2.plate = rnorm(5000, 110, 11)
#) 

#HolisticPELOD2Frame(pelod2.data)
#NeuroPELOD2Frame(pelod2.data)
#ProbMortality(HolisticPELOD2Frame(pelod2.data)$pelod2.score)

#File names:
#admit1picu_deid_unencrypt.xlsx
#PELOD2_2015-2016_deid_unencrypt.xlsx
pelod2.data <- read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", 14)
pelod2.data2 <- read.xlsx("PELOD2_2015-2016_deid_unencrypt.xlsx", 15)
print(pelod2.data$Clinical.Event)