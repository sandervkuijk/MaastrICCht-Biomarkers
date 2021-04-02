################################################################################
### MaastrICCht cohort, analyses ECG
### Sander van Kuijk
###
### Doel: data uit verschillende bronnen inlezen en koppelen, naar LONG
###
### Start: 31/03/2021
### Laatste aanpassing: 31/03/2021
###
### sessionInfo()
###
### R version 4.0.4 (2021-02-15)
### Platform: x86_64-w64-mingw32/x64 (64-bit)
### Running under: Windows 10 x64 (build 19042)
###
### Matrix products: default
###
### locale:
### [1] LC_COLLATE=English_Netherlands.1252  LC_CTYPE=English_Netherlands.1252
### [3] LC_MONETARY=English_Netherlands.1252 LC_NUMERIC=C
### [5] LC_TIME=English_Netherlands.1252
###
### attached base packages:
### [1] stats     graphics  grDevices utils     datasets  methods   base
###
### loaded via a namespace (and not attached):
### [1] compiler_4.0.4
################################################################################

rm(list = ls())

library(foreign)
library(reshape2)
library(data.table)

# Functie om uit SPSS geimporteerde datum variabelen te transformeren
spss2date <- function(x) as.Date(x/86400, origin = "1582-10-14")

# Verwijzen naar map waarin data staan, en data uit SPSS-bestand inlezen
setwd("C:/Users/sande/Documents/Werk/MaastrICCht/DATA/PROJECT ECG")
d <- read.spss("COVID_19MUMC_CARDIAC BIOMARKERS_BROAD_1.sav", to.data.frame = TRUE)

# Inspectie data: unieke patientnummers
length(d$RecordId)
length(unique(d$RecordId))

# Datum variabelen transformeren
d[, grep("ReportCreationDate.", names(d), value = TRUE)]  <-
    lapply(d[, grep("ReportCreationDate.", names(d), value = TRUE)],
    FUN = function (x) as.Date(substr(x, 1, 10), format = c("%d-%m-%Y")))

# Verwijderen redundante variabelen
drops <- c("V10", grep("Echtedatum.", names(d), value = TRUE),
           grep("ReportNameCustom.", names(d), value = TRUE),
           grep("ReportParent.", names(d), value = TRUE))
d <- d[, !names(d) %in% drops]
rm(drops, spss2date)

# Data van wide naar long
d_long <- melt(setDT(d), id.vars = "RecordId",
               measure.vars = list(grep("ReportCreationDate.",  names(d), value = TRUE),
                                   grep("CK_daily.", names(d), value = TRUE),
                                   grep("hstnt.", names(d), value = TRUE),
                                   grep("CKmb.", names(d), value = TRUE),
                                   grep("NTproBNP.", names(d), value = TRUE)),
               value.name = c("Date", "CK", "hstnt", "CKmb", "NTproBNP"),
               na.rm = FALSE)[order(RecordId)]
d <- data.frame(d)
d_long <- data.frame(d_long)
names(d_long)[2] <- "meetmoment"

# Eerste datum per patient, om relatieve tijd te berekenen
first <- data.frame(d$RecordId, d$ReportCreationDate.1)
names(first) <- c("RecordId", "Day1")
d_long <- merge(d_long, first, by = "RecordId")
d_long$meetdag <- as.numeric(difftime(d_long$Date, d_long$Day1, units = "days"))
d_long <- subset(d_long, !is.na(d_long$Date))
rm(first)

# Koppelen met SOFA data voor ajustment variabelen
setwd("C:/Users/sande/Documents/Werk/MaastrICCht/DATA")

t <- read.csv("COVID19_ICU_250620_SOFA.csv", sep = ",")
names(t)[1] <- "RecordId"
t <- t[!duplicated(t$RecordId),]
t <- t[t$RecordId %in% unique(d$RecordId), ]

t$dialysis[t$dialysis == ""] <- NA

keep <- c("RecordId", "gender", "age", "BMI", "APACHE_II",
          "creatinine_admission", "dialysis",
          "cvrm.Hypertension", "cvrm.Diabetes_Mellitus",
          "cvrm.Smoking", "cvrm.Obesity", "ICU_mortality")
t <- t[, names(t) %in% keep]

d_long <- merge(d_long, t, by = "RecordId")
rm(t, keep)

# Data opslaan om te modelleren
setwd("c:/Users/sande/Documents/Werk/MaastrICCht/DATA/PROJECT ECG")
save(d_long, file = "biomarkerdata.Rda")

### Einde file.
