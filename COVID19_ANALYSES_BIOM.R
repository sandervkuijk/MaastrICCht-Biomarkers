################################################################################
### MaastrICCht cohort, analyses ECG
### Sander van Kuijk
###
### Doel: associaties Biomarkers IC-mortaliteit op ECG-cohort
###
### Start: 31/03/2021
### Laatste aanpassing: 02/04/2021
###
### Model 1 Crude
### Model 2 Sex, age
### Model 3 Sex, age, apache
### Model 4 Sex, age apache kreat dialyse
### Model 5 4 + cvd risk, roken, diabetes, obesitas
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

library(nlme)

## Databestand gecreeert met COVID19_DATABE_BIOM.R inlezen
setwd("c:/Users/sande/Documents/Werk/MaastrICCht/DATA/PROJECT ECG")
load("biomarkerdata.Rda")

## Databewerking
d_long$ICU_mortality <- factor(d_long$ICU_mortality, levels = c("Death", "Alive"))
d_long$meetmoment <- as.numeric(as.character(d_long$meetmoment))

## Check verdeling confounder
table(d_long$ckd_status.Dialysis_dependent)
table(d_long$cvrm.Smoking)

## Check descriptives biomarkers
summary(d_long$CK)
summary(d_long$CKmb)
summary(d_long$hstnt)
summary(d_long$NTproBNP)

## Transformatie biomarkers
d_long$logCK       <- log(d_long$CK)
d_long$logCKmb     <- log(d_long$CKmb)
d_long$loghstnt    <- log(d_long$hstnt)
d_long$logNTproBNP <- log(d_long$NTproBNP)

### Analyses CK
CK1 <- lme(logCK ~ ICU_mortality + meetmoment, data = d_long,
           random = ~1 | RecordId,
           na.action = na.omit, method = "ML",
           control = list(maxIter = 500, msMaxIter = 500, msVerbose = TRUE))
summary(CK1)

CK2 <- lme(logCK ~ ICU_mortality + meetmoment, data = d_long,
           random = ~1 + meetmoment | RecordId,
           na.action = na.omit, method = "ML",
           control = list(maxIter = 500, msMaxIter = 500, msVerbose = TRUE))
summary(CK2)

CK3 <- lme(logCK ~ ICU_mortality + meetmoment, data = d_long,
           random = ~1 + meetmoment | RecordId,
           correlation = corCAR1(form = ~ meetmoment | RecordId),
           na.action = na.omit, method = "ML",
           control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
           msVerbose = TRUE))
summary(CK3)
plot(CK3)

anova(CK1, CK2, CK3) # Meest complexe model het beste
# rm(CK1) #Emacs hangt bij rm?

## Interactie mortaliteit*meetmoment toevoegen
CKi <- lme(logCK ~ ICU_mortality*meetmoment, data = d_long,
           random = ~1 + meetmoment | RecordId,
           correlation = corCAR1(form = ~ meetmoment | RecordId),
           na.action = na.omit, method = "ML",
           control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
           msVerbose = TRUE))
summary(CKi)
round(intervals(CKi, which = "fixed")$fixed, 2)
plot(CKi) # Na transformatie niets meer op aan te merken

## Adjusted models, logCK
CKa1 <- lme(logCK ~ ICU_mortality*meetmoment + gender + age, data = d_long,
            random = ~1 + meetmoment | RecordId,
            correlation = corCAR1(form = ~ meetmoment | RecordId),
            na.action = na.omit, method = "ML",
            control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
            msVerbose = TRUE))
summary(CKa1)
round(intervals(CKa1, which = "fixed")$fixed, 2)

CKa2 <- lme(logCK ~ ICU_mortality*meetmoment +
            gender + age + APACHE_II, data = d_long,
            random = ~1 + meetmoment | RecordId,
            correlation = corCAR1(form = ~ meetmoment | RecordId),
            na.action = na.omit, method = "ML",
            control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
            msVerbose = TRUE))
summary(CKa2)
round(intervals(CKa2, which = "fixed")$fixed, 2)

CKa3 <- lme(logCK ~ ICU_mortality*meetmoment +
            gender + age + APACHE_II + creatinine_admission, data = d_long,
            random = ~1 + meetmoment | RecordId,
            correlation = corCAR1(form = ~ meetmoment | RecordId),
            na.action = na.omit, method = "ML",
            control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
            msVerbose = TRUE))
summary(CKa3)
round(intervals(CKa3, which = "fixed")$fixed, 2)

CKa4 <- lme(logCK ~ ICU_mortality*meetmoment +
            gender + age + APACHE_II + creatinine_admission +
            cvrm.Hypertension + cvrm.Diabetes_Mellitus +
            cvrm.Smoking + cvrm.Obesity, data = d_long,
            random = ~1 + meetmoment | RecordId,
            correlation = corCAR1(form = ~ meetmoment | RecordId),
            na.action = na.omit, method = "ML",
            control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
            msVerbose = TRUE))
summary(CKa4)
round(intervals(CKa4, which = "fixed")$fixed, 2)

### Analyses CKmb
CKmbi <- lme(logCKmb ~ ICU_mortality*meetmoment, data = d_long,
             random = ~1 + meetmoment | RecordId,
             correlation = corCAR1(form = ~ meetmoment | RecordId),
             na.action = na.omit, method = "ML",
             control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
             msVerbose = TRUE))
summary(CKmbi)
round(intervals(CKmbi, which = "fixed")$fixed, 2)
plot(CKmbi)

## Adjusted models
CKmba1 <- lme(logCKmb ~ ICU_mortality*meetmoment + gender + age, data = d_long,
              random = ~1 + meetmoment | RecordId,
              correlation = corCAR1(form = ~ meetmoment | RecordId),
              na.action = na.omit, method = "ML",
              control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
              msVerbose = TRUE))
summary(CKmba1)
round(intervals(CKmba1, which = "fixed")$fixed, 2)

CKmba2 <- lme(logCKmb ~ ICU_mortality*meetmoment +
              gender + age + APACHE_II, data = d_long,
              random = ~1 + meetmoment | RecordId,
              correlation = corCAR1(form = ~ meetmoment | RecordId),
              na.action = na.omit, method = "ML",
              control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
              msVerbose = TRUE))
summary(CKmba2)
round(intervals(CKmba2, which = "fixed")$fixed, 2)

CKmba3 <- lme(logCKmb ~ ICU_mortality*meetmoment +
              gender + age + APACHE_II + creatinine_admission, data = d_long,
              random = ~1 + meetmoment | RecordId,
              correlation = corCAR1(form = ~ meetmoment | RecordId),
              na.action = na.omit, method = "ML",
              control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
              msVerbose = TRUE))
summary(CKmba3)
round(intervals(CKmba3, which = "fixed")$fixed, 2)

CKmba4 <- lme(logCKmb ~ ICU_mortality*meetmoment +
              gender + age + APACHE_II + creatinine_admission +
              cvrm.Hypertension + cvrm.Diabetes_Mellitus +
              cvrm.Smoking + cvrm.Obesity, data = d_long,
              random = ~1 + meetmoment | RecordId,
              correlation = corCAR1(form = ~ meetmoment | RecordId),
              na.action = na.omit, method = "ML",
              control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
              msVerbose = TRUE))
summary(CKmba4)
round(intervals(CKmba4, which = "fixed")$fixed, 2)

### Analyses hstnt
hstnti <- lme(loghstnt ~ ICU_mortality*meetmoment, data = d_long,
              random = ~1 + meetmoment | RecordId,
              correlation = corCAR1(form = ~ meetmoment | RecordId),
              na.action = na.omit, method = "ML",
              control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
              msVerbose = TRUE))
summary(hstnti)
round(intervals(hstnti, which = "fixed")$fixed, 2)
plot(hstnti)

## Adjusted models
hstnta1 <- lme(loghstnt ~ ICU_mortality*meetmoment + gender + age, data = d_long,
               random = ~1 + meetmoment | RecordId,
               correlation = corCAR1(form = ~ meetmoment | RecordId),
               na.action = na.omit, method = "ML",
               control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
               msVerbose = TRUE))
summary(hstnta1)
round(intervals(hstnta1, which = "fixed")$fixed, 2)

hstnta2 <- lme(loghstnt ~ ICU_mortality*meetmoment +
               gender + age + APACHE_II, data = d_long,
               random = ~1 + meetmoment | RecordId,
               correlation = corCAR1(form = ~ meetmoment | RecordId),
               na.action = na.omit, method = "ML",
               control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
               msVerbose = TRUE))
summary(hstnta2)
round(intervals(hstnta2, which = "fixed")$fixed, 2)

hstnta3 <- lme(loghstnt ~ ICU_mortality*meetmoment +
               gender + age + APACHE_II + creatinine_admission +
               creatinine + dialysis_d, data = d_long,
               random = ~1 + meetmoment | RecordId,
               correlation = corCAR1(form = ~ meetmoment | RecordId),
               na.action = na.omit, method = "ML",
               control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
               msVerbose = TRUE))
summary(hstnta3)
round(intervals(hstnta3, which = "fixed")$fixed, 2)

hstnta4 <- lme(loghstnt ~ ICU_mortality*meetmoment +
               gender + age + APACHE_II + creatinine_admission +
               creatinine + dialysis_d +
               cvrm.Hypertension + cvrm.Diabetes_Mellitus +
               cvrm.Smoking + cvrm.Obesity, data = d_long,
               random = ~1 + meetmoment | RecordId,
               correlation = corCAR1(form = ~ meetmoment | RecordId),
               na.action = na.omit, method = "ML",
               control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
               msVerbose = TRUE))
summary(hstnta4)
round(intervals(hstnta4, which = "fixed")$fixed, 2)

### Analyses NTproBNP
ntproi <- lme(logNTproBNP ~ ICU_mortality*meetmoment, data = d_long,
              random = ~1 + meetmoment | RecordId,
              correlation = corCAR1(form = ~ meetmoment | RecordId),
              na.action = na.omit, method = "ML",
              control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
              msVerbose = TRUE))
summary(ntproi)
round(intervals(ntproi, which = "fixed")$fixed, 2)
plot(ntproi)

## Adjusted models
ntproa1 <- lme(logNTproBNP ~ ICU_mortality*meetmoment + gender + age, data = d_long,
               random = ~1 + meetmoment | RecordId,
               correlation = corCAR1(form = ~ meetmoment | RecordId),
               na.action = na.omit, method = "ML",
               control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
               msVerbose = TRUE))
summary(ntproa1)
round(intervals(ntproa1, which = "fixed")$fixed, 2)

ntproa2 <- lme(logNTproBNP ~ ICU_mortality*meetmoment +
               gender + age + APACHE_II, data = d_long,
               random = ~1 + meetmoment | RecordId,
               correlation = corCAR1(form = ~ meetmoment | RecordId),
               na.action = na.omit, method = "ML",
               control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
               msVerbose = TRUE))
summary(ntproa2)
round(intervals(ntproa2, which = "fixed")$fixed, 2)

ntproa3 <- lme(logNTproBNP ~ ICU_mortality*meetmoment +
               gender + age + APACHE_II + creatinine_admission, data = d_long,
               random = ~1 + meetmoment | RecordId,
               correlation = corCAR1(form = ~ meetmoment | RecordId),
               na.action = na.omit, method = "ML",
               control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
               msVerbose = TRUE))
summary(ntproa3)
round(intervals(ntproa3, which = "fixed")$fixed, 2)

ntproa4 <- lme(logNTproBNP ~ ICU_mortality*meetmoment +
               gender + age + APACHE_II + creatinine_admission +
               cvrm.Hypertension + cvrm.Diabetes_Mellitus +
               cvrm.Smoking + cvrm.Obesity, data = d_long,
               random = ~1 + meetmoment | RecordId,
               correlation = corCAR1(form = ~ meetmoment | RecordId),
               na.action = na.omit, method = "ML",
               control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
               msVerbose = TRUE))
summary(ntproa4)
round(intervals(ntproa4, which = "fixed")$fixed, 2)

## Resultaten modellen opslaan voor figuren
setwd("c:/Users/sande/Documents/Werk/MaastrICCht/DATA/PROJECT ECG")
save(CKi, CKmbi, hstnti, ntproi, file = "crude_models.RData")

### Einde
