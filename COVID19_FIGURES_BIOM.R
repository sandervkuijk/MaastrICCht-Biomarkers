################################################################################
### MaastrICCht cohort, analyses ECG
### Sander van Kuijk
###
### Doel: figuren associaties mortaliteit-biomarkers
###
### Start: 01/04/2021
### Laatste aanpassing: 30/07/2021
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

## Inlezen resultaten linear mixed-effects model
setwd("c:/Users/sande/Documents/Werk/MaastrICCht/DATA/PROJECT ECG")
load("crude_models.RData")
load("biomarkerdata.Rda")

## Transformatie biomarkers
d_long$logCK       <- log(d_long$CK)
d_long$logCKmb     <- log(d_long$CKmb)
d_long$loghstnt    <- log(d_long$hstnt)
d_long$logNTproBNP <- log(d_long$NTproBNP)

## Figuren opslaan onder ~MANUSCRIPTEN
setwd("c:/Users/sande/Documents/Werk/MaastrICCht/MANUSCRIPTEN/ECG")

## CK
png("logCK.png", width = 800, height = 600, pointsize = 18)

par(mar=c(5,4,1,2))

plot(jitter(d_long$meetdag[d_long$ICU_mortality == "Death"], 0.8),
     jitter(d_long$logCK[d_long$ICU_mortality == "Death"], 0.8),
     xlim = c(0, 5*7), ylim = c(0, max(d_long$logCK, na.rm = TRUE)),
     pch = 1, cex = 0.5, xlab = "Follow-up time (weeks)",
     ylab = "log-CK (U/L)", col = "grey50", xaxt = "n")
points(jitter(d_long$meetdag[d_long$ICU_mortality == "Alive"], 0.5),
       jitter(d_long$logCK[d_long$ICU_mortality == "Alive"], 0.5), pch = 3,
       col = "grey50", cex = 0.5)
axis(side = 1, at = seq(0, 35, 7), labels = 0:5)

b  <- fixef(CKi)
vb <- vcov(CKi)

days <- seq(0, 35, 7)
x <- t(sapply(days, function(x) c(1, 0, x, 0)))
pred1 <- x%*%b
x <- t(sapply(days, function(x) c(1, 1, x, x)))
pred0 <- x%*%b

lines(days, pred1, col = "black", lwd = 1, lty = 1)
lines(days, pred0, col = "black", lwd = 1, lty = 2)

legend("topright", col=c("black", "black"), lty=c(1, 2, 0, 0),
       lwd=c(1, 1, 1, 1), pch = c(NA, NA, 1, 3),
       legend=c("Non-survivor, estimated", "Survivor, estimated",
                "Non-survivor, raw data", "Survivor, raw data"), inset=.02,
       cex = 0.7)

dev.off()

## CKmb

png("logCKmb.png", width = 800, height = 600, pointsize = 18)

par(mar=c(5,4,1,2))

plot(jitter(d_long$meetdag[d_long$ICU_mortality == "Death"], 0.8),
     jitter(d_long$logCKmb[d_long$ICU_mortality == "Death"], 0.8),
     xlim = c(0, 5*7), ylim = c(0, max(d_long$logCKmb, na.rm = TRUE)),
     pch = 1, cex = 0.5, xlab = "Follow-up time (weeks)",
     ylab = expression(paste("log-CKmb (", mu, "g/L)")), col = "grey50", xaxt = "n")
points(jitter(d_long$meetdag[d_long$ICU_mortality == "Alive"], 0.5),
       jitter(d_long$logCKmb[d_long$ICU_mortality == "Alive"], 0.5), pch = 3,
       col = "grey50", cex = 0.5)
axis(side = 1, at = seq(0, 35, 7), labels = 0:5)

b  <- fixef(CKmbi)
vb <- vcov(CKmbi)

days <- seq(0, 35, 7)
x <- t(sapply(days, function(x) c(1, 0, x, 0)))
pred1 <- x%*%b
x <- t(sapply(days, function(x) c(1, 1, x, x)))
pred0 <- x%*%b

lines(days, pred1, col = "black", lwd = 1, lty = 1)
lines(days, pred0, col = "black", lwd = 1, lty = 2)

legend("topright", col=c("black", "black"), lty=c(1, 2, 0, 0),
       lwd=c(1, 1, 1, 1), pch = c(NA, NA, 1, 3),
       legend=c("Non-survivor, estimated", "Survivor, estimated",
                "Non-survivor, raw data", "Survivor, raw data"), inset=.02,
       cex = 0.7)

dev.off()

## hstnt

png("loghstnt.png", width = 800, height = 600, pointsize = 18)

par(mar=c(5,4,1,2))

plot(jitter(d_long$meetdag[d_long$ICU_mortality == "Death"], 0.8),
     jitter(d_long$loghstnt[d_long$ICU_mortality == "Death"], 0.8),
     xlim = c(0, 5*7), ylim = c(0, max(d_long$loghstnt, na.rm = TRUE)),
     pch = 1, cex = 0.5, xlab = "Follow-up time (weeks)",
     ylab = "log-hsTnT (ng/L)", col = "grey50", xaxt = "n")
points(jitter(d_long$meetdag[d_long$ICU_mortality == "Alive"], 0.5),
       jitter(d_long$loghstnt[d_long$ICU_mortality == "Alive"], 0.5), pch = 3,
       col = "grey50", cex = 0.5)
axis(side = 1, at = seq(0, 35, 7), labels = 0:5)

b  <- fixef(hstnti)
vb <- vcov(hstnti)

days <- seq(0, 35, 7)
x <- t(sapply(days, function(x) c(1, 0, x, 0)))
pred1 <- x%*%b
x <- t(sapply(days, function(x) c(1, 1, x, x)))
pred0 <- x%*%b

lines(days, pred1, col = "black", lwd = 1, lty = 1)
lines(days, pred0, col = "black", lwd = 1, lty = 2)

legend("topright", col=c("black", "black"), lty=c(1, 2, 0, 0),
       lwd=c(1, 1, 1, 1), pch = c(NA, NA, 1, 3),
       legend=c("Non-survivor, estimated", "Survivor, estimated",
                "Non-survivor, raw data", "Survivor, raw data"), inset=.02,
       cex = 0.7)

dev.off()

## NTproBNP

png("logNTproBNP.png", width = 800, height = 600, pointsize = 18)

par(mar=c(5,4,1,2))

plot(jitter(d_long$meetdag[d_long$ICU_mortality == "Death"], 0.8),
     jitter(d_long$logNTproBNT[d_long$ICU_mortality == "Death"], 0.8),
     xlim = c(0, 5*7), ylim = c(0, max(d_long$loghstnt, na.rm = TRUE)),
     pch = 1, cex = 0.5, xlab = "Follow-up time (weeks)",
     ylab = "log-NT-proBNP (pg/mL)", col = "grey50", xaxt = "n")
points(jitter(d_long$meetdag[d_long$ICU_mortality == "Alive"], 0.5),
       jitter(d_long$logNTproBNP[d_long$ICU_mortality == "Alive"], 0.5), pch = 3,
       col = "grey50", cex = 0.5)
axis(side = 1, at = seq(0, 35, 7), labels = 0:5)

b  <- fixef(ntproi)
vb <- vcov(ntproi)

days <- seq(0, 35, 7)
x <- t(sapply(days, function(x) c(1, 0, x, 0)))
pred1 <- x%*%b
x <- t(sapply(days, function(x) c(1, 1, x, x)))
pred0 <- x%*%b

lines(days, pred1, col = "black", lwd = 1, lty = 1)
lines(days, pred0, col = "black", lwd = 1, lty = 2)

legend("topright", col=c("black", "black"), lty=c(1, 2, 0, 0),
       lwd=c(1, 1, 1, 1), pch = c(NA, NA, 1, 3),
       legend=c("Non-survivor, estimated", "Survivor, estimated",
                "Non-survivor, raw data", "Survivor, raw data"), inset=.02,
       cex = 0.7)

dev.off()

### Einde
