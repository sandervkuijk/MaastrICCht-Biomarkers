################################################################################
### MaastrICCht cohort, analyses ECG
### Sander van Kuijk, februari 2021
### R version 3.6.1
################################################################################

rm(list = ls())

library(nlme)
library(lme4)
library(lattice)
library(foreign)
library(reshape2)
library(data.table)

setwd("C:/Users/sande/Documents/Werk/MaastrICCht/DATA/SPSS data Moedi")

db <- read.spss("COVID ECG BASELINE.sav", to.data.frame = TRUE)
d1 <- read.spss("COVID ECG DAG1-3.sav", to.data.frame = TRUE)
d2 <- read.spss("COVID ECG DAG4-6.sav", to.data.frame = TRUE)
d3 <- read.spss("COVID ECG DAG7-10.sav", to.data.frame = TRUE)
d4 <- read.spss("COVID ECG DAG11-15.sav", to.data.frame = TRUE)
d5 <- read.spss("COVID ECG DAG16-20.sav", to.data.frame = TRUE)

colnames(d1)[-1] <- paste(colnames(d1)[-1], "dag1-3", sep = "_")
colnames(d2)[-1] <- paste(colnames(d2)[-1], "dag4-6", sep = "_")
colnames(d3)[-1] <- paste(colnames(d3)[-1], "dag7-10", sep = "_")
colnames(d4)[-1] <- paste(colnames(d4)[-1], "dag11-15", sep = "_")
colnames(d5)[-1] <- paste(colnames(d5)[-1], "dag16-20", sep = "_")

d  <- merge(db, d1, by = "Subject_ID", all.x = TRUE)
d  <- merge(d,  d2, by = "Subject_ID", all.x = TRUE)
d  <- merge(d,  d3, by = "Subject_ID", all.x = TRUE)
d  <- merge(d,  d4, by = "Subject_ID", all.x = TRUE)
d  <- merge(d,  d5, by = "Subject_ID", all.x = TRUE)

rm(db, d1, d2, d3, d4, d5)

names(d)

d_long <- melt(setDT(d), id.vars = "Subject_ID",
               measure.vars = list(grep("RheightV1",  names(d), value = TRUE),
                                   grep("RheightV2",  names(d), value = TRUE),
                                   grep("SdepthV1",   names(d), value = TRUE),
                                   grep("RRinterval", names(d), value = TRUE),
                                   grep("QTc",        names(d), value = TRUE)),
               value.name = c("RheightV1", "RheightV2", "SdepthV1",
                              "RRinterval", "QTc"),
               na.rm = TRUE)[order(Subject_ID)]
d_long <- data.frame(d_long)
names(d_long)[2] <- "meetmoment"
d_long$meetmoment <- factor(d_long$meetmoment, levels = c(1:6),
                     labels = c("baseline", "dag 1-3", "dag 4-6", "dag 7-10",
                                "dag 11-15", "dag 16-20"))
# View(d_long)

setwd("C:/Users/sande/Documents/Werk/MaastrICCht/DATA")

t <- read.csv("COVID19_ICU_250620_SOFA.csv", sep = ",")
t <- t[!duplicated(t$Record.Id),]
t <- t[t$Record.Id %in% levels(d$Subject_ID), ]

d_longt <- merge(d_long, t, by.x = "Subject_ID", by.y = "Record.Id")

d_longt$meetcont <- ifelse(d_longt$meetmoment == "baseline", 0,
                     ifelse(d_longt$meetmoment == "dag 1-3", 2,
                      ifelse(d_longt$meetmoment == "dag 4-6", 5,
                       ifelse(d_longt$meetmoment == "dag 7-10", 8.5,
                        ifelse(d_longt$meetmoment == "dag 11-15", 13,
                         ifelse(d_longt$meetmoment == "dag 16-20", 18, NA))))))
rm(t)

################################################################################
#### LME models
################################################################################

d_longt <- d_longt[!duplicated(d_longt[, c("Subject_ID", "meetmoment")]), ]
d_longt$meetcont_s <- scale(d_longt$meetcont, center = TRUE, scale = FALSE)

### Linear mixed-effects models, RheightV1

RheightV1.uni <- lme(RheightV1 ~ ICU_mortality + meetcont, data = d_longt,
                 random = ~ 1 + meetcont | Subject_ID,
                 correlation = corCAR1(form = ~ meetcont | Subject_ID),
                 na.action = na.omit, method = "ML",
                 control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                 msVerbose = TRUE))
summary(RheightV1.uni)

RheightV1.unii <- lme(RheightV1 ~ ICU_mortality*meetcont, data = d_longt,
                 random = ~ 1 + meetcont | Subject_ID,
                 correlation = corCAR1(form = ~ meetcont | Subject_ID),
                 na.action = na.omit, method = "ML",
                 control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                 msVerbose = TRUE))
summary(RheightV1.unii)
anova(RheightV1.uni, RheightV1.unii) # Model met interactie beter

RheightV1.adj <- lme(RheightV1 ~ ICU_mortality*meetcont + age + gender, data = d_longt,
                      random = ~ 1 + meetcont | Subject_ID,
                      correlation = corCAR1(form = ~ meetcont | Subject_ID),
                      na.action = na.omit, method = "ML",
                      control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                     msVerbose = TRUE))
summary(RheightV1.adj)
intervals(RheightV1.adj, which = "fixed") # Verschil op T = 0

# RheightV1.adj <- lme(RheightV1 ~ ICU_mortality*meetcont_s + age + gender, data = d_longt,
#                      random = ~ 1 + meetcont_s | Subject_ID,
#                      correlation = corCAR1(form = ~ meetcont_s | Subject_ID),
#                      na.action = na.omit, method = "ML",
#                      control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
#                                     msVerbose = TRUE))
# summary(RheightV1.adj)
# intervals(RheightV1.adj, which = "fixed") # Verschil op T = 9 (helft 0 - 18)

### Linear mixed-effects models, RheightV2

RheightV2.uni <- lme(RheightV2 ~ ICU_mortality + meetcont, data = d_longt,
                     random = ~ 1 + meetcont | Subject_ID,
                     correlation = corCAR1(form = ~ meetcont | Subject_ID),
                     na.action = na.omit, method = "ML",
                     control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                    msVerbose = TRUE))
summary(RheightV2.uni)

RheightV2.unii <- lme(RheightV2 ~ ICU_mortality*meetcont, data = d_longt,
                      random = ~ 1 + meetcont | Subject_ID,
                      correlation = corCAR1(form = ~ meetcont | Subject_ID),
                      na.action = na.omit, method = "ML",
                      control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                     msVerbose = TRUE))
summary(RheightV2.unii)
anova(RheightV2.uni, RheightV2.unii) # Model met interactie niet beter

RheightV2.adj <- lme(RheightV2 ~ ICU_mortality*meetcont + age + gender, data = d_longt,
                      random = ~ 1 + meetcont | Subject_ID,
                      correlation = corCAR1(form = ~ meetcont | Subject_ID),
                      na.action = na.omit, method = "ML",
                      control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                     msVerbose = TRUE))
summary(RheightV2.adj)
intervals(RheightV2.adj, which = "fixed")

### Linear mixed-effects models, SdepthV1

SdepthV1.uni <- lme(SdepthV1 ~ ICU_mortality + meetcont, data = d_longt,
                     random = ~ 1 + meetcont | Subject_ID,
                     correlation = corCAR1(form = ~ meetcont | Subject_ID),
                     na.action = na.omit, method = "ML",
                     control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                    msVerbose = TRUE))
summary(SdepthV1.uni)

SdepthV1.unii <- lme(SdepthV1 ~ ICU_mortality*meetcont, data = d_longt,
                      random = ~ 1 + meetcont | Subject_ID,
                      correlation = corCAR1(form = ~ meetcont | Subject_ID),
                      na.action = na.omit, method = "ML",
                      control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                     msVerbose = TRUE))
summary(SdepthV1.unii)
anova(SdepthV1.uni, SdepthV1.unii) # Model met interactie niet beter

SdepthV1.adj <- lme(SdepthV1 ~ ICU_mortality*meetcont + age + gender, data = d_longt,
                     random = ~ 1 + meetcont | Subject_ID,
                     correlation = corCAR1(form = ~ meetcont | Subject_ID),
                     na.action = na.omit, method = "ML",
                     control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                    msVerbose = TRUE))
summary(SdepthV1.adj)
intervals(SdepthV1.adj, which = "fixed")

### Linear mixed-effects models, RRinterval

RRinterval.uni <- lme(RRinterval ~ ICU_mortality + meetcont, data = d_longt,
                    random = ~ 1 + meetcont | Subject_ID,
                    correlation = corCAR1(form = ~ meetcont | Subject_ID),
                    na.action = na.omit, method = "ML",
                    control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                   msVerbose = TRUE))
summary(RRinterval.uni)

RRinterval.unii <- lme(RRinterval ~ ICU_mortality*meetcont, data = d_longt,
                     random = ~ 1 + meetcont | Subject_ID,
                     correlation = corCAR1(form = ~ meetcont | Subject_ID),
                     na.action = na.omit, method = "ML",
                     control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                    msVerbose = TRUE))
summary(RRinterval.unii)
anova(RRinterval.uni, RRinterval.unii) # Model met interactie niet beter

RRinterval.adj <- lme(RRinterval ~ ICU_mortality*meetcont + age + gender, data = d_longt,
                    random = ~ 1 + meetcont | Subject_ID,
                    correlation = corCAR1(form = ~ meetcont | Subject_ID),
                    na.action = na.omit, method = "ML",
                    control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                   msVerbose = TRUE))
summary(RRinterval.adj)
intervals(RRinterval.adj, which = "fixed")

### Linear mixed-effects models, QTc

QTc.uni <- lme(QTc ~ ICU_mortality + meetcont, data = d_longt,
                      random = ~ 1 + meetcont | Subject_ID,
                      correlation = corCAR1(form = ~ meetcont | Subject_ID),
                      na.action = na.omit, method = "ML",
                      control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                     msVerbose = TRUE))
summary(QTc.uni)

QTc.unii <- lme(QTc ~ ICU_mortality*meetcont, data = d_longt,
                     random = ~ 1 + meetcont | Subject_ID,
                     correlation = corCAR1(form = ~ meetcont | Subject_ID),
                     na.action = na.omit, method = "ML",
                     control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                    msVerbose = TRUE))
summary(QTc.unii)
anova(QTc.uni, QTc.unii) # Model met interactie niet beter

QTc.adj <- lme(QTc ~ ICU_mortality*meetcont + age + gender, data = d_longt,
                      random = ~ 1 + meetcont | Subject_ID,
                      correlation = corCAR1(form = ~ meetcont | Subject_ID),
                      na.action = na.omit, method = "ML",
                      control = list(opt = "optim", maxIter = 500, msMaxIter = 500,
                                     msVerbose = TRUE))
summary(QTc.adj)
intervals(QTc.adj, which = "fixed")

################################################################################
### Figures
################################################################################

setwd("C:/Users/sande/Documents/Werk/MaastrICCht/MANUSCRIPTEN/ECG")

png("RheightV1.png", width = 800, height = 600, pointsize = 18)
par(mfrow = c(1, 1))
par(mar=c(5,4,1,2))

plot(jitter(d_longt$meetcont[d_longt$ICU_mortality == "Death"], 0.8),
     jitter(d_longt$RheightV1[d_longt$ICU_mortality == "Death"], 0.8),
     xlim = c(0, max(d_longt$meetcont)), ylim = c(0, max(d_longt$RheightV1)),
     pch = 1, cex = 0.5, xlab = "Follow-up time (days)",
     ylab = "Height R-wave V1 (mV)", col = "grey50", xaxt = "n")
points(jitter(d_longt$meetcont[d_longt$ICU_mortality == "Alive"], 0.5),
       jitter(d_longt$RheightV1[d_longt$ICU_mortality == "Alive"], 0.5), pch = 3,
       col = "grey50", cex = 0.5)
axis(side = 1, at = unique(d_longt$meetcont))

b  <- fixef(RheightV1.adj)
vb <- vcov(RheightV1.adj)

days <- seq(0, 18, by = 1)
x <- t(sapply(days, function(x) c(1, 0, x, 0, 0, 0)))
pred0 <- x%*%b
x <- t(sapply(days, function(x) c(1, 1, x, 0, 0, x)))
pred1 <- x%*%b

lines(days, pred1, col = "black", lwd = 1, lty = 1)
lines(days, pred0, col = "black", lwd = 1, lty = 2)

legend("topright", col=c("black", "black"), lty=c(1, 2, 0, 0),
       lwd=c(1, 1, 1, 1), pch = c(NA, NA, 1, 3),
       legend=c("Non-survivor, estimated", "Survivor, estimated",
                "Non-survival, raw data", "Survivor, raw data"), inset=.02,
       cex = 0.7)

dev.off()


###

png("RheightV2.png", width = 800, height = 600, pointsize = 18)
par(mfrow = c(1, 1))
par(mar=c(5,4,1,2))

plot(jitter(d_longt$meetcont[d_longt$ICU_mortality == "Death"], 0.8),
     jitter(d_longt$RheightV2[d_longt$ICU_mortality == "Death"], 0.8),
     xlim = c(0, max(d_longt$meetcont)), ylim = c(0, 0.8),
     pch = 1, cex = 0.5, xlab = "Follow-up time (days)",
     ylab = "Height R-wave V2 (mV)", col = "grey50", xaxt = "n")
points(jitter(d_longt$meetcont[d_longt$ICU_mortality == "Alive"], 0.5),
       jitter(d_longt$RheightV2[d_longt$ICU_mortality == "Alive"], 0.5), pch = 3,
       col = "grey50", cex = 0.5)
axis(side = 1, at = unique(d_longt$meetcont))

b  <- fixef(RheightV2.adj)
vb <- vcov(RheightV2.adj)

days <- seq(0, 18, by = 1)
x <- t(sapply(days, function(x) c(1, 0, x, 0, 0, 0)))
pred0 <- x%*%b
x <- t(sapply(days, function(x) c(1, 1, x, 0, 0, x)))
pred1 <- x%*%b

lines(days, pred1, col = "black", lwd = 1, lty = 1)
lines(days, pred0, col = "black", lwd = 1, lty = 2)

legend("topright", col=c("black", "black"), lty=c(1, 2, 0, 0),
       lwd=c(1, 1, 1, 1), pch = c(NA, NA, 1, 3),
       legend=c("Non-survivor, estimated", "Survivor, estimated",
                "Non-survival, raw data", "Survivor, raw data"), inset=.02,
       cex = 0.7)

dev.off()

###

png("SdepthV1.png", width = 800, height = 600, pointsize = 18)
par(mfrow = c(1, 1))
par(mar=c(5,4,1,2))

plot(jitter(d_longt$meetcont[d_longt$ICU_mortality == "Death"], 0.8),
     jitter(d_longt$SdepthV1[d_longt$ICU_mortality == "Death"], 0.8),
     xlim = c(0, max(d_longt$meetcont)), ylim = c(0, 1.5),
     pch = 1, cex = 0.5, xlab = "Follow-up time (days)",
     ylab = "Depth S-wave in V1 (mV)", col = "grey50", xaxt = "n")
points(jitter(d_longt$meetcont[d_longt$ICU_mortality == "Alive"], 0.5),
       jitter(d_longt$SdepthV1[d_longt$ICU_mortality == "Alive"], 0.5), pch = 3,
       col = "grey50", cex = 0.5)
axis(side = 1, at = unique(d_longt$meetcont))

b  <- fixef(SdepthV1.adj)
vb <- vcov(SdepthV1.adj)

days <- seq(0, 18, by = 1)
x <- t(sapply(days, function(x) c(1, 0, x, 0, 0, 0)))
pred0 <- x%*%b
x <- t(sapply(days, function(x) c(1, 1, x, 0, 0, x)))
pred1 <- x%*%b

lines(days, pred1, col = "black", lwd = 1, lty = 1)
lines(days, pred0, col = "black", lwd = 1, lty = 2)

legend("topright", col=c("black", "black"), lty=c(1, 2, 0, 0),
       lwd=c(1, 1, 1, 1), pch = c(NA, NA, 1, 3),
       legend=c("Non-survivor, estimated", "Survivor, estimated",
                "Non-survival, raw data", "Survivor, raw data"),
       inset =.02, cex = 0.7)

dev.off()

###

png("RRinterval.png", width = 800, height = 600, pointsize = 18)
par(mfrow = c(1, 1))
par(mar=c(5,4,1,2))

plot(jitter(d_longt$meetcont[d_longt$ICU_mortality == "Death"], 0.8),
     jitter(d_longt$RRinterval[d_longt$ICU_mortality == "Death"], 0.8),
     xlim = c(0, max(d_longt$meetcont)), ylim = c(min(d_longt$RRinterval), max(d_longt$RRinterval)),
     pch = 1, cex = 0.5, xlab = "Follow-up time (days)",
     ylab = "RR interval (ms)", col = "grey50", xaxt = "n")
points(jitter(d_longt$meetcont[d_longt$ICU_mortality == "Alive"], 0.5),
       jitter(d_longt$RRinterval[d_longt$ICU_mortality == "Alive"], 0.5), pch = 3,
       col = "grey50", cex = 0.5)
axis(side = 1, at = unique(d_longt$meetcont))

b  <- fixef(RRinterval.adj)
vb <- vcov(RRinterval.adj)

days <- seq(0, 18, by = 1)
x <- t(sapply(days, function(x) c(1, 0, x, 0, 0, 0)))
pred0 <- x%*%b
x <- t(sapply(days, function(x) c(1, 1, x, 0, 0, x)))
pred1 <- x%*%b

lines(days, pred1, col = "black", lwd = 1, lty = 1)
lines(days, pred0, col = "black", lwd = 1, lty = 2)

legend("topright", col=c("black", "black"), lty=c(1, 2, 0, 0),
       lwd=c(1, 1, 1, 1), pch = c(NA, NA, 1, 3),
       legend=c("Non-survivor, estimated", "Survivor, estimated",
                "Non-survival, raw data", "Survivor, raw data"),
       inset =.02, cex = 0.7)

dev.off()

###

png("QTc.png", width = 800, height = 600, pointsize = 18)
par(mfrow = c(1, 1))
par(mar=c(5,4,1,2))

plot(jitter(d_longt$meetcont[d_longt$ICU_mortality == "Death"], 0.8),
     jitter(d_longt$QTc[d_longt$ICU_mortality == "Death"], 0.8),
     xlim = c(0, max(d_longt$meetcont)),
     ylim = c(min(d_longt$QTc), max(d_longt$QTc)),
     pch = 1, cex = 0.5, xlab = "Follow-up time (days)",
     ylab = "QTc time (ms)", col = "grey50", xaxt = "n")
points(jitter(d_longt$meetcont[d_longt$ICU_mortality == "Alive"], 0.5),
       jitter(d_longt$QTc[d_longt$ICU_mortality == "Alive"], 0.5), pch = 3,
       col = "grey50", cex = 0.5)
axis(side = 1, at = unique(d_longt$meetcont))

b  <- fixef(QTc.adj)
vb <- vcov(QTc.adj)

days <- seq(0, 18, by = 1)
x <- t(sapply(days, function(x) c(1, 0, x, 0, 0, 0)))
pred0 <- x%*%b
x <- t(sapply(days, function(x) c(1, 1, x, 0, 0, x)))
pred1 <- x%*%b

lines(days, pred1, col = "black", lwd = 1, lty = 1)
lines(days, pred0, col = "black", lwd = 1, lty = 2)

legend("topright", col=c("black", "black"), lty=c(1, 2, 0, 0),
       lwd=c(1, 1, 1, 1), pch = c(NA, NA, 1, 3),
       legend=c("Non-survivor, estimated", "Survivor, estimated",
                "Non-survival, raw data", "Survivor, raw data"),
       inset =.02, cex = 0.7)

dev.off()

################################################################################