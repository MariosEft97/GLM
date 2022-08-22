rm(list=ls())

#load libraries
#install.packages("DescTools")
library(DescTools) # PseudoR2
library(MASS)
library(plyr)
library(dplyr)
library(car)
library(DHARMa)
library(lmtest) # coeftest
library(sandwich) # coeftest argument vcov
#install.packages("countreg", repos="http://R-Forge.R-project.org") 
library(countreg) # rootogram

glm_RR = function(glm, digits=3){
  
  COEF=stats::coef(glm)
  CONFINT=stats::confint(glm)
  TABLE=cbind(coef=COEF, CONFINT)
  TABLE.EXP=round(exp(TABLE),digits)
  
  colnames(TABLE.EXP)[1]="RR"
  
  TABLE.EXP
}

# Data are imported into RStudio.
filepath = "C:/Users/35799/Desktop/KU Leuven/Semester 2/Generalized Linear Models/Project/16-victim.txt"
initial_data = read.csv(filepath, header=TRUE, sep = "\t", skip = 20)
head(initial_data)

# black: 0, white: 1
homicide_data=mutate(initial_data,race=unclass(as.factor(initial_data$race))-1)
head(homicide_data)

# Exploratory Data Analysis
hist(homicide_data$resp)
freq = table(homicide_data$race, homicide_data$resp)
barplot(freq[,2:7], col=c('black','grey'), legend=c('black','white'), 
        main='Barplot of victims known per race',
        ylab='number of respondents',
        xlab='number of victims known')
# Although black respondents are less than white ones they know more homicide victims.

# number of rows in the dataset per race
black_respondents=nrow(initial_data[initial_data$race=="black",])
white_respondents=nrow(initial_data[initial_data$race=="white",])

rm(initial_data)
attach(homicide_data)

# Scientific Question: Does race help explain how many homicide victims a person knows?

# 1. Fit a Poisson Model.
psm = glm(resp~race, family=poisson(link="log"), data=homicide_data)
summary(psm)
# both coefficients significant and negative

# 2. Calculate the risk ratio and the corresponding confidence interval.
glm_RR(psm)
#             RR    2.5 %  97.5 %
# (Intercept) 0.522 0.418  0.643
# race        0.177 0.133  0.236

# or
RR_black = exp(psm$coefficients[1])
RR_black
RR_white = exp(psm$coefficients[2])
RR_white

# 3. Calculate the ratio of the means of the response for each race.
mr_black = round(mean(resp[race==0]), 3)
mr_white = round(mean(resp[race==1]), 3)
mean_ratio = round(exp(psm$coefficients[1])/exp(psm$coefficients[1]+psm$coefficients[2]),3)
report = cbind(mr_black, mr_white, mean_ratio)
colnames(report)=c("black_mean", "white_mean", "means_ratio")
rownames(report)=c(1)
report
# The mean number of homicide victims known is 5.674 times greater for black people compared to that of
# white people. However we have to also take into account that the number of respondents from each race
# is different (159 black and 1149 white).

# 4. Calculate the predictions of the models for each race.
new_data = homicide_data[c("race")]
homicide_data_pred = homicide_data
homicide_data_pred$predicted = exp(predict(psm, new_data))
head(homicide_data_pred)
# the prediction of the model is equal to the relative risk of each race

# 5. Analyze the GOF of the model.

# Pearson test
x2_psm=sum(residuals(psm, type = "pearson")^2)
n_psm=dim(homicide_data)[1]
p_psm=length(coef(psm))
data.frame(x2s=x2_psm,pvalue=(1-pchisq(x2_psm,n_psm-p_psm)))
# p-value is zero and thus null hypothesis is rejected.

# Deviance tests
dev_psm=summary(psm)$deviance
df_psm=summary(psm)$df.residual
data.frame(dev=dev_psm, df=df_psm, pvalue=(1-pchisq(dev_psm,df_psm)))
# No evidence against the null hypothesis.

# Likelihood ratio test
Anova(psm, test="LR",type=3)
# the impact of race is highly significant

# Wald test
Anova(psm, test="Wald",type=3)
# the impact of race is highly significant
anova(psm, test ="Rao")
# all tests for the nested models confirm that indeed, race is important for the analysis

# Residuals
psm_res=simulateResiduals(psm, plot = T)
# ks.test borderline

# Uniformity of Residuals and Dispersion
hist(psm_res)
# Residuals are clearly not following a uniform distribution.

rootogram(psm, ylab="Root Square of Frequency", main="Poisson")

testUniformity(psm)

testDispersion(psm)
# We have clear indications of over-dispersion.

# 6. Fit a negative binomial model
nbm = glm.nb(resp~race, data=homicide_data)
summary(nbm)
testDispersion(nbm)
rootogram(nbm, ylab="Root Square of Frequency", main="Negative Binomial")

# Calculate model based variance
variance_black = exp(nbm$coefficients[1]) + (1/nbm$theta)*exp(nbm$coefficients[1])^2
variance_white = exp(nbm$coefficients[1] + nbm$coefficients[2]) + (1/nbm$theta)*exp(nbm$coefficients[1] + nbm$coefficients[2])^2

# Calculate observed variances
obs_var_black = var(homicide_data[which(homicide_data$race == 0),1])
obs_var_white = var(homicide_data[which(homicide_data$race == 1),1])

comparison = cbind(rbind(round(obs_var_black,3), round(obs_var_white,3)),rbind(round(variance_black,3),round(variance_white,3)), rbind(mr_black,mr_white))
colnames(comparison) = c("observed","pred_nb","pred_poisson")
rownames(comparison)=c("black", "white")
comparison

# 7. Fit a Quasi-likelihood model
qlm = glm(resp~race, family=quasipoisson, data=homicide_data)
summary(qlm)
#same coefficients as before st.errors changed
#variance is 74% larger than the mean
sqrt(1.745694)
#standard errors are about 32% larger than before

# 8. Discuss all results

round(data.frame(
  Po=coef(psm),QL=coef(qlm),NB=coef(nbm),
  se.Po=summary(psm)$coefficients[, 2],
  se.QL=summary(qlm)$coefficients[, 2],
  se.NB=summary(nbm)$coefficients[, 2]),4)


psm_metrics = PseudoR2(psm, which = 'all')
nbm_metrics =PseudoR2(nbm, which ='all')

comparison = cbind(rbind(psm_metrics[10], nbm_metrics[10]),rbind(psm_metrics[11],psm_metrics[11]))
colnames(comparison) = c("AIC","BIC")
rownames(comparison)=c("Poisson", "NB")
comparison
