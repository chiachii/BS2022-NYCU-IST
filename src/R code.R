# LIBARIES ----
library(epitools)
library(VGAM)
library(corrplot)
library(VIM)
library(missMDA)
library(naniar)
library(epiDisplay)
library(DescTools)
library(class)
library(caret)
library(sampling)
library(dplyr)
library(DataExplorer)
library(explore)
library(tidyverse)
library(SmartEDA)
library(flextable)
library(dlookr)
library(performance)
library(factoextra)
library(car)

# DATA ----
data = read.csv("Desktop/NYCU Courses/111-1 Biostatistic/term project/data/data.csv", header=T)
data = tibble(data)
data

# Exploratory Data Analysis ----
## EDA Report ----
data %>% create_report(
  output_file   = "eda_report",
  output_dir    = "Desktop/NYCU Courses/111-1 Biostatistic/term project/result/",
  report_title  = "EDA report - In-hospital-mortality"
)

## outcome qqplot
plot_qq(data$outcome)


## Correlation ----
cor_df = correlate(data, outcome, method = "spearman")
cor_df = arrange(cor_df, desc(coef_corr))
print(cor_df, n=51)

## Outliers ----
## Missing value ----
gg_miss_var(data)
res.comp = imputePCA(
  subset(data, select = -c(group, ID, Admission.Date, Discharge.Date)),
  ncp=2
) # iterativePCA algorithm
imp_data = res.comp$completeObs # the imputed data

# Principal Component Analysis ----
pca = prcomp(
  formula = ~ .,
  data = subset(data, select = -c(group, ID, outcome, Admission.Date, Discharge.Date)),
  scale=T
)
pca

imp_data = data.frame(imp_data)
imp_pca = prcomp(
  formula = ~ .,
  data = subset(imp_data, select = -c(outcome)),
  scale=T
)
imp_pca

## Plot PCA (Scree plot & Pareto plot) ----
plot(pca$sdev, type="b", xlim=c(0,20), main="In-hospital-mortality (pca)")
abline(h=1, col="blue") # Kaiser eigenvalue-greater-than-one rule

plot(imp_pca$sdev, type="b", xlim=c(0,20), main="In-hospital-mortality (imp_pca)")
abline(h=1, col="blue") # Kaiser eigenvalue-greater-than-one rule

## 計算每個主成分的解釋比例＝各個主成分的特徵值/總特徵值
vars = (pca$sdev)^2
props = vars/sum(vars)
props

imp_vars = (imp_pca$sdev)^2
imp_props = imp_vars/sum(imp_vars)
imp_props

## Plot cumulative variance ----
### PCA (with missing value) ----
plot(cumsum(props), type="b")
abline(v=18, col="blue")
abline(h=cumsum(props[1:18])[18], col="blue")
legend(x=0, 
       y=cumsum(props[1:18])[18], 
       legend=sprintf("%3.3f", round(cumsum(props[1:18])[18], 3)),
       bty="n",
       cex=0.9,
       text.col="blue"
)

abline(v=31, col="red")
abline(h=cumsum(props[1:31])[31], col="red")
legend(x=0, 
       y=cumsum(props[1:31])[31], 
       legend=sprintf("%3.3f", round(cumsum(props[1:31])[31], 3)),
       bty="n",
       cex=0.9,
       text.col="red"
)

### PCA (without missing value) ----
plot(cumsum(imp_props), type="b")
abline(v=16, col="blue")
abline(h=cumsum(imp_props[1:16])[16], col="blue")
legend(x=0, 
       y=cumsum(imp_props[1:16])[16], 
       legend=sprintf("%3.3f", round(cumsum(imp_props[1:16])[16], 3)),
       bty="n",
       cex=0.9,
       text.col="blue"
)

abline(v=32, col="red")
abline(h=cumsum(imp_props[1:32])[32], col="red")
legend(x=0, 
       y=cumsum(imp_props[1:32])[32], 
       legend=sprintf("%3.3f", round(cumsum(imp_props[1:32])[32], 3)),
       bty="n",
       cex=0.9,
       text.col="red"
)

## Top18,31,49 PCA data
top18_pca.data <- pca$x[, 1:18]
top18_pca.data
top31_pca.data <- pca$x[, 1:31]
top31_pca.data
top49_pca.data <- pca$x[, 1:49]
top49_pca.data


# Categorical Data Analysis ----
## GLM (with missing value) ----
### logit regression analysis ----
m1 = glm(
  formula = outcome ~ .,
  data = subset(data, select = -c(group, ID, Admission.Date, Discharge.Date)),
  family = binomial
)
summary(m1) # significant level, AIC
PseudoR2(m1) # R-square
lroc(m1, grid.col="grey", lwd=3, auc.coords=c(0.4, 0.1)) # ROC, AUC

## pca18
idx = as.numeric(rownames(top18_pca.data))
score = as_tibble(factoextra::get_pca_ind(pca)$coord) 
mod18 <- cbind(data[idx,][3], score[1:18]) 
m1_pca18 = glm(
  formula = outcome ~ .,
  data = mod18,
  family = binomial
)
summary(m1_pca18)

## pca31
idx = as.numeric(rownames(top31_pca.data))
score = as_tibble(factoextra::get_pca_ind(pca)$coord) 
mod31 <- cbind(data[idx,][3], score[1:31]) 
m1_pca31 = glm(
  formula = outcome ~ .,
  data = mod31,
  family = binomial
)
summary(m1_pca31)

## pca49
idx = as.numeric(rownames(top49_pca.data))
score = as_tibble(factoextra::get_pca_ind(pca)$coord) 
mod49 <- cbind(data[idx,][3], score[1:49]) 
m1_pca49 = glm(
  formula = outcome ~ .,
  data = mod49,
  family = binomial
)
summary(m1_pca49)


### probit regression analysis ----
m2 = glm(
  formula = outcome ~ .,
  data = subset(data, select = -c(group, ID, Admission.Date, Discharge.Date)),
  family = binomial(link="probit")
)
summary(m2) # significant level, AIC
PseudoR2(m2) # R-square
lroc(m2, grid.col="grey", lwd=3, auc.coords=c(0.4, 0.1)) # ROC, AUC

## GLM (without missing value) ----
imp_m2 = glm(
  formula = outcome ~ .,
  data = imp_data,
  family = binomial(link="probit")
)
summary(imp_m2) # significant level, AIC
PseudoR2(imp_m2) # R-square
lroc(imp_m2, grid.col="grey", lwd=3, auc.coords=c(0.4, 0.1)) # ROC, AUC

## Features Selection ----
## training glm with removing missing value
data_omit = na.omit(data)
data_omit = subset(data_omit, select = -c(group, ID, Admission.Date, Discharge.Date))
omit_m2 = glm(
  formula = outcome ~ .,
  data = data_omit,
  family = binomial(link="probit")
)
summary(omit_m2) # significant level, AIC
PseudoR2(omit_m2) # R-square
lroc(omit_m2, grid.col="grey", lwd=3, auc.coords=c(0.3, 0.1)) # ROC, AUC

### Backward Selection ----
backward_m2 = step(omit_m2, direction="backward")
backward_m2
backward_m2$formula
summary(backward_m2)
lroc(backward_m2, grid.col="grey", lwd=3, auc.coords=c(0.3, 0.1)) # ROC, AUC

#### Correlation of Backward Selection ----
backward_df = subset(
  data_omit, 
  select = c(
    outcome, age, deficiencyanemias, Hyperlipemia, 
    Renal.failure, COPD, heart.rate, Diastolic.blood.pressure,
    Respiratory.rate, temperature, Platelets, Lymphocyte,
    Creatinine, Urea.nitrogen, Blood.sodium, Blood.calcium, 
    Chloride, Anion.gap, Magnesium.ion, PCO2
  )
)

'''
probit_outcome = as.list(data_omit[,1])
probit_outcome = probitlink(probit_outcome$outcome, bvalue = .Machine$double.eps)
backward_df[,1] = data.frame(probit_outcome)
'''
'''
logit_outcome = as.list(data_omit[,1])
logit_outcome = logitlink(logit_outcome$outcome, bvalue = .Machine$double.eps)
'''
backward_cor_df = correlate(backward_df, outcome, method = "spearman")
backward_cor_df = arrange(backward_cor_df, desc(coef_corr))
print(backward_cor_df)

M=cor(backward_df)
corrplot(M)

### Forward Selection ----
null_m2 = glm(
  formula = outcome ~ 1,
  data = data_omit,
  family = binomial(link="probit")
)
forward_m2 = step(
  null_m2, 
  # 從空模型開始，一個一個丟變數，
  # 最大不會超過完整的線性迴歸
  # (一定要加上界 upper=full，不可以不加) 
  scope=list(lower=null_m2, upper=omit_m2), 
  direction="forward"
)
summary(forward_m2)
lroc(forward_m2, grid.col="grey", lwd=3, auc.coords=c(0.3, 0.1))

### Forward Stepwise Selection ----
stepwise_m2 = step(omit_m2, direction="both")

## VIF ----
### GLM for VIF (test) ----
vif(backward_m2) # non-exist interaction

test = glm(
  formula = outcome ~ data_omit$hematocrit+
    data_omit$RBC+
    data_omit$MCH+
    data_omit$MCHC+
    data_omit$MCV+
    data_omit$PT+
    data_omit$INR+
    data_omit$Blood.sodium+
    data_omit$Chloride+
    data_omit$Anion.gap+
    data_omit$Bicarbonate+
    data_omit$PCO2,
  data = data_omit,
  family = binomial(link="probit")
)
vif(test)

### Interaction Terms (test) ----
interaction_backward_m2 = glm(
  #formula = outcome ~ . + 
    #Urea.nitrogen*Creatinine + 
    #Renal.failure*Creatinine +
    #Creatinine*Anion.gap +
    #Blood.sodium*Chloride +
    #Anion.gap*PCO2,
  formula = outcome ~ age + deficiencyanemias + Hyperlipemia + Renal.failure + 
    COPD + heart.rate + Diastolic.blood.pressure + Respiratory.rate + 
    temperature + Platelets + Lymphocyte + Creatinine + Urea.nitrogen + 
    Blood.sodium + Blood.calcium + Chloride + Anion.gap + Magnesium.ion + 
    PCO2 + Anion.gap*Blood.sodium*Chloride*PCO2,
    #hematocrit + RBC + MCH + MCHC + MCV + PT + INR + Blood.sodium + Chloride +
    #Anion.gap:Blood.sodium:Chloride:PCO2 + 
    #Blood.sodium:Chloride + Blood.sodium:PCO2 + Chloride:PCO2,
  #+ hematocrit*RBC*MCH*MCHC*MCV*PT*INR*Blood.sodium*Chloride*Anion.gap*Bicarbonate*PCO2,
  data = data_omit,
  family = binomial(link="probit")
) # choose vif>=10
summary(interaction_backward_m2)
lroc(interaction_backward_m2, grid.col="grey", lwd=3, auc.coords=c(0.4, 0.1)) # ROC, AUC

### Spearman Correlation Test for Interaction terms ----
cor.test(data_omit$Chloride, data_omit$PCO2, method=c("spearman"))
plot(data_omit$Chloride, data_omit$PCO2, col="orange", pch=19)

### Result ----
# backward + interaction + backward
ibm2 = step(interaction_backward_m2, direction="backward")
summary(ibm2)
lroc(ibm2, grid.col="grey", lwd=3, auc.coords=c(0.3, 0.1)) # ROC, AUC

# backward + interaction + forward
ifm2 = step(
  null_m2,
  scope=list(lower=null_m2, upper=interaction_backward_m2),
  direction="forward"
)
summary(ifm2)
lroc(ifm2, grid.col="grey", lwd=3, auc.coords=c(0.3, 0.1)) # ROC, AUC

# forward + interaction + backward
fullinter_m2 = glm(
  formula = outcome ~ . + Anion.gap*Blood.sodium*Chloride*PCO2,
  data = data_omit,
  family = binomial(link="probit")
)

finter_m2 = step(
  null_m2, 
  scope=list(lower=null_m2, upper=fullinter_m2), 
  direction="forward"
)
ifbm2 = step(finter_m2, direction="backward")
summary(ifbm2)
lroc(ifbm2, grid.col="grey", lwd=3, auc.coords=c(0.3, 0.1)) # ROC, AUC

# forward + interaction + forward
iffm2 = step(
  null_m2,
  scope=list(lower=null_m2, upper=finter_m2),
  direction="forward"
)
summary(iffm2)
lroc(iffm2, grid.col="grey", lwd=3, auc.coords=c(0.3, 0.1)) # ROC, AUC

# ————— ----
## PCA for backward selection ----
backward_pca = prcomp(
  formula = ~ .,
  data = subset(backward_df, select = -c(outcome)),
  scale=T
)
backward_pca

### Plot PCA (Scree plot & Pareto plot) ----
plot(backward_pca$sdev, type="b", xlim=c(0,20), main="In-hospital-mortality (backward_pca)")
abline(h=1, col="blue") # Kaiser eigenvalue-greater-than-one rule

## 計算每個主成分的解釋比例＝各個主成分的特徵值/總特徵值
backward_vars = (backward_pca$sdev)^2
backward_props = backward_vars/sum(backward_vars)
backward_props

### Plot cumulative variance ----
plot(cumsum(backward_props), type="b")
abline(v=7, col="blue")
abline(h=cumsum(backward_props[1:7])[7], col="blue")
legend(x=0, 
       y=cumsum(backward_props[1:7])[7], 
       legend=sprintf("%3.3f", round(cumsum(backward_props[1:7])[7], 3)),
       bty="n",
       cex=0.9,
       text.col="blue"
)

abline(v=14, col="red")
abline(h=cumsum(backward_props[1:14])[14], col="red")
legend(x=0, 
       y=cumsum(backward_props[1:14])[14], 
       legend=sprintf("%3.3f", round(cumsum(backward_props[1:14])[14], 3)),
       bty="n",
       cex=0.9,
       text.col="red"
)

## Top7, 14, 19 PCA Data ----
backward_top7_pca.data <- backward_pca$x[, 1:7]
backward_top14_pca.data <- backward_pca$x[, 1:14]
backward_top19_pca.data <- backward_pca$x[, 1:19]

## Model with backward selection PCA ----
## backward_top7_pca ----
idx = as.numeric(rownames(backward_top7_pca.data))
score = as_tibble(factoextra::get_pca_ind(backward_pca)$coord) 
backward_mod7 <- cbind(data[idx,][3], score[1:7]) 
backward_m2_pca7 = glm(
  formula = outcome ~ .,
  data = backward_mod7,
  family = binomial(link="probit")
)
summary(backward_m2_pca7)
PseudoR2(backward_m2_pca7) # R-square
lroc(backward_m2_pca7, grid.col="grey", lwd=3, auc.coords=c(0.4, 0.1)) # ROC, AUC

## backward_top14_pca ----
idx = as.numeric(rownames(backward_top14_pca.data))
score = as_tibble(factoextra::get_pca_ind(backward_pca)$coord) 
backward_mod14 <- cbind(data[idx,][3], score[1:14]) 
backward_m2_pca14 = glm(
  formula = outcome ~ .,
  data = backward_mod14,
  family = binomial(link="probit")
)
summary(backward_m2_pca14)
PseudoR2(backward_m2_pca14) # R-square
lroc(backward_m2_pca14, grid.col="grey", lwd=3, auc.coords=c(0.4, 0.1)) # ROC, AUC

## backward_top19_pca ----
idx = as.numeric(rownames(backward_top19_pca.data))
score = as_tibble(factoextra::get_pca_ind(backward_pca)$coord) 
backward_mod19 <- cbind(data[idx,][3], score[1:19]) 
backward_m2_pca19 = glm(
  formula = outcome ~ .,
  data = backward_mod19,
  family = binomial(link="probit")
)
summary(backward_m2_pca19)
PseudoR2(backward_m2_pca19) # R-square
lroc(backward_m2_pca19, grid.col="grey", lwd=3, auc.coords=c(0.4, 0.1)) # ROC, AUC


# ———Advances——— ----
## 2x2 Tables ----
### prep dataframe
counts <- c(
  sum(data_omit$outcome == 0 & data_omit$Renal.failure == 0),
  sum(data_omit$outcome == 0 & data_omit$Renal.failure == 1),
  sum(data_omit$outcome == 1 & data_omit$Renal.failure == 0),
  sum(data_omit$outcome == 1 & data_omit$Renal.failure == 1)
)
state <- gl(n = 2, k = 1, length = 4, labels = c("Alive", "Death"))
Renal_failure <- gl(n = 2, k = 2, length = 4, labels = c("No", "Yes"))
tab_df = data.frame(counts, state, Renal_failure)

conting_table <- xtabs(counts ~ state + Renal_failure)
conting_table

oddsratio(conting_table)
