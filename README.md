# BS2022(NYCU-IST)
Biostatistics (2022) at National Yang Ming Chiao Tung University (NYCU) IST (生物統計, 陽明交大統計所). 

## Overview
This project aims to predict in-hospital mortality for ICU-admitted heart failure patients using the MIMIC-III dataset. We employed survival analysis and categorical data analysis to identify key factors influencing patient outcomes. The dataset contains information on 1177 ICU patients, including demographic, clinical, and laboratory features. Various statistical models, including Cox regression and generalized linear models (GLMs), were applied to discover potential relationships and interactions between the variables and patient mortality.

## Contents
- **Dataset**: We utilized data from the MIMIC-III database, focusing on 1177 ICU heart failure patients. After dealing with missing values, the final dataset consists of 425 patients.
- **Survival Analysis**: We constructed Cox proportional hazards models to estimate the hazard ratios of various risk factors and assessed model performance using metrics such as concordance and Akaike Information Criterion (AIC).
- **Categorical Data Analysis**: Generalized Linear Models (GLMs) were built to assess the importance of specific variables, with forward and backward selection processes used to refine model performance.

## Results
We obtained two GLMs with better adaptability on this dataset and the following analysis results accordingly.
- Identified significant 19 factors such as Age, PL (血小板), NTproBNP (利鈉肽), K+, and Ca2+ have medically reasonable results on mortality in patients with heart failure.
- The coefficient of Renal.failure (腎衰竭), COPD (慢性阻塞性肺病), and Creatinine (肌酸酐) have unreasonable relationship with mortality in our model. In further analysis, we found that this unreasonable phenomenon comes from the source of this data, because this data is collected from ICU patients, it is possible that the patient was admitted to the ICU because of renal failure, so focus on the relationship between Renal.failure and target variables is meaningless.

<!--
## Term Project
-->
