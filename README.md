# Predicting socioeconomic status in childhood and adulthood from midlife intrinsic connectivity

### Code repository for:
Larsen et al. *The functional organization of the brain in midlife differentially reflects socioeconomic status in childhood and adulthood*.

Code used to tune, train, test, and interpret regression models that predict childhood and adulthood SES from midlife intrinsic connectivity, derived from task and resting-state fMRI.

### Background:
Emerging evidence suggests that socioeconomic differences are associated with the functional organization of the brain as it develops. This raises the hypothesis that alterations in the development of functional brain organization may be one means through which disadvantage becomes biologically embedded to later affect mental health, or a mechanism through which higher socioeconomic status confers advantages. The persistence of such differences in functional brain organization remains untested through longitudinal evaluation, as do potential differences across scales of socioeconomic status (SES; i.e., individual, neighborhood/area). 

Among members of a population-representative birth-cohort followed to midlife (the New Zealand-based Dunedin Study), we tested the association of both familial/individual and neighborhood socioeconomic status (SES) in childhood (birth–15y) and adulthood (26y–45y) with fMRI-assessed intrinsic whole-brain connectivity at age 45y (*N*=769; 49% female). 

### File directory:

#### 1. parameter_tuning_main_analyses.R
This file runs regularized regression models 100 times, varying the model regularization hyperparameter used to control the penalty strength each time (lambda). This code loops across the six SES/timepoint combinations (childhood and adulthood, individual- and neighborhood-level SES. 

* *Inputs*: GFC edges per Study member (matrix), reliability (ICC) values per edge (vector), framewise displacement values per Study member (dataframe), behavioral dataframe with SES and other sociodemographic data (dataframe).

* *Outputs*: For each variable, outputs are created for a) base models and b) models with covariates added. Outputs are 100 saved enet objects (model) for each variable. Additionally, within each iteration (per variable), the splits for 90/10 training/test are saved so they can be reloaded when the models are run. From the model outputs, the optimal lambda is chosen per variable as the one most frequently selected by the cross-validation using *caret()*.


#### 2. parameter_tuning_additional_analyses.R
This file contains additional tuning loops for two models: one predicting adult individual-level SES covarying for childhood individual-SES, and a second predicting adult neighborhood-level SES covarying for childhood neighborhood-level SES. Otherwise identical to the parameter_tuning_main_analyses file.

#### 3. predict_fc_full_loop.R
This file contains code for the regularized regression training and testing, using the hyperparameter (lambda) selected during the tuning steps above. This code loops across the six SES/timepoint combinations (childhood and adulthood, individual- and neighborhood-level SES. 

* *Inputs*: GFC edges per Study member (matrix), reliability (ICC) values per edge (vector), framewise displacement values per Study member (dataframe), behavioral dataframe with SES and other sociodemographic data (dataframe), 90/10 train/test splits generated during tuning, and lambda values selected during tuning.

* *Outputs*: For each variable, outputs are created for a) base models and b) models with covariates added. Outputs are: 1) a dataframe of model performance metrics extracted from the model output from *caret* function *predict()* in the test data. Performance metrics include the RMSE, R-squared value, MAE, and *r* (the correlation between observed SES values and SES values predicted from the model. 2) Haufe-transformed feature importance scores.

#### 4. predict_fc_additional_analyses.R
This file contains code for the regularized regression training and testing for adult SES models covarying for childhood SES within the same level (individual- or neighborhood-). Otherwise identical to predict_fc_full_loop.R.

#### 5. inspecting_performance.R
This file contains two loops, one for the main analyses, and one for the additional analyses predicting adult SES covarying for childhood SES. The loops load model output from each folder and extract model performance statistics.
