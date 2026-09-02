# Predicting socioeconomic status in childhood and adulthood from midlife intrinsic connectivity

### Code repository for:
Larsen et al. *The functional organization of the brain in midlife differentially reflects socioeconomic status in childhood and adulthood*.

Code used to tune, train, test, and interpret regression models that predict childhood and adulthood SES from midlife intrinsic connectivity, derived from task and resting-state fMRI.

### Background:
Emerging evidence suggests that socioeconomic differences are associated with the functional organization of the brain as it develops. This raises the hypothesis that alterations in the development of functional brain organization may be one means through which disadvantage becomes biologically embedded to later affect mental health, or a mechanism through which higher socioeconomic status confers advantages. The persistence of such differences in functional brain organization remains untested through longitudinal evaluation, as do potential differences across scales of socioeconomic status (SES; i.e., individual, neighborhood/area). 

Among members of a population-representative birth-cohort followed to midlife (the New Zealand-based Dunedin Study), we tested the association of both familial/individual and neighborhood socioeconomic status (SES) in childhood (birth–15y) and adulthood (26y–45y) with fMRI-assessed intrinsic whole-brain connectivity at age 45y (*N*=769; 49% female). 

### File directory:

#### 1. parameter_tuning_main_analyses.R
This file runs predictive models 100 times, varying model parameters (lambda) each time. This code loops across the six SES/time-    point combinations (childhood and adulthood, individual- and neighborhood-level SES. 

* *Inputs*: GFC edges per Study member (matrix), reliability (ICC) values per edge (vector), framewise displacement values per Study member (dataframe), behavioral dataframe with SES and other sociodemographic data (dataframe).

* *Outputs*: For each variable, outputs are created for a) base models and b) models with covariates added. Outputs are 100 saved enet objects (model).

2. 
