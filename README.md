# IVQR-GMM-computation-codes
Matlab codes for exact computation of the IVQR GMM estimator. Description of this estimator and its computation details can be found in the paper:

Chen, Le-Yu and Lee, Sokbae (2018), "Exact computation of GMM estimators for instrumental variable quantile regression models".

The paper has been published at Journal of Applied Econometrics. See https://onlinelibrary.wiley.com/doi/full/10.1002/jae.2619. The latest working paper version of this work can be found in this repository.


The matlab function IVQR_GMM defind in IVQR_GMM.m can be used for exact computation of GMM estimators for instrumental variable quantile regression models (IVQR). Implementation of this function requires the Gurobi solver and another three functions (IVQR_MIO, miobnd_fn, and Two_stage_LS) whose implementation can be found in their corresponding matlab files stored in this repository. The Gurobi solver is freely available for academic purposes.  

