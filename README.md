# IVQR-GMM-computation-codes
Matlab codes for exact computation of the IVQR GMM estimator. Description of this estimator and its computation details can be found in the paper:

Chen, Le-Yu and Lee, Sokbae (March 2017), "Exact computation of GMM estimators for instrumental variable quantile regression models".

The paper is available at arXiv: https://arxiv.org/abs/1703.09382


The matlab function IVQR_GMM defind in IVQR_GMM.m can be used for exact computation of GMM estimators for instrumental variable quantile regression models (IVQR). Implementation of this function requires the Gurobi solver and another three functions (IVQR_MIO, miobnd_fn, and Two_stage_LS) whose implementation can be found in their corresponding matlab files stored in this repository. The Gurobi solver is freely available for academic purposes.  

