
% This function computes the coefficient estimates and the estimated asymptotic variance
% for the two-stage least square (2SLS) regression of y on x using z as instruments

% function input :
% y     : the outcome vector
% x     : (n by k) matrix of covariate data where k is the number of covariates   
% z     : (n by p) matrix of instrument data where p is the number of instruments
% robust : set robust = 1 for calculating the estimated heteroskedasticity robust asymptotic variance


% function output :
% bhat : the vector of 2SLS regression coefficient estimates
% avar : the estimated 2SLS asymptotic variance 

function [bhat,avar] = Two_stage_LS(y,x,z,robust)
n=length(y);
k=size(x,2);
P = z*((z'*z)\(z'));
xhat=P*x;
bhat = (xhat'*xhat)\(xhat'*y);
uhat=y-x*bhat;

inv_xhat_xhat=(xhat'*xhat)\eye(k);

if robust==1 
  avar= inv_xhat_xhat*(xhat'*diag(uhat.*uhat)*xhat)*inv_xhat_xhat;
else
  avar= ((uhat'*uhat)/n)*inv_xhat_xhat;
end
end
