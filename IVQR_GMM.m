
% The IVQR_GMM function computes the exact GMM estimator of the IVQR model
% via the MIO approach as described in Chen and Lee (2017).

% function input :
% y : vector of outcomes
% w     : (n by k) matrix of the covariate dataset
% z     : (n by p ) matrix of the instrument variable dataset
% tau   : quantile index
% T     : scalar
%         If T>0, then T is the time limit specified for early termination
%         of the MIO solver. Otherwise, the MIO solver keeps running until
%         convergence.
% abgap : the absolute gap specified for early termination of the MIO solver
% bnd   : (k by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients         

% The arguments T, abgap and bnd are optional. When they are not specified,
% the following default values are used.
% T : set T = 0  ==> solve the MIO problem until convergence
% abgap : set abgap = 0  ==> solve the MIO problem until convergence
% bnd : Calculate the parameter bounds based on the two-stage least square
%       regression results as used in Chen and Lee (2017)

% function output :
% theta_hat  : the vector of the coefficient estimates
% s_hat : the estimated asymptotic standard errors
% obj_v : the value of the GMM objective function
% gap   : the MIO optimization gap value in case of early termination
%         gap = 0 ==> optimal solution is found within the time limit
% rtime : the time used by the MIO solver in the estimation procedure
% ncount: the number of nodes already explored by the MIO solver 

function [theta_hat,s_hat,obj_v,gap,rtime,ncount] = IVQR_GMM(y,w,z,tau,T,abgap,bnd)

switch nargin
    case 4
        T=0; abgap=0; bnd=[];
    case 5
        abgap=0; bnd=[];
    case 6
        bnd = [];
    otherwise
        if nargin~=7
        disp('error in input arguments');
        return;
        end
end

method=1; % Computing the GMM estimator based on the formulation (3.3) of Chen and Lee (2017)
          % set method = 2 for using the formulation (C.1)
          % set method = 3 for using the formulation (C.10)


n=length(y);

% Use of the Hall-Sheath bandwidth choice (Koenker 1994)
q_tau = norminv(tau);
H_S_ratio = (normpdf(q_tau)^2)/(1+2*q_tau*q_tau);
h_Hall_Sheath =(1.96^(2/3))*((1.5/n*H_S_ratio)^(1/3));

Q=z*inv(z'*z/n)*(z')/(tau*(1-tau)); 
% Q is the matrix G*Q_hat*G' stated in the MIQP formulation of the GMM
% estimation problem

k=size(w,2);
theta_hat = zeros(k,1);

if isempty(bnd)
[b,var] = Two_stage_LS(y,w,z,1);
bnd=[b-10*sqrt(diag(var)) b+10*sqrt(diag(var))];
end

[theta_hat,obj_v,gap,rtime,ncount] = IVQR_MIO(y,w,Q,tau,T,abgap,bnd,method);

% compute the estimated asymptotic standard errors based on the method of
% Powell (1986) using Gaussian kernal and the Hall-Sheath bandwidth choice
e_hat = y-w*theta_hat;
kern=normpdf(e_hat/h_Hall_Sheath)/h_Hall_Sheath;
k_x = repmat(kern,1,k).*w;
s_hat = sqrt(diag(inv(k_x'*z*inv(z'*z)*z'*k_x)*tau*(1-tau)));
end


