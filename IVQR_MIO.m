
% function input :
% y : vector of outcomes
% x     : (n by k) matrix of the covariate data
% Q     : (n by n ) matrix equal to (G*Q_hat*G') stated in the MIQP formulation 
% tau   : quantile index
% T     : the time limit specified for early termination of the MIO solver
% abgap : the absolute gap specified for early termination of the MIO solver
% bnd   : (k by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients

% function output :
% bhat  : the vector of the coefficient estimates
% obj_v : the value of the GMM objective function
% gap   : the MIO optimization gap value in case of early termination
%         gap = 0 ==> optimal solution is found within the time limit
% rtime : the time used by the MIO solver in the estimation procedure
% ncount: the number of nodes already explored by the MIO solver 

function [bhat,obj_v,gap,rtime,ncount] = IVQR_MIO(y,x,Q,tau,T,abgap,bnd)

n=length(y);
k=size(x,2);
bhat=zeros(k,1);

gap=0;
rtime=0;
ncount=0;

model.sense = '<';
model.modelsense = 'min';

model.lb = [zeros(n,1); bnd(:,1)];
model.ub = [ones(n,1); bnd(:,2) ];

% 'B' : int code 66
% 'C' : int code 67
model.vtype = char([66*ones(1,n) 67*ones(1,k)]); 

tol=1e-6;

miobnd=miobnd_fn(y,x,bnd);  % miobnd_fn computes the values M(i) defined in (3.6) of Chen and Lee (2017)
miobnd_bar = miobnd+tol;

params.outputflag = 0; 
params.OptimalityTol=tol;
params.FeasibilityTol=tol;
params.IntFeasTol=tol;

if T > 0
params.TimeLimit =T;
end

if abgap > 0
params.MIPGapAbs=abgap;
end

 tau_vec=ones(n,1)*tau;

 model.obj = [-2*Q*tau_vec;zeros(k,1)];
 model.objcon = tau_vec'*Q*tau_vec;
 model.Q=sparse([Q zeros(n,k);zeros(k,n+k)]);
try
    model.A = sparse([diag(miobnd) -x; -diag(miobnd_bar) x]);
    model.rhs = [miobnd*(1-tol)-y;y-tol*miobnd_bar];
    
    result = gurobi(model, params);

    bhat=result.x(n+1:n+k);
    obj_v=result.objval;
    gap=(obj_v-result.objbound);
    rtime=result.runtime;
    ncount=result.nodecount;
    fprintf('Optimization returned status: %s\n', result.status);
   
catch gurobiError
    fprintf('Error reported\n');
end

end

