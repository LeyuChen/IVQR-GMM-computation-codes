
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
% method : set method=1 for solving the MIQP formulation (3.3)
%          set method=2 for solving the MIQP formulation (C.1)
%          set method=3 for solving the MILP formulation (C.10)

% function output :
% bhat  : the vector of the coefficient estimates
% obj_v : the value of the GMM objective function
% gap   : the MIO optimization gap value in case of early termination
%         gap = 0 ==> optimal solution is found
% rtime : the time used by the MIO solver in the estimation procedure
% ncount: the number of nodes already explored by the MIO solver 

function [bhat,obj_v,gap,rtime,ncount] = IVQR_MIO(y,x,Q,tau,T,abgap,bnd,method)

n=length(y);
k=size(x,2);
bhat=zeros(k,1);

gap=0;
rtime=0;
ncount=0;

tau_vec=ones(n,1)*tau;
model.objcon = tau_vec'*Q*tau_vec;

model.modelsense = 'min';
tol=1e-6;

if method==1 % MIQP formulation (3.3)
disp('solving MIQP formulation (3.3)');  

model.sense = '<';
model.lb = [zeros(n,1); bnd(:,1)];
model.ub = [ones(n,1); bnd(:,2) ];

% 'B' : int code 66
% 'C' : int code 67
model.vtype = char([66*ones(1,n) 67*ones(1,k)]); 

miobnd=miobnd_fn(y,x,bnd);  % miobnd_fn computes the values M(i) defined in (3.6)
miobnd_bar = miobnd+tol;

model.obj = [-2*Q*tau_vec;zeros(k,1)];
model.Q=sparse([Q zeros(n,k);zeros(k,n+k)]);
model.A = sparse([diag(miobnd) -x; -diag(miobnd_bar) x]);
model.rhs = [miobnd*(1-tol)-y;y-tol*miobnd_bar];

elseif method==2 % MIQP formulation (C.1)
disp('solving MIQP formulation (C.1)');   

model.sense = [repmat('=',1,2*n) repmat('<',1,n)];
model.lb = [zeros(n,1); bnd(:,1);zeros(3*n,1)];
model.ub = [ones(n,1); bnd(:,2);ones(n,1);1/eps*ones(2*n,1)];

% 'B' : int code 66
% 'C' : int code 67
model.vtype = char([66*ones(1,n) 67*ones(1,k) 66*ones(1,n) 67*ones(1,2*n)]); 

model.obj = [-2*Q*tau_vec;zeros(3*n+k,1)];
model.Q=sparse([Q zeros(n,3*n+k);zeros(3*n+k,4*n+k)]);

model.A = sparse([eye(n) zeros(n,k) eye(n) zeros(n,2*n);zeros(n) x zeros(n) eye(n) -eye(n);zeros(n,2*n+k) -eye(n) -eye(n)]);
model.rhs = [ones(n,1);y;(-1e-5)*ones(n,1)];

params.PreSOS1BigM=0; 

% Specification of the SOS-1 constraints

for j=1:n
    model.sos(j).type   = 1;
    model.sos(j).index  = [j 2*n+k+j]; % (r,e): SOS-1
    model.sos(n+j).type   = 1;
    model.sos(n+j).index  = [n+k+j 3*n+k+j]; % (s,1-e): SOS-1
end
   
elseif method == 3 % MILP formulation (C.10)
disp('solving MILP formulation (C.10)');   

aux_num=n*(n-1)/2;

model.sense = '<';
model.lb = [zeros(n,1); bnd(:,1); zeros(aux_num,1)];
model.ub = [ones(n,1); bnd(:,2); ones(aux_num,1)];

% 'B' : int code 66
% 'C' : int code 67
model.vtype = char([66*ones(1,n) 67*ones(1,k) 66*ones(1,aux_num)]);

miobnd=miobnd_fn(y,x,bnd);  % miobnd_fn computes the values M(i) defined in (3.6)
miobnd_bar = miobnd+tol;

Q_vecl=[];
aux_constr1 = zeros(aux_num,n+k+aux_num);
aux_constr2 = zeros(aux_num,n+k+aux_num);
aux_constr3 = zeros(aux_num,n+k+aux_num);

s=0;
for i=1:n-1
    Q_vecl=[Q_vecl;Q(i+1:n,i)];
    
    % -e_i+x_ij <= 0
    aux_constr1(s+1:s+n-i,i)=-ones(n-i,1);
    aux_constr1(s+1:s+n-i,n+k+s+1:n+k+s+n-i)=eye(n-i);
    
    % -e_j+x_ij <= 0
    aux_constr2(s+1:s+n-i,i+1:n)=-eye(n-i);
    aux_constr2(s+1:s+n-i,n+k+s+1:n+k+s+n-i)=eye(n-i);
    
    % e_i+e_j-x_ij <= 1
    aux_constr3(s+1:s+n-i,i)=ones(n-i,1);
    aux_constr3(s+1:s+n-i,i+1:n)=eye(n-i);
    aux_constr3(s+1:s+n-i,n+k+s+1:n+k+s+n-i)=-eye(n-i);
    
    s=n-i+s;
end
 
 model.obj = [diag(Q)-2*Q*tau_vec;zeros(k,1);2*Q_vecl];
 model.A = sparse([diag(miobnd) -x zeros(n,aux_num); -diag(miobnd_bar) x zeros(n,aux_num);aux_constr1;aux_constr2;aux_constr3]);
 model.rhs = [miobnd*(1-tol)-y;y-tol*miobnd_bar;zeros(2*aux_num,1);ones(aux_num,1)];
 
else
   disp('error in input arguments');
   return;
end

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

try
      
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

