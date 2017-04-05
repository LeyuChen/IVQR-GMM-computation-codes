
% Given (y,x), this function solves the following maximization problem
% for each i, max |y(i)-x(i,:)*b| over b confined to the space
% described by bnd

% function input :
% y     : (n by 1) matrix of outcomes
% x     : (n by k) matrix of covariate data    
% bnd   : (k by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients

% function output :
% value : the value of the maximized objective function

function value = miobnd_fn(y,x,bnd)

n=length(y);

model.modelsense = 'max';
model.sense = '>';
model.lb = bnd(:,1);
model.ub = bnd(:,2);

tol=1e-6;

params.outputflag = 0; 
params.OptimalityTol=tol;
params.FeasibilityTol=tol;
params.IntFeasTol=tol;

value=zeros(n,1);

for i=1:n

v=zeros(2,1);    
    
model.obj = -x(i,:);
model.objcon = y(i);

try
    model.A = sparse(-x(i,:));
    model.rhs = -y(i);
    result= gurobi(model, params);
    %disp(result.x);
    v(1)=result.objval;
    %rtime=result.runtime;
    %disp(v(1));
    %fprintf('Optimization returned status: %s\n', result.status);
  catch gurobiError
    %fprintf('Error reported\n');disp(v(1));
end

model.obj = x(i,:);
model.objcon =-y(i);

try
    model.A = sparse(x(i,:));
    model.rhs = y(i);
    result = gurobi(model, params);
    %disp(result.x);
    v(2)=result.objval;
    %rtime=result.runtime;
    %disp(v(2));
    %fprintf('Optimization returned status: %s\n', result.status);
  catch gurobiError
    %fprintf('Error reported\n');disp(v(2));
end
value(i)=max(v);
end
end

