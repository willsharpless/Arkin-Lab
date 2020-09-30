close all
clear all
clc

%INITIAL GUESS
load PW22_result.mat
Po = param6; 

k = 6;

%Parameter fitting 
laminp = logspace(-5,1,12);
laminp = laminp(k);

%UPPER AND LOWER BOUNDS
lbx = [0.001 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 -10];
lb = repmat(lbx,12,1)';

ubx = [2 -0.01 10 10 10 10 10 10 10 10 10 10 10];
ub = repmat(ubx,12,1)';
A = [];
b = [];
Aeq = [];
beq = [];

lam = laminp;

options = optimoptions('fmincon','MaxFunctionEvaluations',Inf,'MaxIterations',Inf);
tic
[pX,fval,exitflag,output] = fmincon(@(x)objectiveFunction_2(x,lam),Po(:),A,b,Aeq,beq,lb(:),ub(:),[],options);
%[pX,fval,exitflag,output] = fmincon(@(x)objectiveFunction(x,lam),Po(:),A,b,Aeq,beq,lb(:),ub(:),[]);
toc

paramoptM = pX(:);
paramregM = Po(:);



