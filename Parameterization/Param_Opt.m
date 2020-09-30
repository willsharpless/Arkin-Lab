%%% GLV Parameter Optimization of MDSINE data
%willsharpless@berkeley.edu

clc
clear all

numloops = 5; %number of times to run fmincon on subcommunities before taking the optimum scoring (random initial)

% Optimization Hyperparameters
Aopt = []; Aopt_eq = [];
bopt = []; bopt_eq = [];
options = optimoptions('fmincon','MaxFunctionEvaluations',Inf,'MaxIterations',Inf); % ,'Display','iter'


% Optimizing over hyper parameters: noise percentage and num of training subjects %%%%%%%%%%%%%%%%%%%%%%%%%

for nzp = [0 1 5]
for k = 1:3
nzp = 0;
k = 1;
    
disp(strcat({'Parameterizing by training on '},num2str(k),{' subject(s) with '},num2str(nzp),{'% noise bands'}))

% Primary Optimization of 3 Members Subcommunities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 3;
I = eye(n);

% Storing best scoring parameters
Organisms = cell(4*numloops,1);
Mu_A_flat = cell(4*numloops,1);
Score = zeros(4*numloops,1);
Subcomm_scores = table(Organisms, Mu_A_flat, Score);

Org_mat = [1 5 6; 2 3 4; 7 8 9; 10 11 12]; %already done 1/5/6 so a little messy

for j = 1:4

% Three member groups of species to ease burden of fmincon
s = Org_mat(j,:);

% Looping over different initial guess to take the highest scoring
for i = 1:numloops

% Initial guesses and bounds
mu_0 = rand(1,n)./100; %?? [0 0.01]
A_0 = (rand(n)-0.5)./4; %?? [-0.25 0.25]
A_0(I>0) = abs(A_0(I>0)).*-1; % need A_0(i,i) <=0
mu0_A0_flat = [mu_0(:);A_0(:)];
mu_lb = zeros(n,1);
mu_ub = 5*ones(n,1);
A_lb = -5*ones(n);
A_ub = 5*ones(n).*(I<1); % need A(i,i) <=0
lb = [mu_lb(:);A_lb(:)];
ub = [mu_ub(:);A_ub(:)];

% k = 1; % train on only 1 subject
plot_bool = 0;
fun = @(x) GLVobejctive(x,nzp,n,s,plot_bool,k);

[pX,fval] = fmincon(fun, mu0_A0_flat,...
    Aopt,bopt,Aopt_eq,bopt_eq,lb,ub,[],options); %,exitflag,output

% Storing Data
Subcomm_scores.Organisms{i+(j-1)*numloops} = s;
Subcomm_scores.Mu_A_flat{i+(j-1)*numloops} = pX;
Subcomm_scores.Score(i+(j-1)*numloops) = fval;

end % subcomm param looping for best
end % for all subcomms

% save('St.mat','St')

% Parameterization of Subject 1 full community combining 3-member optimizations %%%%%%%%%%%%%%%%%%%%%%%%

% build educated guess on random guess format from above
n = 12;
I = eye(n);
good_mu0 = rand(1,n)./100; %?? [0 0.01]
good_A0 = (rand(n)-0.5)./4; %?? [-0.25 0.25]
good_A0(I>0) = abs(good_A0(I>0)).*-1; % need A_0(i,i) <=0

% good_A0 = zeros(12,12); %debugging

% input subcomm optimized values to educated guess
for j=1:4
    [~, idx] = min(Subcomm_scores.Score(((j-1)*numloops+1):(j*numloops))); %identify best scoring 3-mem opt
    idx = idx + (j-1)*numloops;
    params = Subcomm_scores.Mu_A_flat{idx};
    good_mu0(Subcomm_scores.Organisms{idx}) = params(1:3);
    good_A0(Subcomm_scores.Organisms{idx},Subcomm_scores.Organisms{idx}) = reshape(params(4:end),3,3);
end

good_mu0_A0_flat = [good_mu0(:);good_A0(:)];

% adjust bounds for full 12
mu_lb = zeros(n,1);
mu_ub = 5*ones(n,1);
A_lb = -5*ones(n);
A_ub = 5*ones(n).*(I<1); % need A(i,i) <=0
lb = [mu_lb(:);A_lb(:)];
ub = [mu_ub(:);A_ub(:)];

% full comm optimization
s = 1:12;
% nzp = 5; % 5% noise?
% k = 3; % train on only 3 subjects?
plot_bool = 0;
fun = @(x) GLVobejctive(x,nzp,n,s,plot_bool,k);

[pX_full_sj1,fval] = fmincon(fun, good_mu0_A0_flat,...
    Aopt,bopt,Aopt_eq,bopt_eq,lb,ub,[],options); % ,exitflag,output

plot_bool = 1;
GLVobejctive(pX_full_sj1,nzp,n,s,plot_bool,k) %plotting optimized solution

end % num training subjects
end % noise percentage

% save('pX_full_sj1.mat','pX_full_sj1')
% Fits to Subject 1 well! Mediocre prediction of other subjects though...

%% Training Parameters on Multiple Subjects and Noise, 
% This did not work, could not optimize much further after training on 1
% full community to then move to 3 communities, have to do it with
% subcommunities first

load 'pX_full_opt.mat'
pX_full_sj1 = pX_full;

n=12;
I=eye(12);

mu_lb = zeros(n,1);
mu_ub = 5*ones(n,1);
A_lb = -5*ones(n);
A_ub = 5*ones(n).*(I<1); % need A(i,i) <=0
lb = [mu_lb(:);A_lb(:)];
ub = [mu_ub(:);A_ub(:)];

% Train on Subjects 1 through 3
k = 3;
plot_bool = 0;
s=1:12;

fun = @(x) GLVobejctive(x,nzp,n,s,plot_bool,k);
[pX_full_allsj,fval] = fmincon(fun, pX_full_sj1,...
    Aopt,bopt,Aopt_eq,bopt_eq,lb,ub,[],options); %,exitflag,output

% Predict Trajectories of 4 and 5

plot_bool = 1;
GLVobejctive(pX_full_allsj,nzp,n,s,plot_bool,k)