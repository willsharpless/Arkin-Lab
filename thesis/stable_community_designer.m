%% GA for designing communities with positive-valued steady states
% Created by Will Sharpless on October 17, 2019
% willsharpless@berkeley.edu

%% Import & Define Parameters

% Import Interaction coefficients and growth rates from Venturelli
S = xml2struct('/Users/willsharpless/Documents/MATLAB/arkin lab/thesis/twelveSpeciesT4.xml');

%% Reformat Ophelia's interaction/growth parameters
% This works from the xml2struct
param_struc_ca = S.sbml.model.listOfParameters.parameter;

id = strings(168,1);
value = zeros(168,1);
for i=1:168
    id(i) = param_struc_ca{i}.Attributes.id;
    value(i) = str2double(param_struc_ca{i}.Attributes.value);
end

id_table = array2table(cellstr(id),'VariableNames',{'id'}); value_table = array2table(value,'VariableNames',{'value'});
T = [id_table value_table];

mu_T = T(end-11:end,:);
alpha_T = T(1:end-24,:);

members = replace(string(table2cell(mu_T(:,'id'))),"u",""); %library of members
L = numel(members); %size of library

%% Randomly Choose Some Constituents

% Choose the desired size of the community (N <= L)
N = 11; 

% in future it might be more powerful to fluc the size (which is not
% important compared to containing interactions)
% could also remove useless members with parameter scan/sensitivity
% analyses

% Optimization parameters
S = 10; % string-guesses for genetic algorithm (each string is a set of organisms)
G = 3; % generations

%randomly choose members for community
[~, out] = sort(rand(S,L),2);
lambda = out(:,1:N); % all strings (S x N)

%should experiment with choosing unique sets of organisms via other
%functions

%each row of u correspond to a community's growth rate
mu = mu_T.value(lambda); %mu vectors for all strings (S x N)

A = zeros(N,N,S); % alpha matrix corresponding to lambda (N x N x S)

for i = 1:S
    Ni = 1;
    for ii = lambda(i,:)
        alphas = alpha_T.value(contains(alpha_T.id,["a"+members(ii)])); % one member to all others
        A(Ni,:,i) = alphas(lambda(i,:));
        Ni = Ni+1;
    end
end

% A(:,:,i), interaction mat for community i corresponding to row i of lambda:
%       _x1__x2_...__xn
%   x1 |a11 a12 ... a1n|
%   x2 |a21 a22 ... a2n|
%   ...|... ... ... ...|
%   xn |an1 ... ... ann|

%% Define the Jacobian
% Depends on the size of the community

syms x [N 1]
syms u [N 1]
syms a [N N]
syms gamma

u=0;
% x'.*(u'+a*x')
dX(:,1) = x.*(u+a*x);
C = gamma*x + u;


J = matlabFunction(jacobian(dX(:,1),x),'Vars',{x, a, u});
% returns function J(x,a,u)
% x - [N,1]
% a - [N,N]
% u - [N,1]

%% Influence Factor
I = cell(L,1); % Data of I's

% for i = 1:L
%     alphas_i_id = alpha_T.id(contains(alpha_T.id,members(i)));
%     alphas_i_val = alpha_T.value(contains(alpha_T.id,members(i)));
%     alphas_ij = alphas_i_val(contains(alphas_i_id,"a"+members(i))); % j intrxn on i
%     alphas_ji = alphas_i_val(~contains(alphas_i_id,"a"+members(i))); % i intrxn on j
%     mu_i = mu_T.value(contains(mu_T.id,members(i)));
%     
%     %how to organize the alphas correctly
%     Ifn = @(x) abs(mu_i*x(1)+sum(x.*alphas_ij)) + sum(abs(x.*alphas
%     
% end


%% Optimation and Simulation % UNDER CONSTRUCTION ----------------------------------

for g = 1:G %ga generation loop
    
 % reset values?
  
for s=1:S % string loop
 
 As = A(:,:,s); mus = mu(s,:)';
 fp = As\(-mus); % corresponds to species lambda(i,:)'
 
 eigs = eig(J(fp,As,mus));
 
 
 % plot gridspace of points within 1,1,...n
 % solve dx/dt for all points in space and store those which score under
 % threshold
 
 % calculate eig of jacobian for those met stable points NOTE SLOW
 % CALCULATION, IS IT NECESSARY?
 
 % score based on number of points/regions with dx/dt = 0 (meta-stable)
 % - I think this might be done by taking most minimal or all points below norm(grad)<0.1
 % - hmmmmm maybe I should cluster these regions with kmeans, NNMF or other thing Dr. Marshall mentioned in talk
 
 % clustering would allow us to detect the number of stable points as well
 % as their distance
 
 % score based on number of points/regions(clusters) with neg eig jacobians (stable)
 % score based on which of those points have the most non-zero constituent
 % populations
 
 % values/farthers from zero values
  
end %end string loop

 % sort parents
 % track score
 % breed
 % repeat

end %generation loop