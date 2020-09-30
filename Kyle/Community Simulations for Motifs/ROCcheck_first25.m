%%% Native Community Simulation - Radius of Convergence Check and Desolve Validation
% Will Sharpless willshaprless@berkeley.edu
% June, 2020

%% Import the Params

df = readtable('/Users/willsharpless/Documents/MATLAB/arkin lab/Kyle/Community Simulations for Motifs/NativeCommunities_first25_alldata.xls');

mu = df(:,1:5);
alpha = df(:,6:30);
z_stode = df(:,31:35);
z_desolve = df(:,36:40);

df.fp = cell(size(df,1),1); %this will be the nontrivial fixed point -A\mu
df.fp_positive = zeros(size(df,1),1);
df.ROC = zeros(size(df,1),1);
df.stode_est_in_ROC = zeros(size(df,1),1);
df.desolve_est_in_ROC = zeros(size(df,1),1);
df.sto_des_pdist_stonorm = zeros(size(df,1),1);
df.sto_des_pdist_desnorm = zeros(size(df,1),1);
df.r_desolve_p = zeros(size(df,1),1);
df.sto_des_state_pdist = cell(size(df,1),1);

tic
for i=1:size(df,1)
    
    mui = mu{i,:}';
    alphai = alpha{i,:};
    Ai = reshape(alphai,5,5)';
    df.fp{i} = -Ai\mui;
    df.fp_positive(i) = all(-Ai\mui > 0);
    
    % the eig(J) are rly small, is it possible to maximize them without
    % opening up the interaction constraints?
    
    J_z = Ai.*df.fp{i};
    P_z = lyap(J_z,eye(5));
    spec_norm_P = sqrt(max(eig(P_z'*P_z))); % aka the induced two norm of P
    
    r_z = 1/(2*sqrt(5)*spec_norm_P*norm(Ai,'fro'));
    df.ROC(i) = r_z;
    C = min(eig(P_z))*r_z^2; % this is really fuckin small
    
    Vz_stode = z_stode{i,:}*P_z*z_stode{i,:}'; % this is like 10^6 x greater than C
    Vz_desolve = z_desolve{i,:}*P_z*z_desolve{i,:}'; % this is like 10^6 x greater than C
    
    r_stode = norm(df.fp{i}-z_stode{i,:}'); % both of these are within the frickin radius urg
    
    % percent distances of the three different points
    df.sto_des_pdist_stonorm(i) = norm(z_stode{i,:}' - z_desolve{i,:}')/norm(z_stode{i,:}');
    df.sto_des_pdist_desnorm(i) = norm(z_stode{i,:}' - z_desolve{i,:}')/norm(z_desolve{i,:}');
    df.r_desolve_p(i) = norm(df.fp{i}-z_desolve{i,:}')./norm(max(df.fp{i}));
    % sto and fp calculation are within 10^-5 of one another, desolve
    % unpredictably varies in distance (poor measure)
    
    % individual state measures of accuracy of desolve w/ stode
    df.sto_des_state_pdist{i} = abs(z_stode{i,:}' - z_desolve{i,:}')./sum(z_stode{i,:});
    % does this correlate with percentage biomass?
    df.des_pbiomass{i} = z_desolve{i,:}'./sum(z_desolve{i,:});
    df.corr(i) = corr2(df.sto_des_state_pdist{i},df.des_pbiomass{i});
    
    df.stode_est_in_ROC(i) = (Vz_stode <= C);
    df.desolve_est_in_ROC(i) = (Vz_desolve <= C); 
end
toc
%write table to csv