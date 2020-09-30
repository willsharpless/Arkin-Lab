%% Reduced Modeling attempts
% Includes visualization with Comm_Grpaher

% Author: Will Sharpless, Feb 27, 20
% Last edited: May 28, 2020

% The first half of the code was written for Ophelia's data and parameters from her paper (Deciphering microbial comm....)
% incomplete without her unpublished OD data

% The second half of the code was written to use the data and parameters from the 2016 MDSINE paper

% Parameters of the MDSINE data were unable to predict the dynamics of the
% community for some reason ...


%% Import & Reformat Ophelia's interaction/growth parameters

% % Import Interaction coefficients and growth rates from Venturelli
% S = xml2struct('/Users/willsharpless/Documents/MATLAB/arkin lab/thesis/twelveSpeciesT4.xml');

% This works from the xml2struct
% param_struc_ca = S.sbml.model.listOfParameters.parameter;
% 
% id = strings(168,1);
% value = zeros(168,1);
% for i=1:168
%     id(i) = param_struc_ca{i}.Attributes.id;
%     value(i) = str2double(param_struc_ca{i}.Attributes.value);
% end
% 
% id_table = array2table(cellstr(id),'VariableNames',{'id'}); 
% value_table = array2table(value,'VariableNames',{'value'});
% T = [id_table value_table];
% 
% mu_T = T(end-11:end,:);
% alpha_T = T(1:end-24,:);
% 
% members = replace(string(table2cell(mu_T(:,'id'))),"u",""); %library of members
% L = numel(members); %size of library

%% Importing MDSINE CDIFF (first 13 timepoints) data

S = readtable('/Users/willsharpless/Documents/MATLAB/arkin lab/MDSINE-source-code-v1.2/output_cdiff_13/BVS.results.parameters_13.txt');
        
mu_T_pre = S(strcmp(S.parameter_type,{'growth_rate'}),:);
alpha_T_pre = S(strcmp(S.parameter_type,{'interaction'}),:);

mu_id = cellfun(@(s) strcat('u',s(1),s(strfind(s,'-')+1)), mu_T_pre.target_taxon, 'UniformOutput', false);
abbrev_id = cellfun(@(s) strcat(s(1),s(strfind(s,'-')+1)), mu_T_pre.target_taxon, 'UniformOutput', false);
alpha_id_source = cellfun(@(s) strcat(s(1),s(strfind(s,'-')+1)), alpha_T_pre.source_taxon, 'UniformOutput', false);
alpha_id_target = cellfun(@(s) strcat(s(1),s(strfind(s,'-')+1)), alpha_T_pre.target_taxon, 'UniformOutput', false);
alpha_id = strcat('a',alpha_id_source,alpha_id_target);

mu_T = [array2table(mu_id,'VariableNames',{'id'}) mu_T_pre(:,{'value'})];
alpha_T = [array2table(alpha_id,'VariableNames',{'id'}) alpha_T_pre(:,{'value'})];

members = replace(string(table2cell(mu_T(:,'id'))),"u",""); %library of members
L = numel(members); %size of library

%% Create MU and A matricies for the DropOut Communities
% 
% % Choose the desired size of the community (N <= L)
% N = 11; 
% S = factorial(L)/(factorial(N)*factorial(L-N)); % qty of all possible communities, useful if N =/= 11
% 
% out = 1:L; temp1 = [out; zeros(S-1,L)];
% for i=1:(S-1)
%     temp2 = circshift(out,[0,-i]);
%     temp1(i+1,:) = temp2;
% end
% 
% lambda = temp1(:,1:N); % all strings (S x N)
% 
% %each row of u corresponds to a community's growth rate
% mu = mu_T.value(lambda); %mu vectors for all strings (S x N)
% muid = mu_T.id(lambda);
% 
% A = zeros(N,N,S);% alpha matrix corresponding to lambda (N x N x S)
% Aid = cell(N,N,S);
% 
% for i = 1:S
%     Ni = 1;
%     for ii = lambda(i,:)
%         alphas = alpha_T.value(contains(alpha_T.id,["a"+members(ii)])); % one member to all others
%         alphasid = alpha_T.id(contains(alpha_T.id,["a"+members(ii)]));
%         A(Ni,:,i) = alphas(lambda(i,:));
%         Aid(Ni,:,i) = alphasid(lambda(i,:));
%         Ni = Ni+1;
%     end
% end
% 
% % A(:,:,i), interaction mat for community i corresponding to row i of lambda and mu:
% %       _x1__x2_...__xn
% %   x1 |a11 a12 ... a1n|
% %   x2 |a21 a22 ... a2n|
% %   ...|... ... ... ...|
% %   xn |an1 ... ... ann|

%% Creat MU and A for the full community 
N=12;
lambda = 1:N;

mu12 = mu_T.value(lambda); %mu vectors for all strings (S x N)
muid12 = mu_T.id(lambda);

A12 = zeros(N,N);% alpha matrix corresponding to lambda (N x N x S)
Aid12 = cell(N,N);

for i = lambda
    A12(i,:) = alpha_T.value(contains(alpha_T.id,["a"+members(i)])); % one member to all others
    Aid12(i,:) = alpha_T.id(contains(alpha_T.id,["a"+members(i)]));
end

%% Community Graph (Test)

spec_names = strip(mu_T.id(lambda(1,:)),'left','u');
alphas = reshape(A12(:,:,1),[1 L^2]) + 0.0001; %have to make nonzero values
Comm_Grapher(L,ones(1,L),alphas,spec_names,max(abs(A12(:))));

%% Loading dataset EV3, drop out and full community

% rows correspond to t1 thru t7 and columns correspond to species 1.BH, 2.CA, 3.BU, 4.PC, 5.BO, 6.BV, 7.BT, 8.EL, 9.FP, 10.CH, 11.DP, 12.ER
data12_badorder = [0.0138    0.0080    0.0642    0.1774    0.1252    0.0468    0.1004    0.0208    0.1943    0.0083    0.0530    0.1878;
    0.0097    0.0032    0.1913    0.0091    0.2654    0.0340    0.0375    0.0329    0.1663    0.0953    0.0306    0.1248;
    0.0077    0.0126    0.2257    0.0045    0.3779    0.0307    0.0446    0.0384    0.0879    0.0656    0.0638    0.0406;
    0.0129    0.0078    0.2567    0.0002    0.3713    0.0055    0.0868    0.0290    0.1338    0.0642    0.0041    0.0277;
    0.0090    0.0247    0.2528    0.0004    0.3927    0.0214    0.0784    0.0204    0.1536    0.0212    0.0109    0.0145;
    0.0071    0.0403    0.3763    0.0001    0.2144    0.0079    0.1439    0.0246    0.1674    0.0099    0.0004    0.0077;
    0.0100    0.0249    0.3649    0.0000    0.1533    0.0140    0.1452    0.0168    0.2393    0.0212    0.0019    0.0086];

% column order now matches mu and A order: 1.BH, 5.BO, 7.BT, 3.BU, 6.BV, 2.CA, 10.CH, 11.DP, 8.EL, 12.ER, 9.FP, 4.PC
vent_relabu = [data12_badorder(:,1) data12_badorder(:,5) data12_badorder(:,7) data12_badorder(:,3) data12_badorder(:,6) ...
    data12_badorder(:,2) data12_badorder(:,10) data12_badorder(:,11) data12_badorder(:,8) data12_badorder(:,12) ...
    data12_badorder(:,9) data12_badorder(:,4)]';

save('vent_relabu.mat','vent_relabu')

% need to multiply by OD at each time point
%data12_abs = data12.*[od1;od2 ...; od7]

%% Loading Raw Biomass and Counts Data from MDSINE cdiff

Bm = readtable('/Users/willsharpless/Documents/MATLAB/arkin lab/MDSINE-source-code-v1.2/data_cdiff/biomass_13_cutout.txt');
Cts = readtable('/Users/willsharpless/Documents/MATLAB/arkin lab/MDSINE-source-code-v1.2/data_cdiff/counts_13_cutout.txt', 'HeaderLines', 1);

Bm.avg = mean(Bm{:,:},2);

Cts = Cts(ismember(Cts.Var1,mu_T_pre.target_taxon),:); %species MDSINE gave params for (significant counts throughout the experiment)
Cts.Properties.RowNames = Cts.Var1;
Cts.Var1 = [];

% Abundance table is counts * biomass
Abu = Cts;
Abu{:,1:end} = Cts{:,1:end}.*Bm.avg';

Rel_abu = Cts;
total_cts_ti = sum(Cts{:,:},1);
Rel_abu{:,:} = Cts{:,:}./total_cts_ti;
Rel_abu.Avg = mean(Rel_abu{:,:},2);

% Each mouse subject abundance
i = 13;
Abu1 = Abu(:,1:13);
Abu2 = Abu(:,i+1:2*i);
Abu3 = Abu(:,2*i+1:3*i);
Abu4 = Abu(:,3*i+1:4*i);
Abu5 = Abu(:,4*i+1:5*i);

Rel_abu1 = Rel_abu(:,1:13);
Rel_abu2 = Rel_abu(:,i+1:2*i);
Rel_abu3 = Rel_abu(:,2*i+1:3*i);
Rel_abu4 = Rel_abu(:,3*i+1:4*i);
Rel_abu5 = Rel_abu(:,4*i+1:5*i);

Cts1 = Cts(:,1:13);
Cts2 = Cts(:,i+1:2*i);
Cts3 = Cts(:,2*i+1:3*i);
Cts4 = Cts(:,3*i+1:4*i);
Cts5 = Cts(:,4*i+1:5*i);

%time vector
meta = readtable('/Users/willsharpless/Documents/MATLAB/arkin lab/MDSINE-source-code-v1.2/data_cdiff/metadata_13.txt');
tt = meta.measurementid(1:13);

mdsine_data = {Abu, tt};
  
%% Writing a file
%save('mdsine_data.mat','mdsine_data');

%% Plotting Raw Biomass for Mouse Subjects

figure
hold on
for i=1:5
    plot(1:13,Bm.avg((i-1)*13+1:i*13))
end

%% Plotting Abundance for Mouse Subjects
tabs = {Abu1, Abu2, Abu3, Abu4, Abu5};
tabs = {Cts1, Cts2, Cts3, Cts4, Cts5};
tabs = {Rel_abu1, Rel_abu2, Rel_abu3, Rel_abu4, Rel_abu5};


% for i=1:5
i=1;

%     abu = tabs{i};
%     x0 = abu{:,1};
    x0 = rand(1,N);
    [t, x] = ode45(@(t,x) GLV_plant(t,x,mu12,A12), [0.75 100], x0);
    
    figure('color','w')
    grid on
%     plot(tt,abu{:,:}')
    hold on
    plot(t,x,'*')
    hold off
    
%     legend([abu.Properties.RowNames; strcat({'Sim-'},abu.Properties.RowNames)],'Location','North')
%     legend(abu.Properties.RowNames,'Location','North')
    title(strcat('Subject ',sprintf('%.1d',i)))
    xlabel('Time (days)')
    ylabel('Abundnace')
%     xlim([0 tt(end)])
%     ylim([min(Abu{:,:},[],'all') max(Abu{:,:},[],'all')])
    
% end

%% Reduced Modeling

%Iterate through each of the 12 dropout communities and the full (start with one for debug)
%for i =1:S
% i = 1;

%load data into temps for ease
% mui = mu(i,:);
% muidi = muid(i,:);
% Ai = A(:,:,i);
% Aidi = Aid(:,:,i);

% %load data
% mui = mu12(i,:);
% muidi = muid12(i,:);
% Ai = A12(:,:,i);
% Aidi = Aid12(:,:,i);
% tdf = tt;

% With Ophelia's dropout data and the parameterized interactions
% infl_avg = zeros(length(muid),1);
% for i=1:length(mui)
%     infltemp = zeros(length(muid),1);
%     for t=1:7
%         infltemp(t) = tdf(t,i)*(abs(mui(i)));
%     end
%     infl_avg(i) = mean(infltemp);
% end
% calculate will's influence factor from abundances @ each time point then average across all except the first (OR simply at 72 hr)
% calculate ophelia's impact score from interactions and 72 hr relative abundances (as she did)

% % Model Full community first
% x0 = 
% [t, x] = ode23(@(t,y) GLV_plant(t,x,mui,Ai), [0.75 28], y0tt(1));

varnames = {'Model_id',abbrev_id{:},'Size','RMSE','Rel_Abu_Sum'};
% rown = 4096 - 66 - 12 - 1; % sum of all combinations of 12 = 2^12, we do not care about n=2,n=1,n=0 case
rown = 220; %for debugging

% All data will be stored in the "Reduced" Data Frame
Red_df = array2table(zeros(rown,length(varnames)),'VariableNames',varnames);
Red_df.Model_id = [1:rown]';
Red_df.Traj = cell(rown,1);

c = 0; % idx counter
%iterate through through size of reduced community
for i = 3 %:12
    C = combnk(1:12,i);
    
    %iterate through which x of the community to model in reduced community
    for ii = 1:length(C)
        
        %Discern and Detail Redcomm
        c = c+1;
        set = C(ii,:); %numerical set of orgs
        Red_df(c,abbrev_id(set)) = {1}; %storing in Red_df
        Red_df(c,{'Size'}) = {i};
        
        %Compute the sum Avg Relative Abunance for Redcomm
        Red_df(c,{'Rel_Abu_Sum'})={sum(Rel_abu.Avg(set))};
        
        %Locate Redcomm Parameters
        muii_T = mu_T(set,:);
        alphaii = [];
        for iii=set
            alphaii = [alphaii ((iii-1)*L + set)];
        end
        alphaii_T = alpha_T(alphaii);
        
        %Model
        [t, y] = ode23(@(t,y) GLV_plant(t,y,mui,Ai), [0.75 28], t0);
        
        %store trajectory!
        
        
        %Calculate the RMSE for appropriate pops
        
        
        %might involve interpolation
        
        
    end
    
end


% model with initial condition, 









%end
%


