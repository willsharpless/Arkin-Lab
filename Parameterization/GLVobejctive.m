function cost = GLVobejctive(mu_A_flat,nzp,n,s,plot_bool,k)
% Objective Function for Parameterizing a GLV system
% willsharpless@berkeley.edu

% Inputs:
% mu_A_flat = [mu(:);A(:)];
%   mu - n by 1 vector of positive growth rates
%   A - n by b matrix of intrxns (diagonal must be negative!)
% nzp - noise percentage to inhibit overfitting
% n - size of the community
% s - organism list, s = [1 4:5 8:10 7] parameterizes comm of 1,4,5,8,9,10,7
% k - number of subjects to train parameterization on

mu = mu_A_flat(1:n);
A_flat = mu_A_flat(n+1:end);
A = reshape(A_flat,n,n);

% % FOR PARAMETERIZING MDSINE DATA
mdsine_data = load('mdsine_data.mat');
Abu_unscaled = mdsine_data.mdsine_data{1};
tspan = mdsine_data.mdsine_data{2};

% Community Reduction to Specified List s
Abu_unscaled = Abu_unscaled(s,:);

%bionumber conversion for numerical stability
scaled = Abu_unscaled{:,:}./1e14; %~absolute units
Abu = Abu_unscaled;
Abu{:,:} = scaled;

% % Apply noise to raw subject data (old method)
% noise = 0.01.*nzp.*Abu{:,:}.*(2*rand(size(Abu))-1); % percentage of each data * number in [-1,1]
% Abu{:,:} = Abu{:,:} + noise;

% Individual Mats for the Subject
i = 13; Abu1 = Abu(:,1:13); Abu2 = Abu(:,i+1:2*i); 
Abu3 = Abu(:,2*i+1:3*i); Abu4 = Abu(:,3*i+1:4*i); Abu5 = Abu(:,4*i+1:5*i);
tabs = {Abu1, Abu2, Abu3, Abu4, Abu5};
% 
% % Relative Abundance fitting hmmmmmmm
% load('mdsine_Rel_Abu.mat');
% Rel_abu = Rel_abu(s,:);
% i=13; Rel_abu1 = Rel_abu(:,1:13); Rel_abu2 = Rel_abu(:,i+1:2*i); 
% Rel_abu3 = Rel_abu(:,2*i+1:3*i); Rel_abu4 = Rel_abu(:,3*i+1:4*i); Rel_abu5 = Rel_abu(:,4*i+1:5*i);
% tabs = {Rel_abu1, Rel_abu2, Rel_abu3, Rel_abu4, Rel_abu5};

% FOR PARAMETERIZING OPHELIA's DATA
% load('growth_L1O7.mat');
% od_all = [s0 s1 s2]; % matrix of ophelia's od data for (141 time points (every half hour), 13 different communities {'ER-'},{'BH-'},{'EL-'},{'FP-'},{'DP-'},{'PC-'},{'CA-'},{'CH-'},{'BO-'},{'BV-'},{'BT-'},{'BU-'},{'NONE'})
% tspan = [0:12:72]'; %ophelia took rel_abu measurements every 12 hours for 72 hours
% adjusted_tspan = 1+[2*tspan(1:end-1);140]; %grabbing tp where rel_abu was measured (ophelia's last rel_abu tp is not in the OD so I am just using the 141's (while the 144th is correct)
% od = od_all(13,adjusted_tspan); % od of the full community
% load('vent_relabu.mat');
% vent_abu = vent_relabu.*od;

sum_norm_error = zeros(1,k);

% flag = 0; % for catching infinite cases

% for each subject, simulate from their x0 and score RMSE
for i=1:k
    
    % MDSINE
    abu = tabs{i};
    x0 = abu{:,1};
    xx = abu{:,:}';

    % Ophelia
%     abu = vent_abu;
%     x0 = abu(:,1);
%     xx = vent_abu';
    
    [~,xhat] = ode45(@(t,xhat) GLV_plant(t,xhat,mu,A), tspan, x0);
    
    % negative values are evaluated as 0 in the ode but need to be corrected
    xhat(xhat<0) = 0;
    
    %scoring
    if size(xhat,1) ~= size(xx,1) % can't simulate because params break ode45 (usually because going to inf)
        sum_norm_error(i) = 10000; %explode score to avoid
    else
        % implementing noise band score protection (scores within nzp percent ring of data are not penalized)
        xx_minus_xhat = abs(xx-xhat) - abs(xx).*(0.01*nzp);
        xx_minus_xhat(xx_minus_xhat < 0) = 0;
        sum_norm_error(i) = sum(sqrt(sum((xx_minus_xhat).^2,2)));
    end  
    % Might be worth trying the mean norm instead...    
end

% Cost is chosen to be the average sum norm error over all subjects to
cost = mean(sum_norm_error);

% % Cost is chosen to be the sum_norm over all subjects
% cost = sum(sum_norm_error);

%% Plotting (NOT USED IN OPTIMIZATION)

if plot_bool == 1

    % 12 Color Palette for the 12 species
    p = [146, 43, 33; 203, 67, 53; 155, 89, 182; 125, 60, 152; 26, 82, 118; 52, 152, 219; 26, 188, 156; 22, 160, 133; 39, 174, 96; 241, 196, 15; 243, 156, 18; 211, 84, 0]./256;  
    p = p(randperm(size(p, 1)), :);
    figure('color','w','Units','normalized','OuterPosition', [0 0 1 1]);
    % tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for i=1:5
        subplot(2,3,i);
    %     subaxis(2,3,i, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
    %     nexttile
        set(gca, 'ColorOrder', p, 'NextPlot', 'replacechildren');

        abu = tabs{i};
        x0 = abu{:,1};
        [tt, xx] = ode45(@(t,x) GLV_plant(t,x,mu,A), [0.75 28], x0);
        grid on
        plot(tspan,abu{:,:}',':','LineWidth',2.5)
        hold on
        plot(tt,xx,'-','LineWidth',2.5)
        hold off
        
        if i<=k
            title(strcat({'Trained: Subject - '},sprintf('%.1d',i)), 'FontSize', 15)
        else
            %cacualte pearson's for multi d
            [~,xhat] = ode45(@(t,xhat) GLV_plant(t,xhat,mu,A), tspan, x0);
            xbar = mean(abu{:,:}');
            r2 = 1 - (sum(sum((abu{:,:}'-xhat).^2,2))/sum(sum((abu{:,:}'-xbar).^2,2))); % 1 - SS_res/SS_tot (using 2-norm not diff)
            title(strcat({'Predicted: Subject - '},num2str(i)), 'FontSize', 15) %,{', r^2 = '},num2str(r2)
        end
        
        xlabel('Time (days)', 'FontSize', 10)
        ylabel('Abundnace', 'FontSize', 10)
        xlim([0 tspan(end)])
        ylim([min(Abu{:,:},[],'all') max(Abu{:,:},[],'all')])
        delete(findobj('type','legend'))
        
    end

    legend([abu.Properties.RowNames; strcat({'Sim-'}, abu.Properties.RowNames)], ...
        'Position', [0.7088 0.2 0.1778 0.1957], ...
        'FontSize', 12, 'NumColumns', 2)

    suptitle(strcat({'Trained on Subjects 1 to '},num2str(k),{', with '},num2str(nzp),{'% noise added, ST3'}))
end

