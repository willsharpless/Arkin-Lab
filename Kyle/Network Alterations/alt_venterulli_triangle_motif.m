%% Alteration of Ophelia's Community In Silico
% willsharpless@berkeley.edu 
% July 27, 2020

u_a_og = readtable('/Users/willsharpless/Documents/MATLAB/arkin lab/Kyle/Network Alterations/venturelli2018_parameters.csv');
a_og = readtable('/Users/willsharpless/Documents/MATLAB/arkin lab/Kyle/Network Alterations/venturelli2018_original_adj_matrix.csv');
a_alt1 = readtable('/Users/willsharpless/Documents/MATLAB/arkin lab/Kyle/Network Alterations/venturelli2018_modified_adj_matrix_1.csv');
a_alt2 = readtable('/Users/willsharpless/Documents/MATLAB/arkin lab/Kyle/Network Alterations/venturelli2018_modified_adj_matrix_2.csv');
a_alt3 = readtable('/Users/willsharpless/Documents/MATLAB/arkin lab/Kyle/Network Alterations/venturelli2018_modified_adj_matrix_3.csv');
L = 12;

mu = u_a_og{end-L+1:end,2};
A(:,:,1) = a_og{:,2:end};
A(:,:,2) = a_alt1{:,2:end};
A(:,:,3) = a_alt2{:,2:end};
A(:,:,4) = a_alt3{:,2:end};

%%  Static Community Graphs
L=12;mx_a=max(abs(A(:)));pops = ones(1,L);
spec_names = strip(u_a_og.Var1(end-L+1:end),'left','u');
colors = ones(L,3);slfloop = 1;

for i=1:1
alphas = reshape(A(:,:,i),[1 L^2]) + 0.0001; %have to make nonzero values
alphas_d = abs(alphas - (max(alphas)+0.0001));
figure
% Comm_Grapher(L,pops,alphas_d,spec_names,mx_a,colors);
Comm_Grapher(L,pops,alphas,spec_names,mx_a,colors,slfloop);
end

%% Simulations

x0 = 0.1*ones(1,12); %get from ophelias pub
tspan = linspace(0,120);
xhat = zeros(100,12,4);

% Simulations
for i = 1:4
[~,xhat_i] = ode45(@(t,xhat) GLV_plant(t,xhat,mu,A(:,:,i)), tspan, x0);
% negative values are evaluated as 0 in the ode but need to be corrected
xhat_i(xhat_i<0) = 0;
xhat(:,:,i) = xhat_i;
end



%% Plotting Community Grpahs and Trajectories
set(0,'defaultTextInterpreter','latex')

% 12 Color Palette for the 12 species
p = [146, 43, 33; 203, 67, 53; 155, 89, 182; 125, 60, 152; 26, 82, 118; 52, 152, 219; 26, 188, 156; 22, 160, 133; 39, 174, 96; 241, 196, 15; 243, 156, 18; 211, 84, 0]./256;  
p = p(randperm(size(p, 1)), :);

filename = 'Vent_Comm_Alt_allalterations_120hrs.gif';
nImages = length(tspan);
% numsubplots = 4;
fig = figure('color','w','Units','normalized','OuterPosition', [0 0 1 1]);
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% plot at each time point for gif
for ti = 1:length(tspan)

t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% plot at each time point
for i=1:4 %numsubplots
    % Commgraph plot (left) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    set(gca, 'ColorOrder', p, 'NextPlot', 'replacechildren');
    set(gca,'XColor', 'none','YColor','none');
    
    pops = xhat(ti,:,i);
    Ai_weighted = A(:,:,i).*pops.*pops'; % magnify the edge by both giver and receiver size (aij*xi(t)*xj(t))
    alphas = reshape(Ai_weighted*40,[1 L^2]) + 0.0001;
    
    % need to highlight the motifs still
    cg = Comm_Grapher(L,abs(22*pops),alphas,spec_names,mx_a,p,slfloop);
    
%     if i==1
%         title({'Original'}, 'FontSize', 10,'Interpreter','Latex')
%     else
%         title(strcat({'Alteration '},num2str(i-1)), 'FontSize', 13,'Interpreter','Latex')
%     end
    
    % Trajectories plot (right) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    set(gca, 'ColorOrder', p, 'NextPlot', 'replacechildren');
    grid on
    plot(tspan,xhat(:,:,i),'-','LineWidth',2.5)
    xline(tspan(ti));
    
    if i==1
        title({'Original'}, 'FontSize', 13)
    else
        title(strcat({'Alteration '},num2str(i-1)), 'FontSize', 13)
    end
    xlabel('Time (hours)', 'FontSize', 10)
    ylabel('Abundnace', 'FontSize', 10)
    xlim([0 tspan(end)])
    ylim([0 1])
    delete(findobj('type','legend'))

end % plot for each ti

% need to include legend
% [abu.Properties.RowNames; strcat({'Sim-'}, abu.Properties.RowNames)]
legend([spec_names], ...
    'Position', [0.02 0.38 0.10 0.2], ...
    'Location', 'northwest', ...
    'FontSize', 12, 'NumColumns', 2)

% suptitle(strcat({'\bf{Community Graphs and Trajectories at Time '},num2str(ti),{'}'}))
title(t,strcat({'\bf{Community Graphs and Trajectories at Time '},num2str(ti),{'}'}),'FontSize', 20)

% gif stuff
drawnow
frame = getframe(fig);
im = frame2im(frame);
[imind, map] = rgb2ind(im,256);    
if ti == 1
    imwrite(imind, map,filename,'gif','LoopCount',Inf,'DelayTime',0.075);
else
    imwrite(imind, map,filename,'gif','WriteMode','append','DelayTime',0.075);
end

end % full gif iter
