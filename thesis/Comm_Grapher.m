function [cgplot] = Comm_Grapher(n,pops,alphas,spec_names,maxabsalpha,colors,selfloops)
% Community Graph Visualizer
% willsharpless@berkeley.edu

% n = 6; %num of community members
% pops = 50*rand(1,n); %size of the nodes (populations x_i)
% alphas = 2*rand(1,n^2)-1; %weights of the graph (alpha_ij)
% spec_names = {'Ab', 'Dc', 'Ef', 'Gh', 'Ij', 'Kl'}; %species names in cell

out = []; in = [];

for i=1:n
    out = [out i*ones(1,n)];
    in = [in 1:n];
end

if selfloops == 0
    dg = digraph(out, in, alphas,'omitselfloops'); %directed graph object
    idx = ones(1,n^2);
    for i=1:n
        idx(n*(i-1)+i) = 0;
    end
    alphas(~idx) = [];
    c = plot(dg,'Layout','circle', 'NodeColor', colors, 'NodeLabel',spec_names);
    c.EdgeCData=alphas;
    c.LineWidth=2.5*abs(alphas);
else
    dg = digraph(out, in, alphas); %directed graph object
    c = plot(dg,'Layout','circle', 'NodeColor', colors, 'NodeLabel',spec_names);
    c.EdgeCData=alphas;
    c.LineWidth=2.5*abs(alphas);
end

for i=1:n
   highlight(c, i,'MarkerSize',pops(i)) %rand(1,3)
end

%Defining a red to white to green colormap
dred = [182 0 0];
lred = [233 65 65];
white = [255 255 255];
lgreen = [88 189 98];
dgreen = [42 129 51];
raw = [dgreen; lgreen; white; lred; dred];
vec = linspace(0,100,5);
map = interp1(vec,raw,linspace(100,0,128),'pchip');


colormap(gca,map./255) 
% c.NodeLabel={};
caxis([-maxabsalpha maxabsalpha]) % I fixed this to max and min of ophelias data
% colorbar

%still need to reorient the Spec Names

cgplot = c;

