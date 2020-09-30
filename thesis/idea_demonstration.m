%% Calculation of the Jacobian
% 
% syms BU FP BO aBUBU aBUFP aBUBO aFPFP aFPBU aFPBO aBOBO aBOBU aBOFP uBU uFP uBO
% dBU = BU*(uBU + aBUBU*BU + aBUFP*FP + aBUBO*BO);
% dFP = FP*(uFP + aFPFP*FP + aFPBU*BU + aFPBO*BO);
% dBO = BO*(uBO + aBOBO*BO + aBOBU*BU + aBOFP*FP);
% % N = sqrt((dBU^2) + (dFP^2) + (dBO^2));
% % N_f = matlabFunction(N, 'vars', {'BU','FP','BO'});
% 
% % G = matlabFunction(gradient(N,[BU, FP, BO]));
% % H = matlabFunction(hessian(N,[BU, FP, BO]));
% J = matlabFunction(jacobian([dBU, dFP, dBO],[BU, FP, BO]));

%% Gradient Field
%%%%%%% DONT FORGET TO RENAME THE GIF BEFORE RUNNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Validating Perspective on Communities with Venturelli Data

uBU = 0.5841; %growth rate of Bacteroides uniformis
uFP = 0.2192; %growth rate of Faecalibacterium prausnitzii
uBO = 0.4777; %growth rate of Bacteroides ovatus

% inter and intra species interaction coefficients
aBUBU = -0.8804;
aBUFP = -0.7822;
aBUBO = -0.9208;

aFPFP = -1.0382;
aFPBU = 0.2313;
aFPBO = -0.0993;

aBOBO = -0.7339;
aBOBU = -0.6317;
aBOFP = -0.2087;

A = [aBUBU aBUFP aBUBO; aFPBU aFPFP aFPBO; aBOBU aBOFP aBOBO];
u=[uBU; uFP; uBO];

fig=figure('Position', [10 10 1200 800]);
set(gcf,'color','w');

% Uncomment to make a rotating gif
for zz=1:180
    
fp_real = A\(-u);

res_q = 8;

%Dynamics
% [BU, FP, BO] = meshgrid(resolution_quiver);
[BU, FP, BO] = meshgrid(0:0.55/res_q:0.55,0:0.4/res_q:0.4,0:0.7/res_q:0.7);

dBU = BU.*(uBU + aBUBU.*BU + aBUFP.*FP + aBUBO.*BO);
dFP = FP.*(uFP + aFPFP.*FP + aFPBU.*BU + aFPBO.*BO);
dBO = BO.*(uBO + aBOBO.*BO + aBOBU.*BU + aBOFP.*FP);

% Colored Quiver
q = quiver3(BU, FP, BO, dBU, dFP, dBO, 'LineWidth', 1.5,'AutoScaleFactor', 3,'HandleVisibility','off'); %, 'Marker', 'o'
colormap
hold on

norm = sqrt((dBU.^2) + (dFP.^2) + (dBO.^2));

% currentColormap = colormap(jet);

%Custom Color Map
custom_colormatrix = [255,255,49; 255,253,52; 254,250,55; 254,248,58; 253,245,61; 253,243,63; 252,241,66; 252,238,69; 251,236,72; 251,234,74; 250,231,77; 250,229,80; 249,226,83; 249,224,86; 249,222,88; 248,219,91; 248,217,94; 247,215,97; 247,212,99; 246,210,102; 246,207,105; 245,205,108; 245,203,111; 244,200,113; 244,198,116; 243,196,119; 243,193,122; 243,191,124; 242,188,127; 242,186,130; 241,184,133; 241,181,136; 240,179,138; 240,177,141; 239,174,144; 239,172,147; 238,169,149; 238,167,152; 237,165,155; 237,162,158; 237,160,160; 236,158,163; 236,155,166; 235,153,169; 235,150,172; 234,148,174; 234,146,177; 233,143,180; 233,141,183; 232,139,185; 232,136,188; 232,134,191; 231,131,194; 231,129,197; 230,127,199; 230,124,202; 229,122,205; 229,120,208; 228,117,210; 228,115,213; 227,112,216; 227,110,219; 226,108,222; 226,105,224; 225,104,227; 221,104,227; 218,105,227; 214,105,228; 211,106,228; 207,106,229; 204,107,229; 200,107,230; 196,108,230; 193,108,231; 189,109,231; 186,109,231; 182,110,232; 179,110,232; 175,111,233; 172,111,233; 168,112,234; 165,112,234; 161,113,235; 158,113,235; 154,114,235; 150,114,236; 147,115,236; 143,115,237; 140,116,237; 136,116,238; 133,117,238; 129,117,239; 126,118,239; 122,118,240; 119,119,240; 115,119,240; 111,120,241; 108,120,241; 104,121,242; 101,121,242; 97,122,243; 94,122,243; 90,123,244; 87,123,244; 83,124,244; 80,124,245; 76,125,245; 73,125,246; 69,126,246; 65,126,247; 62,127,247; 58,127,248; 55,128,248; 51,128,248; 48,129,249; 44,129,249; 41,130,250; 37,130,250; 34,131,251; 30,131,251; 27,132,252; 23,132,252; 19,132,253; 16,133,253; 12,133,253; 9,134,254; 5,134,254; 2,135,255; 0,134,253; 0,132,250; 0,130,246; 0,128,243; 0,126,240; 0,124,236; 0,121,233; 0,119,229; 0,117,226; 0,115,223; 0,113,219; 0,111,216; 0,109,212; 0,107,209; 0,104,205; 0,102,202; 0,100,199; 0,98,195; 0,96,192; 0,94,188; 0,92,185; 0,90,182; 0,88,178; 0,85,175; 0,83,171; 0,81,168; 0,79,164; 0,77,161; 0,75,158; 0,73,154; 0,71,151; 0,68,147; 0,66,144; 0,64,140; 0,62,137; 0,60,134; 0,58,130; 0,56,127; 0,54,123; 0,51,120; 0,49,117; 0,47,113; 0,45,110; 0,43,106; 0,41,103; 0,39,99; 0,37,96; 0,34,93; 0,32,89; 0,30,86; 0,28,82; 0,26,79; 0,24,76; 0,22,72; 0,20,69; 0,17,65; 0,15,62; 0,13,58; 0,11,55; 0,9,52; 0,7,48; 0,5,45; 0,3,41; 0,1,38; 0,0,37; 0,0,36; 0,0,35; 0,0,35; 0,0,34; 0,0,34; 0,0,33; 0,0,33; 0,0,32; 0,0,31; 0,0,31; 0,0,30; 0,0,30; 0,0,29; 0,0,29; 0,0,28; 0,0,27; 0,0,27; 0,0,26; 0,0,26; 0,0,25; 0,0,24; 0,0,24; 0,0,23; 0,0,23; 0,0,22; 0,0,22; 0,0,21; 0,0,20; 0,0,20; 0,0,19; 0,0,19; 0,0,18; 0,0,17; 0,0,17; 0,0,16; 0,0,16; 0,0,15; 0,0,15; 0,0,14; 0,0,13; 0,0,13; 0,0,12; 0,0,12; 0,0,11; 0,0,10; 0,0,10; 0,0,9; 0,0,9; 0,0,8; 0,0,8; 0,0,7; 0,0,6; 0,0,6; 0,0,5; 0,0,5; 0,0,4; 0,0,3; 0,0,3; 0,0,2; 0,0,2; 0,0,1;0,0,1;0,0,0];
currentColormap = colormap(custom_colormatrix/255);

%Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(norm, size(currentColormap, 1));

%Limiting Colormap to data
clims = num2cell(get(gca, 'clim'));
[~, ~, ind] = histcounts(norm, linspace(clims{:}, size(currentColormap, 1)));

%Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%Repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

%Repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');

colorbar

% Figure 2 - Meta Stable Regions
set(gcf,'color','w');

%recalculate with finer resolution
[BU, FP, BO] = meshgrid(0:0.0125:0.55,0:0.0125:0.4,0:0.0125:0.7);
dBU = BU.*(uBU + aBUBU.*BU + aBUFP.*FP + aBUBO.*BO);
dFP = FP.*(uFP + aFPFP.*FP + aFPBU.*BU + aFPBO.*BO);
dBO = BO.*(uBO + aBOBO.*BO + aBOBU.*BU + aBOFP.*FP);
norm = sqrt((dBU.^2) + (dFP.^2) + (dBO.^2)); % redefining because of new resolution

%meta stable points
fp_meta_01 = [BU(norm(:,:,:) <= 0.01) FP(norm(:,:,:) <= 0.01) BO(norm(:,:,:) <= 0.01)];
fp_meta_0025 = [BU(norm(:,:,:) <= 0.0025) FP(norm(:,:,:) <= 0.0025) BO(norm(:,:,:) <= 0.0025)];

h = scatter3(fp_meta_01(:,1),fp_meta_01(:,2),fp_meta_01(:,3),'g','DisplayName','OD$/$hr $<$ 0.01');
alpha = 0.5;
set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
g = scatter3(fp_meta_0025(:,1),fp_meta_0025(:,2),fp_meta_0025(:,3),'g','DisplayName','OD$/$hr $<$ 0.0025','MarkerFaceColor',[0, 0.5, 0]);

%copied from the paper (approximate)
o1 = [0.08 0.2 0.125];
o2 = [0.1875 0.15 0.25];
o3 = [0.22 0.1 0.3245];
o4 = [0.25 0.13 0.324];
o5 = [0.25 0.14 0.325];
o6 = [0.326 0.15 0.23];
o7 = [0.325 0.24 0.13];
o = [o1; o2; o3; o4; o5; o6; o7];
plot3(o(:,1),o(:,2),o(:,3),':ks','DisplayName','Ophelia raw time points, approx','LineWidth',1)

[~, y] = ode23(@(t,y) GLV_plant(t,y,u,A), [0 50], o1');
plot3(y(:,1), y(:,2), y(:,3),'b','DisplayName','ode45(x0=ophelia t1, t=50 hrs)','LineWidth',2)

[~, y] = ode23(@(t,y) GLV_plant(t,y,u,A), [0 50], o3');
plot3(y(:,1), y(:,2), y(:,3),'y','DisplayName','ode45(x0=ophelia t3, t=50 hrs)','LineWidth',2)

[~, y] = ode23(@(t,y) GLV_plant(t,y,u,A), [0 50], o7');
plot3(y(:,1), y(:,2), y(:,3),'m','DisplayName','ode45(x0=ophelia t7, t=50 hrs)','LineWidth',2)

% this plots stable positions w (j eigvals<0), make sure latter code has
% been run once
% for k=1:length(stable)
%     plot3(stable(1,k),stable(2,k),stable(3,k),'b*')
%     hold on
% end
grid on
hold off

xlim([0 0.45])
ylim([0 0.35])
zlim([0 0.7])
xlabel('Bacteriodes uniformis','FontSize',15,'Interpreter','latex')
ylabel('Faecalibacterium prausnitzii','FontSize',15,'Interpreter','latex')
zlabel('Bacteroides ovatus','FontSize',15,'Interpreter','latex')
title('Reduced Community State Space','FontSize',30,'Interpreter','latex')
legend('FontSize',15,'Interpreter','latex','Location','best')

% Uncomment to make a gif
view(-10-zz/2, 15)
  
drawnow

frame = getframe(fig);
im = frame2im(frame);
[imind, map] = rgb2ind(im,256);    


%%%%%%% DONT FORGET TO RENAME THE GIF BEFORE RUNNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!
if zz == 1
    imwrite(imind, map,'OV_reduced_style2.gif','gif','LoopCount',Inf,'DelayTime',0.125);
else
    imwrite(imind, map,'OV_reduced_style2.gif','gif','WriteMode','append','DelayTime',0.125);
end
hold off

end

 %% Newton Evaluation, maybe should try a couple different i/c
% fps = ones(3,1);
% 
% for x0=1:length(for_nm)
%     [lambda, its, hist] = myNewton_2(G, H, {for_nm(1,x0);for_nm(2,x0);for_nm(3,x0)}, 10^-6, 100);
%     point = cell2mat(lambda);
%     if abs(N_f(point(1),point(2),point(3))) < 0.1 
%         fps = [fps point];
%     end
%     
%     fprintf("Percent Complete: %f \n",(x0/length(for_nm)))
% end
% 
%  disp("Done!")
 
 %Maybe its incorrect to use NM with the normal of the gradient? 
 %% Jacobian Eig Evaluation
 eigvals = zeros(3,length(fp_meta));
 stable = [;;];
 
%  figure('Position', [10 10 900 600])
%  set(gcf,'color','w');
 
 for x0=1:length(fp_meta)
    eigvals(:,x0)=eig(J(fp_meta(1,x0),fp_meta(2,x0),fp_meta(3,x0)));
%     plot3(eigvals(1,x0),eigvals(2,x0),eigvals(3,x0),'ko')
%     hold on
    
    if eigvals(1,x0)<0 && eigvals(2,x0)<0 && eigvals(3,x0)<0
        stable = [stable fp_meta(:,x0)];
%         plot3(eigvals(1,x0),eigvals(2,x0),eigvals(3,x0),'r*')
    end
%     fprintf("Percent Complete: %f \n",(x0/length(for_nm)))
 end

%  grid on
%  hold off
% 
% xlabel('$\lambda$ 1','FontSize',15,'Interpreter','latex')
% ylabel('$\lambda$ 2','FontSize',15,'Interpreter','latex')
% zlabel('$\lambda$ 3','FontSize',15,'Interpreter','latex')
% title('Eigenvalues of Jacobian where Norm $<$ 0.01','FontSize',30,'Interpreter','latex')
 
%% Eigen value plotting

for k=1:length(stable)
    plot3(stable(1,k),stable(2,k),stable(3,k),'b*')
    hold on
end
 hold off
 
 
    