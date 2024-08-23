clc,clear;
close all;
G = 1.0e10;
% Objectname = 'FH_0_nuu_0.3_gamma_0_pflag_3_c_4e-07';
% Objectname = "NewFH_0_nuu_0.262_gamma_0_pflag_3_c_3.7707e-07_factor_1";
% Objectname = "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_2.6982e-07_factor_1";
% Objectname = 'Reduced_gamma_0_pflag_3_c_4e-07';
% Objectname = 'Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1';
% Objectname = 'FluxTime_7.5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1';
Objectname = 'FluxTime_5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1';
V_dyn = 1;

% wide flag, trange
wideflag = 0;

% Xrange
Xrange = [-50, 50];
Xticks = [-30, 0, 30];

crange = [-13, 1];

if wideflag == 1
    Xrange = [-200, 200];
    Xticks = [-120, 0, 120];
end

% Set T ticks and range
Yticks = 0:800:2400;
Trange = [0, 4000];

% Yticks = 920:5:930;
% Trange = [920, 930];

% Initialize names
filename = strcat('../outputMats/', Objectname, '.mat');
videoname = strcat(Objectname, '_logVP_withDim_wide_closeTime_', num2str(wideflag), '.mp4');

% Initialize video
myVideo = VideoWriter(strcat('../mp4files/', videoname), 'MPEG-4');
myVideo.FrameRate = 10;
myVideo.Quality = 100;
open(myVideo);

% Load .mat file
load(filename);
tmax = Trange(2);
si0 = 4.0e6;
ind = find(tsaveplot > tmax);
psave = psave(:, 1:ind(1));
pcsave = pcsave(:, 1:ind(1));
sigrsave = sigrsave(:, 1:ind(1));
Vsave = Vsave(:,1:ind(1));
sigrnsave = 4 * (psave - 1/4 * sigrsave - 1/2 * pcsave);
tsaveplot = tsaveplot(1:ind(1));

% Some constants
fontsize = 20;


% Non-dimensionalize fault length
L_nu = G * L / (b(1) - a(1)) / si0;
% L_nu = 1;


% Non-dimensionalize time, time to diffuse by 1 nucleation length
t_ = L_nu * L_nu / 0.2;

figg = figure(1);
figg.Position = [1000, 597, 500, 700];


% Find the mask of 0.5 Mpa
pc_ = 0.5e6; % 0.5e6;
mask_pc = zeros(2, size(psave, 2));
for iiii =1:1:size(psave, 2)
    id = find(psave(:,iiii) > pc_);
    if isempty(id)
        mask_pc(1, iiii) = 0.;
        mask_pc(2, iiii) = 0.;
    else
        mask_pc(1, iiii) = x(id(1));
        mask_pc(2, iiii) = x(id(end));
    end
end

[XX,YY] = meshgrid(tsaveplot, x);
%% Pre-plot figure 1
subplot(2,1,1);
%% logV plot   

ind = [];

% Remove negative slip rates for log scale plot:
for i = 1:length(tsaveplot)
     V_neg = [];
     V_vector = Vsave(:,i);
     V_neg = V_vector(V_vector<0);
     if isempty(V_neg)
         ind(end+1) = i;
     end
end

V = Vsave(:,ind);
time = tsaveplot(ind);
V_log=log10(V);

% Color V
h=pcolor((YY)', (XX)', V_log');
shading interp;
hold on; 
set(h, 'EdgeColor', 'none')
caxis([-13, 1]);
colormap(gca,brewermap([],'YlOrRd'));
caxis(crange);
ylim(Trange);
xlim(Xrange);
c=colorbar;
ylabel('Time [s]', 'interpreter', 'latex'); xlabel('$X [\mathrm{m}]$', 'interpreter', 'latex'); 
title('Evolution of Log Slip Rate' , 'interpreter', 'latex')
set(gca, 'FontSize', 20);

set(c,'LineWidth',1);
ylabel(c,'$\log(V [\mathrm{m/s}])$','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
set(gca, 'TickLength', [.01 .01],...
'TickDir','in',...
'XMinorTick', 'on','YMinorTick', 'on','FontName',...
'Avenir','FontSize',fontsize) 
box on; 
set(gca,'layer','top')
set(gca,'LineWidth',1);
set(gca,'xtick',Xticks);
set(gca,'ytick',Yticks);

% P mask 0.5 MPa
plot( mask_pc(1,:), tsaveplot,'--k', 'linewidth', 1.5);
hold on; grid on;
plot(mask_pc(2,:), tsaveplot, '--k', 'linewidth', 1.5);

%% Loop to plot the time line and time evolution of 3 p's
% Initialize time handle
jjj = 1;
while jjj < size(pcsave, 2)
%% Draw a line on the Pm plot
    subplot(2,1,1);
    if jjj ~= 1
        delete(lastline);
    end
    lastline = yline(tsaveplot(jjj), 'LineWidth', 2.0);
%% Second plot -- Evolution of p+, p- and p_c  
    subplot(2,1,2);
    plot(x, sigrsave(:, jjj) / 1.e6, 'linewidth', 2.5);

    hold on; grid on;
    plot(x, pcsave(:, jjj) / 1.e6, '--', 'linewidth', 2.0)
    plot(x, sigrnsave(:, jjj) / 1.e6, 'linewidth', 1.5);
    xlim(Xrange);
    ylim([-2.5, 5]);
    xlabel('$X [\mathrm{m}]$', 'interpreter', 'latex');
    ylabel ('$\delta p$ [MPa]', 'interpreter', 'latex');
    hold on; grid on;
    title(strcat('Simulated Time $t = $', num2str(tsaveplot(jjj), '%.4f'), ' s'), 'interpreter', 'latex');
    legend('$\delta p^+$','$\delta p_c$', '$\delta p^-$', 'location', 'best', 'interpreter', 'latex');
    set(gca, 'FontSize', 20);
    xticks(Xticks);
    hold off;
    pause(0.001);
    indd = find(tsaveplot > tsaveplot(jjj) + 10);
    if isempty(indd)
        break;
    end
    jjj = min(indd(1), jjj + 50);
    
    %% Write the video
    set(gcf,'color','w');
    frame = getframe(gcf);
    writeVideo(myVideo, frame);  
end
close(myVideo);