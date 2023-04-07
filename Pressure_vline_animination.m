clc,clear;
close all;
G = 1.0e10;
% Objectname = 'NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_3.7707e-07_factor_1';
Objectname = 'NewSiRegFH_0_nuu_0.35_gamma_0_pflag_3_c_4e-08_factor_1_BC_0.85_4e-08';

% Initialize names
filename = strcat('../outputMats/', Objectname, '.mat');
videoname = strcat(Objectname, '_P.mp4');

% Initialize video
myVideo = VideoWriter(strcat('../mp4files/', videoname), 'Motion JPEG AVI');
myVideo.FrameRate = 15;
myVideo.Quality = 100;
open(myVideo);

% Load .mat file
load(filename);
tmax = 2000;
si0 = 4.0e6;
ind = find(tsaveplot > tmax);
if isempty(ind)
    ind = length(tsaveplot);
end
psave = psave(:, 1:ind(1));
pcsave = pcsave(:, 1:ind(1));
sigrsave = sigrsave(:, 1:ind(1));
Vsave = Vsave(:,1:ind(1));
sigrnsave = 4 * (psave - 1/4 * sigrsave - 1/2 * pcsave);
tsaveplot = tsaveplot(1:ind(1));

% Some constants
fontsize = 24;
% Xrange
% Xrange = [-50, 50];
Xrange = [-250, 250];
Xticks = 0:500:1500;
Trange = [0, 2000];
Yticks = [-30, 0, 30];
crange = [-2, 4];

% Non-dimensionalize time
t_ = L/Vr;

% Non-dimensionalize fault length
% L_nu = G * L / (b(1) - a(1)) / si0;
L_nu = 1;

figg = figure(1);
% figg.Position = [1000, 597, 500, 700]; % Mac
figg.Position = [1000, 597, 800, 1120]; % Linux


% Find the mask of 0.5 Mpa
pc_ = 0.5e6;
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
%% P plot    
% Color P
h=pcolor(XX,YY,psave/1e6);
shading interp;
hold on; 
set(h, 'EdgeColor', 'none')
ylabel('x [m]','FontSize',fontsize)
caxis([-4,4]);
Map = bluewhitered(400);
colormap(gca,Map([1:2:199 200:end],:));
caxis(crange);
xlim(Trange);
ylim(Xrange);
c=colorbar;
xlabel('Time [s]'); ylabel('X [m]'); title('Evolution of Pressure $\delta p_m$')
set(gca, 'FontSize', 20);

set(c,'LineWidth',1);
caxis(crange);
ylabel(c,'Pressure [MPa]','FontName','Avenir','FontSize',fontsize, 'Interpreter', 'latex');
set(gca, 'TickLength', [.01 .01],...
'TickDir','in',...
'XMinorTick', 'on','YMinorTick', 'on','FontName',...
'Avenir','FontSize',fontsize) 
box on; 
set(c, 'ylim', crange);
set(gca,'layer','top')
set(gca,'LineWidth',1);
set(gca,'xtick',Xticks);
set(gca,'ytick',Yticks);

% P mask 0.5 MPa
plot(tsaveplot, mask_pc(1,:), '--k', 'linewidth', 1.5);
hold on; grid on;
plot(tsaveplot, mask_pc(2,:), '--k', 'linewidth', 1.5);

%% Loop to plot the time line and time evolution of 3 p's
% Initialize time handle
jjj = 1;
while jjj < size(pcsave, 2)
%% Draw a line on the Pm plot
    subplot(2,1,1);
    if jjj ~= 1
        delete(lastline);
    end
    lastline = xline(tsaveplot(jjj), 'LineWidth', 2.0);
%% Second plot -- Evolution of p+, p- and p_c  
    subplot(2,1,2);
    plot(x ./ L_nu, sigrsave(:, jjj)/si0, 'linewidth', 2.5);

    hold on; grid on;
    plot(x ./ L_nu, pcsave(:, jjj)/si0, '--', 'linewidth', 2.0)
    plot(x ./ L_nu, sigrnsave(:, jjj)/si0, 'linewidth', 1.5);
    plot(x ./ L_nu, (sigrnsave(:, jjj) / 4 + pcsave(:, jjj) / 2 + sigrsave(:, jjj) / 4) / si0, 'k', 'linewidth', 3.0)
    xlim(Xrange);
    ylim([-0.5, 1]);
    xlabel('X [m]');
    ylabel ('$\delta p / (\sigma_0 - p_0)$', 'interpreter', 'latex');
    hold on; grid on;
    title(strcat('Simulated Time t =  ', num2str(tsaveplot(jjj), '%.1f'), ' s'));
    legend('$\delta p^+$','$\delta p_c$', '$\delta p^-$', '$\delta p_m$', 'location', 'best', 'interpreter', 'latex');
    set(gca, 'FontSize', 20);
    xticks(-50:25:50);
    hold off;
    set(gcf,'color','w');
    pause(0.001);
    indd = find(tsaveplot > tsaveplot(jjj) + 20);
    if isempty(indd)
        break;
    end
    jjj = min(indd(1), jjj + 50);

    %% Write the video
    frame = getframe(gcf);
    % writeVideo(myVideo, frame);  
end
close(myVideo);