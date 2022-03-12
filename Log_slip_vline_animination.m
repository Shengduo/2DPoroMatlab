clc,clear;
close all;
G = 1.0e10;
Objectname = 'FH_0_nuu_0.262_gamma_0_pflag_3_c_4e-07_factor_1';

% Initialize names
filename = strcat(Objectname, '.mat');
videoname = strcat(Objectname, '_V.mp4');

% Initialize video
myVideo = VideoWriter(strcat('mp4files/', videoname), 'MPEG-4');
myVideo.FrameRate = 15;
myVideo.Quality = 100;
open(myVideo);

% Load .mat file
load(filename);
tmax = 2000;
si0 = 4.0e6;
ind = find(tsaveplot > tmax);
psave = psave(:, 1:ind(1));
pcsave = pcsave(:, 1:ind(1));
sigrsave = sigrsave(:, 1:ind(1));
Vsave = Vsave(:,1:ind(1));
sigrnsave = 4 * (psave - 1/4 * sigrsave - 1/2 * pcsave);
tsaveplot = tsaveplot(1:ind(1));

% Some constants
fontsize = 24;
% Xrange
Xrange = [-50, 50];
Xticks = 0:500:1500;
Trange = [0, 2000];
Yticks = [-30, 0, 30];
crange = [-13, 1];

% Non-dimensionalize time
t_ = L/Vr;

% Non-dimensionalize fault length
% L_nu = G * L / (b(1) - a(1)) / si0;
L_nu = 1;

figg = figure(1);
figg.Position = [1000, 597, 500, 700];


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
h=pcolor(XX,YY,V_log);
shading interp;
hold on; 
set(h, 'EdgeColor', 'none')
ylabel('x [m]','FontSize',fontsize)
caxis([-13, 1]);
colormap(gca,brewermap([],'YlOrRd'));
caxis(crange);
xlim(Trange);
ylim(Xrange);
c=colorbar;
xlabel('Time [s]'); ylabel('X [m]'); title('Evolution of Log Slip Rate V [m/s]')
set(gca, 'FontSize', 20);

set(c,'LineWidth',1);
ylabel(c,'Log Slip Rate [m/s]','FontName','Avenir','FontSize',fontsize);
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
    plot(x ./ L_nu, log10(Vsave(:, jjj)), 'linewidth', 2.0);

    hold on; grid on;
    xlim([-50, 50]);
    ylim([-26, 1]);
    xlabel('X [m]');
    ylabel ('Log Slip Rate V [m/s]');
    hold on; grid on;
    title(strcat('Simulated Time t =  ', num2str(tsaveplot(jjj), '%.1f'), ' s'));
    set(gca, 'FontSize', 20);
    xticks(-50:25:50);
    hold off;
    title(strcat('Simulated Time t =  ', num2str(tsaveplot(jjj), '%.1f'), ' s'));
 %% Update next time step
    pause(0.001);
    indd = find(tsaveplot > tsaveplot(jjj) + 20);
    if isempty(indd)
        break;
    end
    jjj = min(indd(1), jjj + 50);
    
    %% Write the video
    frame = getframe(gcf);
    writeVideo(myVideo, frame);  
end
close(myVideo);