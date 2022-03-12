clc,clear;
close all;
G = 1.0e10;
Objectname = 'Reduced_gamma_0_pflag_3_c_4e-07';
% load .mat file
filename = strcat('../outputMats', Objectname, '.mat');
videoname = strcat(Objectname, '.avi');

% Initialize video
myVideo = VideoWriter(videoname);
myVideo.FrameRate = 10;

% Load the .mat file
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

% Non-dimensionalize time
t_ = L/Vr;

% Non-dimensionalize fault length
% L_nu = G * L / (b(1) - a(1)) / si0;
L_nu = 1;
% 1.5* Every time
jjj = 1;

figg = figure(1);
figg.Position = [1000, 597, 560, 700];


% 1.5* Every time
% jjj = 1;
indd = find((Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2) > 7e-7);
jjj = indd(1);
% jjj = min(indd(1), jjj + 50);
while jjj < size(pcsave, 2)
    subplot(2,1,1);
    semilogy(tsaveplot, (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2) , 'linewidth', 3.0);
    hold on; grid on;
    xlabel('Time [s]'); ylabel('Slip Rate [m/s]'); title('Slip Rate at Fault Center')
    xlim([0, tmax]);
    ylim([1e-30, 1e0]);
    set(gca, 'FontSize', 20);
    xline(tsaveplot(jjj), 'LineWidth', 2.0);
    hold off;
    
    subplot(2,1,2);
    plot(x ./ L_nu, sigrsave(:, jjj)/si0, 'linewidth', 2.5);

    hold on; grid on;
    plot(x ./ L_nu, pcsave(:, jjj)/si0, '--', 'linewidth', 2.0)
    plot(x ./ L_nu, sigrnsave(:, jjj)/si0, 'linewidth', 1.5);
    xlim([-50, 50]);
    ylim([-0.5, 1]);
    %ylim([-0.0002, 0.0002])
    % xlabel('X / L_{nu}');
    xlabel('X [m]');
    ylabel ('\delta p / (\sigma_0 - p_0)');
    hold on; grid on;
    title(strcat('Simulated Time $t =  $', num2str(tsaveplot(jjj), '%.1f'), ' s'), 'intepreter', 'latex');
    legend('$\delta p^+$','$\delta p_c$', '$\delta p^-$', 'location', 'best', 'intepreter', 'latex');
    set(gca, 'FontSize', 20);
    xticks(-50:25:50);
    hold off;
    
    pause(0.01);
    indd = find(tsaveplot > tsaveplot(jjj) + 20);
    jjj = min(indd(1), jjj + 50);
   
    %% Write the video
    % frame = getframe(gcf);
    % writeVideo(myVideo, frame);    
end