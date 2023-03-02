clc,clear;
close all;
% Filename to load
prename1 = 'Original_gamma_0_pflag_3_c_0.0004';
filename1 = strcat(prename1, '.mat');
prename2 = 'Original_gamma_0.00017_pflag_3_c_0.0004';
filename2 = strcat(prename2, '.mat');
fontsize = 24;
figure(1);
% load(filename1);
% plot(tsaveplot, (pcsave(size(pcsave, 1)/2,:)+pcsave(size(pcsave, 1)/2 + 1,:))/(2) , 'linewidth', 2.0);
% hold on; grid on;

tsaveplot = 0:0.1:2500;
pc_center = zeros(1, size(tsaveplot, 2));
for i = 1:1:size(tsaveplot, 2)    
    % Prescribed pc at the center of the injection
    if tsaveplot(i) > 2144.2
        pc_center(i) = -191216.9;
    elseif tsaveplot(i) > 1400
        pc_center(i) = -4544.5*(tsaveplot(i))+9553100;
    elseif tsaveplot(i) > 1250
        pc_center(i) = 3190800;            
    else
        pc_center(i) = -2.073*(tsaveplot(i))^2+5144*(tsaveplot(i)); 
    end
end
plot(tsaveplot, (pc_center + 191216.9) / 1e6, 'linewidth', 2.0);
hold on; grid on;
% legend('Simulated Pressure', 'Targeted Pressure');
xlabel('Time [s]', 'interpreter', 'latex');
ylabel('Pressure [MPa]', 'interpreter', 'latex');
ylim([-1, 4]);
xlim([0, 2000]);
title('Prescribed Central Fault Fluid Pressure', 'interpreter', 'latex');
set(gca, 'fontsize', fontsize);


% load(filename2);
% figure(2);
% plot(tsaveplot, (pcsave(size(pcsave, 1)/2,:)+pcsave(size(pcsave, 1)/2 + 1,:))/(2) + 1.912e5 , 'linewidth', 2.0);
% hold on; grid on;
% 
% pc_center = zeros(1, size(tsaveplot, 2));
% for i = 1:1:size(tsaveplot, 2)    
%     % Prescribed pc at the center of the injection
%     if tsaveplot(i) > 2144.2
%         pc_center(i) = 0;
%     elseif tsaveplot(i) > 1400
%         pc_center(i) = -4544.5*(tsaveplot(i))+9.7443e6;
%     elseif tsaveplot(i) > 1250
%         pc_center(i) = 3.382e6;            
%     else
%         pc_center(i) = -2.073*(tsaveplot(i))^2+5144*(tsaveplot(i)) + 1.912e5; 
%     end
% end
% plot(tsaveplot, pc_center, '--', 'linewidth', 4.0);
% legend('Simulated Pressure', 'Targeted Pressure');
% ylim([0, 5e6]);
% xlabel('Time (s)');
% ylabel('Fluid Pressure (Pa)');
% title('Central Fault Fluid Pressure')