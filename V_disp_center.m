clc,clear;
close all;

% Filename to load

prename1 = 'Original_gamma_0_pflag_3_c_4e-08';
filename1 = strcat(prename1, '.mat');
prename2 = 'FH_nuu_0.35_gamma_0_pflag_3_c_4e-08';
filename2 = strcat(prename2, '.mat');
prename3 = 'Original_gamma_0.00017_pflag_3_c_4e-08';
filename3 = strcat(prename1, '.mat');
prename4 = 'Reduced_gamma_0.00017_pflag_3_c_4e-08';
filename4 = strcat(prename2, '.mat');

figure(1);
load(filename1);
semilogy(tsaveplot, (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2) , 'linewidth', 6.0);
hold on; grid on;

figure(2);
plot(tsaveplot, (dsave(size(dsave, 1)/2,:)+dsave(size(dsave, 1)/2 + 1,:))/(2) , 'linewidth', 2.0);
hold on; grid on;

load(filename2);
figure(1);
semilogy(tsaveplot, (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2) , '--', 'linewidth', 3.0);
hold on; grid on;

figure(2);
plot(tsaveplot, (dsave(size(dsave, 1)/2,:)+dsave(size(dsave, 1)/2 + 1,:))/(2) , 'linewidth', 2.0);
hold on; grid on;

figure(1);
xlabel('Time [s]'); ylabel('Slip Rate [m/s]');
title('Slip Rate at X = 0 [m]');
xlim([0, 1400]);
%ylim([1e-8, 1e-1]);
%yticks = [1e-7, 1e-5, 1e-3];
legend('Original, No Dilatancy', 'Flash heating, No Dilatancy', 'location', 'best');
title('Slip Rate at the Center')
%legend('Large Poroelastic Effect', 'Small Poroelastic Effect', 'location', 'best');
set(gca, 'FontSize', 25, 'ytick', yticks);

figure(2);
xlabel('Time [s]'); ylabel('Slip [m]');
title('Slip at X = 0 [m]');
xlim([0, 1400]);
%ylim([1e-8, 1e-1]);
%yticks = [1e-7, 1e-5, 1e-3];
legend('Original, No Dilatancy', 'Flash heating, No Dilatancy', 'location', 'best');
title('Slip at the Center')
%legend('Large Poroelastic Effect', 'Small Poroelastic Effect', 'location', 'best');
set(gca, 'FontSize', 25, 'ytick', yticks);

% saveas(figure(1),[pwd '/plots/' prename1 '_full_Vcenter.png']);

% figure(2);
% load(filename1);
% plot(tsaveplot, (dsave(size(dsave, 1)/2,:)+dsave(size(dsave, 1)/2 + 1,:))/(2) , 'linewidth', 2.0);
% hold on; grid on;
% 
% load(filename2);
% plot(tsaveplot, (dsave(size(dsave, 1)/2,:)+dsave(size(dsave, 1)/2 + 1,:))/(2) , 'linewidth', 2.0);
% hold on; grid on;
% xlabel('Time (s)'); ylabel('Slip (m)');
% title('Nucleation time');
% %legend('Large Poroelastic Effect', 'Small Poroelastic Effect', 'location', 'best');
% set(gca, 'FontSize', 15);

