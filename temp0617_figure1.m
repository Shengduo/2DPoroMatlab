clc,clear; close all; 
% Comparing different Vini
% matfiles = ["Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%             "Elastic_Flag0_FluxTime_0.0001_NewFH_0_nu_nuu_0.24_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_2020_1000", ...
%             "Elastic_Flag0_FluxTime_0.0001_NewFH_0_nu_nuu_0.24_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_2020_1000000", ...
%             "Elastic_Flag0_FluxTime_0.0001_NewFH_0_nu_nuu_0.24_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_2020_1000000000"]; 
% 
% legends = ["$10^{-22}$ [m/s]", "$10^{-19}$ [m/s]", "$10^{-16}$ [m/s]", "$10^{-13}$ [m/s]" ]; 
% lws = [5.0, 4.0, 3.0, 2.0]; 
% tlim = [0, 1000];
%% Comparing elastic with poroelastic, keeping c the same
matfiles = ["Elastic_Flag2_FluxTime_0.0001_NewFH_0_nuu_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
            "Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
            "Elastic_Flag2_FluxTime_0.0001_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"]; 

legends = ["Elastic, $\nu = 0.24$", "Poroelastic", "Elastic, $\nu=0.35$"]; 
lws = [5.0, 4.0, 3.0];
fluxs = [1.0e-4, 1.0e-4, 1.0e-4];
tlim = [850, 1000];
plotTs = [925.2, 925.2, 925.2];

%% Comparing elastic with poroelastic, keeping cmass the same
% matfiles = ["Elastic_Flag2_FluxTime_0.0001_NewFH_0_nu_nuu_0.24_0.24_gamma_0_pflag_3_c_6.1712e-09_factors_1_1_1", ...
%             "Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%             "Elastic_Flag2_FluxTime_0.0001_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_6.1712e-09_factors_1_1_1"]; 
% 
% legends = ["Elastic, $\nu = 0.24$", "Poroelastic", "Elastic, $\nu=0.35$"]; 
% lws = [5.0, 4.0, 3.0];
% fluxs = [1.0e-4, 1.0e-4, 1.0e-4];
% tlim = [850, 1000];
% plotTs = [925.2, 925.2, 925.2];

%% Comparing different injection rates
% matfiles = ["Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%             "FluxTime_7.5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%             "FluxTime_5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"]; 

% matfiles = ["Elastic_Flag2_FluxTime_0.0001_NewFH_0_nuu_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%             "Elastic_Flag2_FluxTime_7.5e-05_NewFH_0_nuu_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%             "Elastic_Flag2_FluxTime_5e-05_NewFH_0_nuu_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1"]; 

% matfiles = ["Elastic_Flag2_FluxTime_0.0001_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%             "Elastic_Flag2_FluxTime_7.5e-05_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%             "Elastic_Flag2_FluxTime_5e-05_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"]; 
% 
% legends = ["baseline flux", "0.75 baseline flux", "0.5 baseline flux"]; 
% fluxs = [1.0e-4, 0.75e-4, 0.5e-4];
% tlim = [0, 2000];
% lws = [5.0, 4.0, 3.0];

%% Intermittent vs. constant injection
% matfiles = ["FluxTime_5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%             "FluxTime_0.0001_1010_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"]; 
% 
% legends = ["Constant", "Intermittent"]; 
% lws = [5.0, 3.0];

fig1 = figure(1);
fig2 = figure(2); 
fig3 = figure(3); 
fig4 = figure(4);
fig5 = figure(5);
fig6 = figure(6);
fig7 = figure(7);
fig8 = figure(8); 
set(fig1, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
set(fig2, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
set(fig3, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
set(fig4, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
set(fig5, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
set(fig6, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
set(fig7, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
set(fig8, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);

% peakV = zeros(1, length(lws));
% VAtFricPeak = zeros(1, length(lws));    
for i = 1:1:length(lws)
    filename = strcat('../outputMats/', matfiles(i), '.mat');
    load(filename, 'pcsave', 'psave', 'thetasave', 'x', 'tsaveplot', ...
        'Vsave', 'dphi0', 'gamma', 'Vr', 'L', 'dsave', 'tauS', 'sisave');
    fontsize = 24;
    
    % Figure 1, V vs. time
    figure(1); 
    semilogy(tsaveplot, (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2) , 'linewidth', lws(i));
    hold on; grid on;
    
    % Figure 2, Fric vs. D
    figure(2); 
    plot((dsave(size(dsave, 1)/2,:)+dsave(size(dsave, 1)/2 + 1,:))/(2e-3), ...
         (tauS(size(tauS, 1)/2,:)+tauS(size(tauS, 1)/2 + 1,:))./(sisave(size(sisave, 1)/2,:)+sisave(size(sisave, 1)/2 + 1,:)), ...
         'linewidth', lws(i));
    hold on; grid on;

    % Figure 3, Fric vs. time
    figure(3); 
    plot(tsaveplot, ...
         (tauS(size(tauS, 1)/2,:)+tauS(size(tauS, 1)/2 + 1,:))./(sisave(size(sisave, 1)/2,:)+sisave(size(sisave, 1)/2 + 1,:)), ...
         'linewidth', lws(i));
    hold on; grid on;

    % Figure 4, theta vs. time
    figure(4)
    semilogy(tsaveplot, ...
             (thetasave(size(tauS, 1)/2,:)+thetasave(size(tauS, 1)/2 + 1,:))./(sisave(size(sisave, 1)/2,:)+sisave(size(sisave, 1)/2 + 1,:)), ...
             'linewidth', lws(i));
    hold on; grid on;

    % Figure 5, First peak slip rate and slip rate at Friction peak
    % thisFric = (tauS(size(tauS, 1)/2,:)+tauS(size(tauS, 1)/2 + 1,:))./(sisave(size(sisave, 1)/2,:)+sisave(size(sisave, 1)/2 + 1,:)); 
    % thisV = (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2);
    % [~, ii] = max(thisFric);
    % VAtFricPeak(i) = thisV(ii);

    % for ii = 1:1:length(thisFric)
    %     if thisV(ii) <= thisV(ii + 1)
    %         continue
    %     else
    %         peakV(i) = thisV(ii);
    %         break
    %     end
    % end
    figure(6)
    semilogy(tsaveplot * fluxs(i), (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2) , 'linewidth', lws(i));
    hold on; grid on;

    figure(7)
    plot(tsaveplot * fluxs(i), (psave(size(psave, 1)/2,:)+psave(size(psave, 1)/2 + 1,:))/(2e6) , 'linewidth', lws(i));
    hold on; grid on;

    figure(8); 
    first_idx = find(tsaveplot >= plotTs(i));
    this_time = tsaveplot(first_idx(1)); 
    plot(x, (psave(:, first_idx(1))) ./ 1.e6, 'linewidth', lws(i));
    hold on; grid on;
end

peakV = [0.0230071, 0.0265035, 0.014258, 0.00210182]; 
figure(1);
legend(legends, 'location', 'best', 'interpreter', 'latex');
xlabel("Time [s]", "Interpreter","latex");
ylabel("V [m/s]", "Interpreter","latex");
xlim(tlim);
ylim([10 ^ (-10), 1.]);
set(gca, 'Fontsize', fontsize);
set(gcf, "color", "white");


figure(2);
xlabel("slip [mm]", "Interpreter","latex");
ylabel("Friction", "Interpreter","latex");
yline(0.5375, "--", "linewidth", lws(1) + 0.5);
legend([legends, "Initial Friction"], 'location', 'best', 'interpreter', 'latex');
xlim([0, 1]);
set(gca, 'Fontsize', fontsize);
set(gcf, "color", "white");

figure(3);
xlabel("Time [s]", "Interpreter","latex");
ylabel("Friction", "Interpreter","latex");
xlim(tlim);
yline(0.5375, "--", "linewidth", lws(1) + 0.5);
% legend([legends, "Initial friction"], 'location', 'best', 'interpreter', 'latex');
set(gca, 'Fontsize', fontsize);
set(gcf, "color", "white");


figure(4);
% legend(legends, 'location', 'best', 'interpreter', 'latex');
xlabel("Time [s]", "Interpreter","latex");
ylabel("$\theta$ [s]", "Interpreter","latex");
xlim(tlim);
set(gca, 'Fontsize', fontsize);
set(gcf, "color", "white");

% figure(5);
% scatter([-22, -19, -16, -13], peakV, 150, "filled");
% ylim([1e-8, 1e0]);
% set(gca,'yscale','log');
% hold on; grid on;
% scatter([-22, -19, -16, -13], VAtFricPeak, 150, "filled");
% xlabel("$\log_{10}(V_{ini})$ [m/s]", "Interpreter", "latex");
% ylabel("$V$ [m/s]", "Interpreter", "latex")
% legend(["$V$ at first peak of $V$", "$V$ at first peak of friction"], "Interpreter", "latex", "location", "best");
% set(gca, 'xtick', [-22, -19, -16, -13]);
% set(gca, 'Fontsize', fontsize);
% set(gcf, "color", "white");

figure(6);
% legend(legends, 'location', 'best', 'interpreter', 'latex');
xlabel("Injected mass $[\mathrm{Kg / m}]$", "Interpreter","latex");
ylabel("V [m/s]", "Interpreter","latex");
xlim(tlim * fluxs(1));
set(gca, 'Fontsize', fontsize);
set(gcf, "color", "white");

figure(7);
legend(legends, 'location', 'best', 'interpreter', 'latex');
xlabel("Injected mass $[\mathrm{Kg / m}]$", "Interpreter","latex");
ylabel("$\delta p_m\ [\mathrm{MPa}]$", "Interpreter","latex");
xlim(tlim * fluxs(1));
set(gca, 'Fontsize', fontsize);
set(gcf, "color", "white");

figure(8);
xlabel("$x\ [\mathrm{m}]$", "Interpreter","latex");
ylabel("$\delta p_m\ [\mathrm{MPa}]$", "Interpreter","latex");
xlim([-200, 200]);
legend(legends, 'location', 'best', 'interpreter', 'latex');
% set(gca, 'xticks', [-30, 0, 30]);
set(gca, 'Fontsize', fontsize);
set(gcf, "color", "white");