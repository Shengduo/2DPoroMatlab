clc,clear; close all; 
matfiles = ["Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
            "Elastic_Flag2_FluxTime_0.0001_NewFH_0_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
            "Elastic_Flag2_FluxTime_0.0001_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_6.1712e-09_factors_1_1_1"]; 

fig = figure();
set(fig, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
lws = [5.0, 3.0, 1.0]; 
legends = ["Poroelastic", "Elastic, $c$", "Elastic, $c_{mass}$"];
for i = 1:1:3
    filename = strcat('../outputMats/', matfiles(i), '.mat');
    load(filename, 'pcsave', 'psave', 'thetasave', 'x', 'tsaveplot', ...
        'Vsave', 'dphi0', 'gamma', 'Vr', 'L', 'dsave', 'tauS', 'sisave');
    fontsize = 24;
    semilogy(tsaveplot, (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2) , 'linewidth', lws(i));
    hold on; grid on;
end
legend(legends, 'location', 'best', 'interpreter', 'latex');
xlabel("Time [s]", "Interpreter","latex");
ylabel("V [m/s]", "Interpreter","latex");
xlim([0, 2000]);
set(gca, 'Fontsize', fontsize);
set(gcf, "color", "white");