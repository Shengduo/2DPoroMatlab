%% Load 
clc,clear;
close all;
filenames = ["FH_0_nuu_0.35_gamma_1.7e-05_pflag_3_c_4e-07_factor_1", ...         % 68
            "FH_0_nuu_0.35_gamma_1.7e-05_pflag_3_c_4e-08_factor_1", ...         % 71
            "FH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_4e-07_factor_1", ...         % 74
            "FH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_4e-08_factor_1"];            % 77

fontsize = 24;

for i = 1:1:size(filenames, 2)
    filename = strcat('../outputMats/', filenames(i), '.mat');
    load(filename, 'pcsave', 'psave', 'x', 'tsaveplot', ...
         'sigrsave', 'sisave', 'epsi', 'kappac', 'rhof0');
    nOfTimeSteps = size(tsaveplot, 2);

    sigrnsave = 4 * psave - 2 * pcsave - sigrsave;
    flux = kappac .* ((pcsave - sigrsave) + (pcsave - sigrnsave)) ./ epsi;
    cum_flux = zeros(1, nOfTimeSteps);
    leaking = zeros(1, nOfTimeSteps);

    for t = 1 : 1 : nOfTimeSteps
        cum_flux(t) = trapz(x, flux(:, t)');
        leaking(t) = trapz(tsaveplot, cum_flux);
    end

    % Change flux into Kg/m
    leaking = leaking * rhof0;
    plot(tsaveplot, leaking, 'linewidth', 2.0); hold on; grid on;

end

Trange = [0, 2000];
Lrange = [0, 3e-2];
xlim(Trange);
ylim(Lrange);
xlabel('Time [s]', 'interpreter', 'latex');
ylabel('Fluid leakage [Kg/m]', 'interpreter', 'latex');
title('Injected fluid mass $\gamma=1.7\times 10^{-5}$', 'interpreter', 'latex');    
legend('$\nu_u = 0.35, c = 4 \cdot 10^{-7}$', ...
       '$\nu_u = 0.35, c = 4 \cdot 10^{-8}$', ...
       '$\nu_u = 0.262, c = 4 \cdot 10^{-7}$', ...
       '$\nu_u = 0.262, c = 4 \cdot 10^{-8}$', ...
       'location', 'best', 'interpreter', 'latex');
set(gca, 'fontsize', fontsize);
savename = strcat(pwd, '/../dsvg_plots1/', 'CumMassAll_gamma1.7e-5.png');
disp(savename);
% saveas(figure(2),savename);
print(figure(1) ,savename, '-dpng', '-r500');