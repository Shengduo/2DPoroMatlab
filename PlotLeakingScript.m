%% Load 
clc,clear;
close all;
% filenames = ["NewFH_0_nuu_0.35_gamma_0.00017_pflag_3_c_2.6982e-07_factor_1", ...         % 68
%             "NewFH_0_nuu_0.35_gamma_0.00017_pflag_3_c_2.6982e-08_factor_1", ...         % 71
%             "NewFH_0_nuu_0.262_gamma_0.00017_pflag_3_c_3.7707e-07_factor_1", ...         % 74
%             "NewFH_0_nuu_0.262_gamma_0.00017_pflag_3_c_3.7707e-08_factor_1"];            % 77

filenames = ["NewFH_0_nuu_0.35_gamma_1.7e-05_pflag_3_c_2.6982e-08_factor_1", ...
             "NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_3.7707e-08_factor_1", ...
             "NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_4e-08_factor_1_BC_0.43206_2.8461e-08"];
fontsize = 24;
lws = linspace(4.0, 1.0, length(filenames));

for i = 1:1:size(filenames, 2)
    filename = strcat('../outputMats/', filenames(i), '.mat');
    load(filename, 'pcsave', 'psave', 'x', 'tsaveplot', ...
         'sigrsave', 'sisave', 'epsi', 'kappac', 'rhof0', ...
         'INprofile', 'kappacx');
     
    nOfTimeSteps = size(tsaveplot, 2);

    sigrnsave = 4 * psave - 2 * pcsave - sigrsave;
    flux = kappac .* ((pcsave - sigrsave) + (pcsave - sigrnsave)) ./ epsi;
    cum_flux = zeros(1, nOfTimeSteps);
    leaking = zeros(1, nOfTimeSteps);
    flux_total = zeros(1, nOfTimeSteps);
    Injector_total = zeros(1, nOfTimeSteps);
    % Find horizontal flux:
    IID = find(INprofile > 0);
    DX = x(IID(1)) - x(IID(1) - 1);
             
    for t = 1 : 1 : nOfTimeSteps
        cum_flux(t) = trapz(x, flux(:, t)');
        leaking(t) = trapz(tsaveplot, cum_flux);
        flux_total(t) = kappacx * (psave(IID(1), t) - psave(IID(1) - 1, t)) / DX * 2 * epsi + ...
                    kappacx * (psave(IID(2), t) - psave(IID(1) + 1, t)) / DX * 2 * epsi + ...
                    trapz(x(IID), flux(IID, t)');
        Injector_total(t) = trapz(tsaveplot, flux_total);
    end

    % Change flux into Kg/m
    leaking = leaking * rhof0;
    injected_total = Injector_total * rhof0;
    plot(tsaveplot, leaking, 'linewidth', lws(i)); hold on; grid on;
    [val, id] = min(abs(tsaveplot - 2000));
    disp(strcat("Case name: ", filenames(i)));
    disp(strcat("leaked / total injected at 2000 s: ", num2str(leaking(id) / injected_total(id), '%2f')));
end

Trange = [0, 2000];
Lrange = [0, 3e-2];
xlim(Trange);
ylim(Lrange);
xlabel('Time [s]', 'interpreter', 'latex');
ylabel('Fluid leakage [Kg/m]', 'interpreter', 'latex');
title('Injected fluid mass $\gamma=1.7\times 10^{-5}$', 'interpreter', 'latex');    
% legend('$\nu_u = 0.35, c = 4 \cdot 10^{-7}$', ...
%        '$\nu_u = 0.35, c = 4 \cdot 10^{-8}$', ...
%        '$\nu_u = 0.262, c = 4 \cdot 10^{-7}$', ...
%        '$\nu_u = 0.262, c = 4 \cdot 10^{-8}$', ...
%        'location', 'best', 'interpreter', 'latex');

legend("Base 3", ...
       "Case 3", ...
       "Case 3'", ...
       'location', 'best', 'interpreter', 'latex');

set(gca, 'fontsize', fontsize);
savename = strcat(pwd, '/../dsvg_plots1/', 'CumMassControlGamma1.7e-5c-8.png');
disp(savename);
% saveas(figure(2),savename);
print(figure(1) ,savename, '-dpng', '-r300');