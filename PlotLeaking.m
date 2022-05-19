function PlotLeaking(prename)
    % pwdd = '/Users/shengduoliu/Documents/Elias_Matlab/cand_simus/stacy_replicate/Poroelastic_properties_comparison';
    %saveflag = 1;
    %pcflag = 0;

    % Filename to load
    %prename = 'Original_gamma_0_pflag_3_c_4e-08';
    %gamma = 1.7e-4;
    %Vr = 1.0e-6;
    filename = strcat('../outputMats/', prename, '.mat');
    load(filename, 'pcsave', 'psave', 'x', 'tsaveplot', ...
         'sigrsave', 'sisave', 'epsi', 'kappac', 'rhof0');
    nOfTimeSteps = size(tsaveplot, 2);
    
    fontsize = 24;
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
    
    Trange = [0, 2000];
    Lrange = [0, 3e-2];
    
    plot(tsaveplot, leaking, 'linewidth', 2.0); hold on; grid on;
    xlim(Trange);
    
    ylim(Lrange);
    xlabel('Time [s]', 'interpreter', 'latex');
    ylabel('Fluid leakage [Kg/m]', 'interpreter', 'latex');
    title('Cumulative fluid mass into the bulk');    
    set(gca, 'fontsize', fontsize);
    
    savename = strcat(pwd, '/../dsvg_plots1/', prename, '_cumMass.png');
    disp(savename);
    % saveas(figure(2),savename);
    print(figure(1) ,savename, '-dpng', '-r500');
end
    