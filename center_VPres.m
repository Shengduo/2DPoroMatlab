function center_VPres(prename, saveflag, BigT)

    filename = strcat('../outputMats/', prename, '.mat');
    load(filename, 'pcsave', 'psave', 'sigrsave', 'tsaveplot', 'thetasave', ...
        'Vsave', 'G', 'si0', 'L', 'Vr', 'dsave', 'sisave', 'tauS');
    fontsize = 34;
    sigrnsave = 4 * (psave - 1/4 * sigrsave - 1/2 * pcsave);
    % Xrange
    % Xrange = [-50, 50];
    
    Trange = [0, BigT];
    Vrange = [1e-25, 1];
    Prange = [-0.6, 1.2]; %[-2, 4];
    % Plot the center slip rate vs. Time
    fig1 = figure(1);
    set(fig1, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
    semilogy(tsaveplot, (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2) , 'linewidth', 2.0);
    xlim(Trange);
    ylim(Vrange);
    xlabel('Time [s]', 'interpreter', 'latex');
    ylabel('Slip Rate [m/s]', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../plots_2/', prename, '_centerV.png');
        disp(savename);
        print(figure(1) ,savename, '-dpng', '-r350');
    end
    
    fig2 = figure(2);
    set(fig2, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
    plot(tsaveplot, (pcsave(size(pcsave, 1)/2,:)+pcsave(size(pcsave, 1)/2 + 1,:))/(2e6), 'linewidth', 4.0);
    hold on; grid on;
    plot(tsaveplot, (sigrsave(size(sigrsave, 1)/2,:)+sigrsave(size(sigrsave, 1)/2 + 1,:))/(2e6), '--', 'linewidth', 3.5);
    plot(tsaveplot, (sigrnsave(size(sigrnsave, 1)/2,:)+sigrnsave(size(sigrnsave, 1)/2 + 1,:))/(2e6), 'linewidth', 3.0);
    plot(tsaveplot, (psave(size(psave, 1)/2,:)+pcsave(size(psave, 1)/2 + 1,:))/(2e6), 'linewidth', 2.5);
    
    xlim(Trange);
    ylim(Prange);
    legend('$\delta p_c$', '$\delta p^+$', '$\delta p^-$', '$\delta p_m$', 'location', 'best', 'interpreter', 'latex');
    xlabel('Time [s]', 'interpreter', 'latex');
    ylabel('Fluid Pressure [MPa]', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../plots_2/', prename, '_centerP.png');
        disp(savename);
        print(figure(2) ,savename, '-dpng', '-r350');
    end
    
    
    % Plot the center theta vs. Time
    fig3 = figure(3);
    set(fig3, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
    thetarange = [10^(-15), 10^15];
    semilogy(tsaveplot, (thetasave(size(thetasave, 1)/2,:)+thetasave(size(thetasave, 1)/2 + 1,:))/(2) , 'linewidth', 2.0);
    xlim(Trange);
    ylim(thetarange);
    xlabel('Time [s]', 'interpreter', 'latex');
    ylabel('$\theta$ [s]', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../plots_2/', prename, '_centerTheta.png');
        disp(savename);
        print(figure(3) ,savename, '-dpng', '-r350');
    end
    
    %% Plots vs. slip
    % Slip rate vs. slip
    drange = [0, 50000];
    fig4 = figure(4);
    set(fig4, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
    semilogy((dsave(size(dsave, 1)/2,:)+dsave(size(dsave, 1)/2 + 1,:))/(2e-6), (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2), 'linewidth', 3.0);
    hold on; grid on;
    xlim(drange);
    ylim(Vrange);
    xlabel('Slip [$\mu$m]', 'interpreter', 'latex');
    ylabel('Slip rate [m/s]', 'interpreter', 'latex');
    title('Slip rate vs. slip at the center', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../plots_2/', prename, '_centerVvsD.png');
        disp(savename);
        print(figure(4) ,savename, '-dpng', '-r350');
    end
    
    % Friction coefficient vs. slip
    frange = [0.4, 1];
    drange = [0, 15];
    fig5 = figure(5);
    set(fig5, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
    plot((dsave(size(dsave, 1)/2,:)+dsave(size(dsave, 1)/2 + 1,:))/(2e-3), (tauS(size(tauS, 1)/2,:)+tauS(size(tauS, 1)/2 + 1,:))./(sisave(size(sisave, 1)/2,:)+sisave(size(sisave, 1)/2 + 1,:)), 'linewidth', 3.0);
    hold on; grid on;
    xlim(drange);
    ylim(frange);
    xlabel('Slip [mm]', 'interpreter', 'latex');
    ylabel('Friction coefficient', 'interpreter', 'latex');
    % title('Friction coefficient vs. slip at the center', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../plots_2/', prename, '_centerfvsD.png');
        disp(savename);
        print(figure(5) ,savename, '-dpng', '-r350');
    end
    
    % Slip rate vs. slip
    fig6 = figure(6);
    set(fig6, 'Units', 'inches', 'Position', [0    10    7.7778    5.8333]);
    semilogy((dsave(size(dsave, 1)/2,:)+dsave(size(dsave, 1)/2 + 1,:))/(2e-6), (thetasave(size(thetasave, 1)/2,:)+thetasave(size(thetasave, 1)/2 + 1,:))/(2), 'linewidth', 3.0);
    hold on; grid on;
    xlim(drange);
    ylim(thetarange);
    xlabel('Slip [$\mu$m]', 'interpreter', 'latex');
    ylabel('$\theta$ [s]', 'interpreter', 'latex');
    title('$\theta$ vs. slip at the center', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../plots_2/', prename, '_centerThetavsD.png');
        disp(savename);
        print(figure(6) ,savename, '-dpng', '-r350');
    end
end