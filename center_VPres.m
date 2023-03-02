function center_VPres(prename, saveflag)

    filename = strcat('../outputMats/', prename, '.mat');
    load(filename, 'pcsave', 'psave', 'sigrsave', 'tsaveplot', 'thetasave', ...
        'Vsave', 'G', 'si0', 'L', 'Vr', 'dsave', 'sisave', 'tauS');
    fontsize = 24;
    sigrnsave = 4 * (psave - 1/4 * sigrsave - 1/2 * pcsave);
    % Xrange
    % Xrange = [-50, 50];
    
    Trange = [0, 2000];
    Vrange = [1e-25, 1];
    Prange = [-2, 4];
    % Plot the center slip rate vs. Time
    figure(1);
    semilogy(tsaveplot, (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2) , 'linewidth', 2.0);
    xlim(Trange);
    ylim(Vrange);
    xlabel('Time [s]', 'interpreter', 'latex');
    ylabel('Slip Rate [m/s]', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../plots/', prename, '_centerV.png');
        disp(savename);
        print(figure(1) ,savename, '-dpng', '-r350');
    end
    
    figure(2);
    plot(tsaveplot, (pcsave(size(pcsave, 1)/2,:)+pcsave(size(pcsave, 1)/2 + 1,:))/(2e6), 'linewidth', 3.0);
    hold on; grid on;
    plot(tsaveplot, (sigrsave(size(sigrsave, 1)/2,:)+sigrsave(size(sigrsave, 1)/2 + 1,:))/(2e6), '--', 'linewidth', 3.0);
    plot(tsaveplot, (sigrnsave(size(sigrnsave, 1)/2,:)+sigrnsave(size(sigrnsave, 1)/2 + 1,:))/(2e6), 'linewidth', 2.0);
    xlim(Trange);
    ylim(Prange);
    legend('$\delta p_c$', '$\delta p^+$', '$\delta p^-$', 'location', 'best', 'interpreter', 'latex');
    xlabel('Time [s]', 'interpreter', 'latex');
    ylabel('Fluid Pressure [MPa]', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../plots/', prename, '_centerP.png');
        disp(savename);
        print(figure(2) ,savename, '-dpng', '-r350');
    end
    
    
    % Plot the center theta vs. Time
    figure(3);
    thetarange = [10^(-15), 10^15];
    semilogy(tsaveplot, (thetasave(size(thetasave, 1)/2,:)+thetasave(size(thetasave, 1)/2 + 1,:))/(2) , 'linewidth', 2.0);
    xlim(Trange);
    ylim(thetarange);
    xlabel('Time [s]', 'interpreter', 'latex');
    ylabel('$\theta$ [s]', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../plots/', prename, '_centerTheta.png');
        disp(savename);
        print(figure(3) ,savename, '-dpng', '-r350');
    end
    
    %% Plots vs. slip
    % Slip rate vs. slip
    drange = [0, 5000];
    figure(4);
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
        savename = strcat(pwd, '/../plots/', prename, '_centerVvsD.png');
        disp(savename);
        print(figure(4) ,savename, '-dpng', '-r350');
    end
    
    % Friction coefficient vs. slip
    frange = [0, 1];
    % drange = [0, 100];
    figure(5);
    plot((dsave(size(dsave, 1)/2,:)+dsave(size(dsave, 1)/2 + 1,:))/(2e-6), (tauS(size(tauS, 1)/2,:)+tauS(size(tauS, 1)/2 + 1,:))./(sisave(size(sisave, 1)/2,:)+sisave(size(sisave, 1)/2 + 1,:)), 'linewidth', 3.0);
    hold on; grid on;
    xlim(drange);
    ylim(frange);
    xlabel('Slip [$\mu$m]', 'interpreter', 'latex');
    ylabel('Friction coefficient', 'interpreter', 'latex');
    title('Friction coefficient vs. slip at the center', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../plots/', prename, '_centerfvsD.png');
        disp(savename);
        print(figure(5) ,savename, '-dpng', '-r350');
    end
    
    % Slip rate vs. slip
    figure(6);
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
        savename = strcat(pwd, '/../plots/', prename, '_centerThetavsD.png');
        disp(savename);
        print(figure(6) ,savename, '-dpng', '-r350');
    end
end