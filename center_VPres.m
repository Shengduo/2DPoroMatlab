function center_VPres(prename, saveflag)

    filename = strcat(prename, '.mat');
    load(filename, 'pcsave', 'psave', 'sigrsave', 'tsaveplot', ...
        'Vsave', 'G', 'si0', 'L', 'Vr');
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
    xlabel('Time [s]');
    ylabel('Slip Rate [m/s]');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/plots/', prename, '_centerV.png');
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
    legend('\delta p_c', '\delta p^+', '\delta p^-', 'location', 'best');
    xlabel('Time [s]');
    ylabel('Fluid Pressure [MPa]');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/plots/', prename, '_centerP.png');
        disp(savename);
        print(figure(2) ,savename, '-dpng', '-r350');
    end
end