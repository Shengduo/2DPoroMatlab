function center_VPres_nodim(prename, saveflag)

    filename = strcat('../outputMats/', prename, '.mat');
    load(filename, 'pcsave', 'psave', 'sigrsave', 'tsaveplot', ...
        'Vsave', 'G', 'si0', 'L', 'Vr', 'a', 'b');
    fontsize = 24;
    sigrnsave = 4 * (psave - 1/4 * sigrsave - 1/2 * pcsave);
    % Xrange
    % Xrange = [-50, 50];
    % Non-dimensionalize fault length
    L_nu = G * L / (b(1) - a(1)) / si0;
    V_dyn = 1;
    % L_nu = 1;


    % Non-dimensionalize time, time to diffuse by 1 nucleation length
    t_ = L_nu * L_nu / 0.2;
    
    Trange = [0, 5];
    Vrange = [1e-25, 1];
    Prange = [-0.5, 2.0];
    % Plot the center slip rate vs. Time
    fig1 = figure(1);
    set(fig1,'Position',[700 700 400 300]);  
    semilogy(tsaveplot ./ t_, (Vsave(size(Vsave, 1)/2,:)+Vsave(size(Vsave, 1)/2 + 1,:))/(2) , 'linewidth', 2.0);
    xlim(Trange);
    ylim(Vrange);
    xlabel('$t / t_{nu}$', 'interpreter', 'latex');
    ylabel('$\log(V/V_{dyn})$', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../dsvg_plots_nodim1/', prename, '_centerV.png');
        disp(savename);
        print(figure(1) ,savename, '-dpng', '-r350');
    end
    
    fig2 = figure(2);
    set(fig2,'Position',[700 700 400 300]);  
    plot(tsaveplot ./ t_, (pcsave(size(pcsave, 1)/2,:)+pcsave(size(pcsave, 1)/2 + 1,:)) ./ si0, 'linewidth', 3.0);
    hold on; grid on;
    plot(tsaveplot ./ t_, (sigrsave(size(sigrsave, 1)/2,:)+sigrsave(size(sigrsave, 1)/2 + 1,:)) ./ si0, '--', 'linewidth', 3.0);
    plot(tsaveplot ./ t_, (sigrnsave(size(sigrnsave, 1)/2,:)+sigrnsave(size(sigrnsave, 1)/2 + 1,:)) ./ si0, 'linewidth', 2.0);
    xlim(Trange);
    ylim(Prange);
    legend('$\delta p_c$', '$\delta p^+$', '$\delta p^-$', 'location', 'best', 'interpreter', 'latex');
    xlabel('$t / t_{nu}$', 'interpreter', 'latex');
    ylabel('$\delta p / (\sigma_0 - p_0)$', 'interpreter', 'latex');
    set(gca, 'Fontsize', fontsize);
    
    % Save the figure
    if saveflag == 1
        savename = strcat(pwd, '/../dsvg_plots_nodim1/', prename, '_centerP.png');
        disp(savename);
        print(figure(2) ,savename, '-dpng', '-r350');
    end
end