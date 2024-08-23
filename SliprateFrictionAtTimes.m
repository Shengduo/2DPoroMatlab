function SliprateFrictionAtTimes(prename, saveflag, pcflag, subtraction_flag, plotT)
    % pwdd = '/Users/shengduoliu/Documents/Elias_Matlab/cand_simus/stacy_replicate/Poroelastic_properties_comparison';
    %saveflag = 1;
    %pcflag = 0;

    % Filename to load
    %prename = 'Original_gamma_0_pflag_3_c_4e-08';
    %gamma = 1.7e-4;
    %Vr = 1.0e-6;
    filename = strcat('../outputMats/', prename, '.mat');
    load(filename, 'G', 'pcsave', 'psave', 'thetasave', 'x', 'tsaveplot', ...
        'Vsave', 'dphi0', 'gamma', 'Vr', 'L', 'dsave', 'tauS', 'sisave', ...
        'f0', 'si0', 'a', 'b', 'nu', 'x');
    fontsize = 24;
    
    % pcsave -- center pressure p_c history
    % psave -- center effective (mean) p history
    % dsave -- displacement delta history
    % tauS -- shear traction history
    % sisave -- effective normal stress history, si = si_0 + si_yy - p_eff

    % A few parameters
    % f0 = 0.5375;
    % si0 = 4e6;
    tau0 = f0*si0; 

    % Nucleation length with time
    L_nusave = G * L / (b(1) - a(1)) ./ sisave ./ (1. - nu);
    % Xrange
    % Xrange = [-5, 5];
    % yticks = [-3, 0, 3];

    Xrange = [-50, 50];
    yticks = [-30, 0, 30];

    if subtraction_flag == 0
        pcsave = pcsave + 1.912e5;
        psave = psave + 1.912e5;
    end
    if pcflag == 1
        psave = pcsave;
    end
    % timeend = find(tsaveplot > 1.05*Trange(2));
    % if size(timeend, 2) == 0
    %     timeend = size(pcsave, 2) - 1;
    % else
    %     timeend = timeend(1);
    % end
    timeend = size(pcsave, 2) - 1;
    pcsave = pcsave(:, 1:timeend);
    psave = psave(:, 1:timeend);
    tsaveplot = tsaveplot(:, 1:timeend);
    thetasave = thetasave(:, 1:timeend);
    sisave = sisave(:, 1:timeend);
    tauS = tauS(:, 1:timeend);
    % Find the mask of 0.5 Mpa
    pc_ = 0.5e6;
    mask_pc = zeros(2, size(psave, 2));
    for iiii =1:1:size(psave, 2)
        id = find(psave(:,iiii) > pc_);
        if isempty(id)
            mask_pc(1, iiii) = 0.;
            mask_pc(2, iiii) = 0.;
        else
            mask_pc(1, iiii) = x(id(1));
            mask_pc(2, iiii) = x(id(end));
        end
    end
    %--------------------------------------------------------------------------
    fig1=figure(1);
    % set(fig1,'Position',[700 700 450 500]);      
    set(fig1, 'Units', 'inches', 'Position', [0    10    6.2500    4 * 2]);

    %% plot the two lines
    % lws = linspace(3.0, 1.0, length(plot_times));
    subplot(2, 1, 1);
    first_idx = find(tsaveplot >= plotT);
    this_time = tsaveplot(first_idx(1)); 
    plot(log10(Vsave(:, first_idx(1))), x, 'linewidth', 1.5, 'color', "#0072BD");
    hold on; grid on;
    ylim(Xrange);
    xlim([-23, 1]);
    xlabel("$\log_{10}(V)\ [\mathrm{m/s}]$", 'interpreter', 'latex');
    ylabel("$x\ [\mathrm{m}]$", 'interpreter', 'latex');
    title("t = " + num2str(this_time) + " [s]",  'interpreter', 'latex')
    set(gca,'ytick',yticks);
    set(gca, 'fontsize', fontsize); 
    set(gcf, 'color', 'w');

    subplot(2, 1, 2)
    plot(tauS(:, first_idx(1))./ sisave(:, first_idx(1)), x, 'linewidth', 1.5, ...
         'color', "#D95319");
    hold on; grid on;
    xlabel("Friction", 'interpreter', 'latex');
    ylabel("$x\ [\mathrm{m}]$", 'interpreter', 'latex');
    ylim(Xrange);
    set(gca,'ytick',yticks);
    set(gca, 'fontsize', fontsize);
    set(gcf, 'color', 'w');
    %% Save the figures
    if saveflag == 1
        savename = strcat(pwd, '/../dsvg_plots2/', prename, '_VandFricRightBefore.png');
        disp(savename);
        % saveas(figure(1),savename);
        print(figure(1) ,savename, '-dpng', '-r500');
    end
end