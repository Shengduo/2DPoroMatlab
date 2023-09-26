function stacy_LnuVsPRegion(prename, saveflag, pcflag, subtraction_flag)
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

    xticks = 0:500:1500;
    Trange = [0, 2000];
    crange = [-2, 4];
    if subtraction_flag == 0
        pcsave = pcsave + 1.912e5;
        psave = psave + 1.912e5;
    end
    if pcflag == 1
        psave = pcsave;
    end
    timeend = find(tsaveplot > 1.05*Trange(2));
    if size(timeend, 2) == 0
        timeend = size(pcsave, 2) - 1;
    else
        timeend = timeend(1);
    end
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

    plot_times = [750, 1000, 1250, 1500];
    %--------------------------------------------------------------------------
    % % First the color map
    [XX,YY] = meshgrid(tsaveplot, x);


    ha = tight_subplot(2,1,0.040,[.15 .075],[.2 .16]);
    fig1=figure(1);
    % set(fig1,'Position',[700 700 450 500]);      
    set(fig1, 'Units', 'inches', 'Position', [0    10    6.2500    6.9444]);
    %% P plot
    axes(ha(1));
    %[X, Y] = meshgrid(time,x);
    h=pcolor(XX,YY,psave/1e6);
    shading interp;
    hold on; 
    set(h, 'EdgeColor', 'none')
    ylabel('x [m]','FontSize',fontsize, 'interpreter', 'latex')
    %colormap(ha(1),brewermap([],'REDBLUE'));
    % imagesc(peaks(250));
    caxis([-4,4]);
    c = colorbar;
    Map = bluewhitered(400);
    colormap(ha(1),Map([1:2:199 200:end],:));
    caxis(crange);
    xlim(Trange);
    ylim(Xrange);
    c=colorbar;

    %--------------------------------------------------------------------------
    P = polyfit(mask_pc(1,:), tsaveplot, 2);
    fittime = polyval(P, mask_pc(1,:));
    %plot(fittime, mask_pc(1,:), '--k', 'linewidth', 1.5);
    plot(tsaveplot, mask_pc(1,:), '--k', 'linewidth', 1.5);
    hold on; grid on;
    P = polyfit(mask_pc(2,:), tsaveplot, 2);
    fittime = polyval(P, mask_pc(2,:));
    % plot(fittime, mask_pc(2,:), '--k', 'linewidth', 1.5);
    plot(tsaveplot, mask_pc(2,:), '--k', 'linewidth', 1.5);
    %--------------------------------------------------------------------------

    %plot(time(1:10:end),r_press_0_5(1:10:end),'--k','Linewidth',1.5);
    %plot(time(1:10:end),-r_press_0_5(1:10:end),'--k','Linewidth',1.5);
    set(c,'LineWidth',1);
    caxis(crange);
    ylabel(c,'Pressure [MPa]','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
    set(gca, 'TickLength', [.01 .01],...
    'TickDir','in',...
    'XMinorTick', 'on','YMinorTick', 'on','FontName',...
    'Avenir','FontSize',fontsize) 
    box on; 
    set(c, 'ylim', crange);
    set(gca,'layer','top')
    set(gca,'LineWidth',1);
    set(gca,'xtick',xticks);
    set(gca,'ytick',yticks);
    % xlabel('Time [s]','FontSize',fontsize);
    set(ha(1),'XTickLabel','');
    % set(gca,'Color','none');
    box on; 

    %% log V plot
    ind = [];

    % Remove negative slip rates for log scale plot:

    for i = 1:length(tsaveplot)
         V_neg = [];
         V_vector = Vsave(:,i);
         V_neg = V_vector(V_vector<0);
         if isempty(V_neg)
             ind(end+1) = i;
         end
    end

    V = Vsave(:,ind);
    time = tsaveplot(ind);
    % r_press_0_5 = r_press_0_5(ind);

    axes(ha(2));
    V_log=log10(V);
    [X, Y] = meshgrid(time,x);
    h=pcolor(X,Y,V_log);
    shading interp;
    hold on; 
    set(h, 'EdgeColor', 'none')
    ylabel('x [m]','FontSize',fontsize, 'interpreter', 'latex');
    colormap(ha(2),brewermap([],'YlOrRd'));

    % Mask for 0.5 MPa
    %--------------------------------------------------------------------------
    P = polyfit(mask_pc(1,:), tsaveplot, 2);
    fittime = polyval(P, mask_pc(1,:));
    %plot(fittime, mask_pc(1,:), '--k', 'linewidth', 1.5);
    plot(tsaveplot, mask_pc(1,:), '--k', 'linewidth', 1.5);
    hold on; grid on;
    P = polyfit(mask_pc(2,:), tsaveplot, 2);
    fittime = polyval(P, mask_pc(2,:));
    % plot(fittime, mask_pc(2,:), '--k', 'linewidth', 1.5);
    plot(tsaveplot, mask_pc(2,:), '--k', 'linewidth', 1.5);
    %--------------------------------------------------------------------------

    ylim(Xrange);
    xlim(Trange);
    c=colorbar;
    set(c,'LineWidth',1);
    caxis([-13 1]);
    ylabel(c,'Log Slip Rate [m/s]','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
    set(gca, 'TickLength', [.01 .01],...
    'TickDir','in',...
    'XMinorTick', 'on','YMinorTick', 'on','FontName',...
    'Avenir','FontSize',fontsize) 
    box on; 
    set(gca,'layer','top')
    set(gca,'xtick',xticks);
    set(gca,'ytick',yticks);
    xlabel('Time [s]','FontSize',fontsize, 'interpreter', 'latex');
    set(gca,'LineWidth',1);
    xtickangle(0);
    box on;  
    
    fig2=figure(2);
    % set(fig1,'Position',[700 700 450 500]);      
    set(fig2, 'Units', 'inches', 'Position', [0    10    6.2500    4 * length(plot_times)]);

    %% plot the two lines
    % lws = linspace(3.0, 1.0, length(plot_times));
    for idx = 1:1:length(plot_times)
        subplot(length(plot_times), 1, idx);
        first_idx = find(tsaveplot > plot_times(idx));
        this_time = tsaveplot(first_idx(1));
        distance1 = 2 * (x - mask_pc(1, first_idx(1))); 
        distance1(distance1 < 0.) = 0.;
        distance2 = 2 * (mask_pc(2, first_idx(1)) - x);
        distance2(distance2 < 0.) = 0.;
        distance = min(distance1, distance2); 
        plot(distance, x, 'linewidth', 1.5, 'color', "#0072BD");
        hold on; grid on;
        plot(L_nusave(:, first_idx(1)), x, 'linewidth', 1.5, ...
             'color', "#D95319");
        if idx == length(plot_times)
            xlabel("Length [m]",  'interpreter', 'latex');
            legend("$L_{pr}$", "$L_{nu}$", "location", "best", "interpreter", "latex")
        end
        title("t = " + num2str(plot_times(idx)) + " [s]",  'interpreter', 'latex')
        ylabel("$x$ [m]", 'interpreter', 'latex');
        set(gca, "fontsize", fontsize);
        xlim([0, 50]);
        ylim(Xrange);
        % set(gca,'LineWidth',1);
        set(gca,'ytick',yticks);
        set(gcf, 'color', 'w');
    end

    fig3=figure(3);
    % set(fig1,'Position',[700 700 450 500]);      
    set(fig3, 'Units', 'inches', 'Position', [0    10    6.2500    4 * length(plot_times)]);
    
    %% plot the two lines
    % lws = linspace(3.0, 1.0, length(plot_times));
    for idx = 1:1:length(plot_times)
        subplot(length(plot_times), 1, idx);
        first_idx = find(tsaveplot > plot_times(idx));
        this_time = tsaveplot(first_idx(1));
        distance1 = 2 * (x - mask_pc(1, first_idx(1))); 
        distance1(distance1 < 0.) = 0.;
        distance2 = 2 * (mask_pc(2, first_idx(1)) - x);
        distance2(distance2 < 0.) = 0.;
        distance = min(distance1, distance2); 
%         plot(distance, x, 'linewidth', 1.5, 'color', "#0072BD");
        hold on; grid on;
        AwayFromSS = 1. - Vsave(:, first_idx(1)) .* thetasave(:, first_idx(1)) ./ L; 
%         plot(L_nusave(:, first_idx(1)), x, 'linewidth', 1.5, ...
%              'color', "#D95319");

%         if idx == length(plot_times)
%             xlabel("Length [m]",  'interpreter', 'latex');
%             legend("$L_{pr}$", "$L_{nu}$", "location", "best", "interpreter", "latex")
%         end
        plot(AwayFromSS, x, 'linewidth', 1.5);
        hold on; grid on; 
        plot(log10(Vsave(:, first_idx(1))), x, 'linewidth', 1.5);
        if idx == length(plot_times)
            legend("$1 - V\theta/D_{RS}$", "$\log V$", 'interpreter', 'latex', "location", "best");
        end
        title("t = " + num2str(plot_times(idx)) + " [s]",  'interpreter', 'latex')
        ylabel("$x$ [m]", 'interpreter', 'latex');
        set(gca, "fontsize", fontsize);
        xlim([-30, 2]);
        ylim(Xrange);
        % set(gca,'LineWidth',1);
        set(gca,'ytick',yticks);
        set(gcf, 'color', 'w');
    end

    % Figure 4, slip zone size vs. nucleation length
    fig4=figure(4);
    % set(fig1,'Position',[700 700 450 500]);      
    set(fig4, 'Units', 'inches', 'Position', [0    10    6.2500    4 * length(plot_times)]);
    
    %% plot the two lines
    % lws = linspace(3.0, 1.0, length(plot_times));
    for idx = 1:1:length(plot_times)
        subplot(length(plot_times), 1, idx);
        first_idx = find(tsaveplot > plot_times(idx));
        
        % Calculate 1 - V\theta/D_{RS}
        AwayFromSS = 1. - Vsave(:, first_idx(1)) .* thetasave(:, first_idx(1)) ./ L; 
        slipZoneIdx = find(AwayFromSS < 0.);
        this_time = tsaveplot(first_idx(1));
        distance1 = 2 * (x - x(slipZoneIdx(1))); 
        distance1(distance1 < 0.) = 0.;
        distance2 = 2 * (x(slipZoneIdx(end)) - x);
        distance2(distance2 < 0.) = 0.;
        distance = min(distance1, distance2); 
        plot(distance, x, 'linewidth', 1.5, 'color', "#0072BD");
        hold on; grid on;
        plot(L_nusave(:, first_idx(1)), x, 'linewidth', 1.5, ...
             'color', "#D95319");

        if idx == length(plot_times)
            xlabel("Length [m]",  'interpreter', 'latex');
            legend("$L_{slip}$", "$L_{nu}$", "location", "best", "interpreter", "latex")
        end
        % plot(AwayFromSS, x, 'linewidth', 1.5);

        title("t = " + num2str(plot_times(idx)) + " [s]",  'interpreter', 'latex')
        ylabel("$x$ [m]", 'interpreter', 'latex');
        set(gca, "fontsize", fontsize);
        xlim([0, 50]);
        ylim(Xrange);
        % set(gca,'LineWidth',1);
        set(gca,'ytick',yticks);
        set(gcf, 'color', 'w');
    end

    %% Save the figures
    if subtraction_flag == 0
        if saveflag == 1
            if pcflag == 1
                savename = strcat(pwd, '/../plots/', prename, '_PcVcontour_long.png');
                disp(savename);
                % saveas(figure(1),savename);
                print('-painters', '-dsvg' ,savename);
                
                savename = strcat(pwd, '/../plots/', prename, '_mfPhicontour_long.png');
                disp(savename);
                % saveas(figure(2),savename);
                print('-painters', '-dsvg' ,savename);
                
                savename = strcat(pwd, '/../plots/', prename, '_VSlipcontour_long.png');
                disp(savename);
                % saveas(figure(2),savename);
                print('-painters', '-dsvg' ,savename);
            else
                savename = strcat(pwd, '/../plots/', prename, '_PmVcontour_long.png');
                disp(savename);
                % saveas(figure(1),savename);
                print('-painters', '-dsvg' ,savename);
                
                savename = strcat(pwd, '/../plots/', prename, '_mfPhicontour_long.png');
                disp(savename);
                print('-painters', '-dsvg' ,savename);
                % saveas(figure(2),savename);
            end
        end
    else
        if saveflag == 1
            if pcflag == 1
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_PcV_LnuVsPRegion.png');
                disp(savename);
                % saveas(figure(1),savename);
                print(figure(1) ,savename, '-dpng', '-r350');
                
            else
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_PmV_LnuVsPRegion.png');
                disp(savename);
                % saveas(figure(1),savename);
                print(figure(2) ,savename, '-dpng', '-r500');
                
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_PmV_Away_logV.png');
                disp(savename);
                % saveas(figure(1),savename);
                print(figure(3) ,savename, '-dpng', '-r500');

                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_PmV_SlipZone_Lnu.png');
                disp(savename);
                % saveas(figure(1),savename);
                print(figure(4) ,savename, '-dpng', '-r500');
            end
        end
    end
end