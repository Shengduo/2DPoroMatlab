function stacy_contours(prename, saveflag, pcflag, subtraction_flag, BigT, wideflag, narrow_flag)
    % pwdd = '/Users/shengduoliu/Documents/Elias_Matlab/cand_simus/stacy_replicate/Poroelastic_properties_comparison';
    %saveflag = 1;
    %pcflag = 0;

    % Filename to load
    %prename = 'Original_gamma_0_pflag_3_c_4e-08';
    %gamma = 1.7e-4;
    %Vr = 1.0e-6;
    filename = strcat('../outputMats/', prename, '.mat');
    load(filename, 'pcsave', 'psave', 'thetasave', 'x', 'tsaveplot', ...
        'Vsave', 'dphi0', 'gamma', 'Vr', 'L', 'dsave', 'tauS', 'sisave');
    fontsize = 24;
    
    % A few parameters
    f0 = 0.5375;
    si0 = 4e6;
    tau0 = f0*si0;    
    
    % Xrange
    % Xrange = [-5, 5];
    % yticks = [-3, 0, 3];
    if wideflag == 0
        if narrow_flag == 1
            Xrange = [-0.50, 0.50];
            yticks = [-0.30, 0, 0.30];
        else
            Xrange = [-50, 50];
            yticks = [-30, 0, 30];
        end
    else
        if narrow_flag == 1
            Xrange = [-2, 2];
            yticks = [-1.2, 0, 1.2];
        else
            Xrange = [-200, 200];
            yticks = [-120, 0, 120];
        end
        prename = prename + "_wide_";
    end
    Vcrange = [-13, 1]; % [-22, -13]; % [-13, 1];
    % Vcrange = [-7, 1];
    % xticks = 0:800:2400; % 
    xticks = 0:500 * (floor(BigT / 2000)):(BigT - 500); % xticks = 0:500:1500;
    Trange = [0, BigT]; % Trange = [0, 2000];
    crange = [-2, 4]; %[-0.6, 1.2]; %[-2, 4];
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
    pc_ = 0.5e6; % 0.05e6; %0.5e6;
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
    % % First the color map
    [XX,YY] = meshgrid(tsaveplot, x);


    ha = tight_subplot(2,1,0.040,[.15 .075],[.2 .16]);
    fig1=figure(1);
    % set(fig1,'Position',[700 700 450 500]);      
    set(fig1, 'Units', 'inches', 'Position', [0   10    6.2500    6.9444]);
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
    caxis([-1.2, 1.2]); %caxis([-4,4]);
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
    clim(Vcrange);
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


    %--------------------------------------------------------------------------
    % Figure 2, dphi and slip rate
    %--------------------------------------------------------------------------
    fig2=figure(2);
    ha = tight_subplot(2,1,0.040,[.15 .075],[.2 .16]);



    % set(fig2,'Position',[100 700 450 500]);      
    set(fig2, 'Units', 'inches', 'Position', [0   10    6.2500    6.9444]);
    %% Friction coefficient plot
    axes(ha(1));
    % Compute friction coefficient: f = \tau / (\sigma - p);
    friction_coeff = tauS() ./ sisave; 
    [X, Y] = meshgrid(tsaveplot,x);
    h=pcolor(X,Y,friction_coeff);
    shading interp;
    hold on; 
    set(h, 'EdgeColor', 'none')
    ylabel('x [m]','FontSize',fontsize, 'interpreter', 'latex');
    colormap(ha(1),brewermap([],'YlOrRd'));

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
    caxis([0, 1]);
    ylabel(c,'Friction Coefficient','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
    set(gca, 'TickLength', [.01 .01],...
    'TickDir','in',...
    'XMinorTick', 'on','YMinorTick', 'on','FontName',...
    'Avenir','FontSize',fontsize) 
    box on; 
    set(gca,'layer','top')
    set(gca,'LineWidth',1);
    set(gca,'xtick',xticks);
    set(gca,'ytick',yticks);
    set(ha(1),'XTickLabel','');
    set(gca,'Color','none');
    box on; 


    %% dphi plot
    axes(ha(2));
    %[X, Y] = meshgrid(time,x);
    h=pcolor(XX,YY,1e3 * (dphi0 - gamma * log(Vr * thetasave ./ L)));
    shading interp;
    hold on; 
    set(h, 'EdgeColor', 'none')
    ylabel('x [m]','FontSize',fontsize, 'interpreter', 'latex');
    % colormap(ha(2),brewermap([],'YlOrRd'));
    caxis([-4, 4]);
    Map = bluewhitered(400);
    colormap(ha(2),Map([1:2:199 200:end],:));
    caxis([-2, 4]);
    xlim(Trange);
    ylim(Xrange);
    c=colorbar;
    % c.Ruler.TickLabelFormat = '%.3f';
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

    set(c,'LineWidth',1);
    % caxis([0, 5e-3]);
    %caxis([min(min(pcsave / 1.0e6)), 4 - min(min(pcsave / 1.0e6))]);
    ylabel(c,'$\delta \phi^{pl}\times 10^{-3}$','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
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
    
    
    
    %--------------------------------------------------------------------------
    % Figure 3, slip and slip rate
    %--------------------------------------------------------------------------
    fig3=figure(3);
    ha = tight_subplot(2,1,0.040,[.15 .075],[.2 .16]);

    % set(fig3,'Position',[100 700 450 500]);      
    set(fig3, 'Units', 'inches', 'Position', [0   10    6.2500    6.9444]);
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

    axes(ha(1));
    V_log=log10(V);
    [X, Y] = meshgrid(time,x);
    h=pcolor(X,Y,V_log);
    shading interp;
    hold on; 
    set(h, 'EdgeColor', 'none')
    ylabel('x [m]','FontSize',fontsize, 'interpreter', 'latex');
    colormap(ha(1),brewermap([],'YlOrRd'));

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
    clim(Vcrange);
    ylabel(c,'Log Slip Rate [m/s]','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
    set(gca, 'TickLength', [.01 .01],...
    'TickDir','in',...
    'XMinorTick', 'on','YMinorTick', 'on','FontName',...
    'Avenir','FontSize',fontsize) 
    box on; 
    set(gca,'layer','top')
    set(gca,'LineWidth',1);
    set(gca,'xtick',xticks);
    set(gca,'ytick',yticks);
    set(ha(1),'XTickLabel','');
    set(gca,'Color','none');
    box on; 


    %% Slip plot
    axes(ha(2));
    %[X, Y] = meshgrid(time,x);
    h=pcolor(XX,YY,dsave(:, 1:size(thetasave, 2)));
    shading interp;
    hold on; 
    set(h, 'EdgeColor', 'none')
    ylabel('x [m]','FontSize',fontsize, 'interpreter', 'latex');
    colormap(ha(2),brewermap([],'YlOrRd'));
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
    caxis([0, 20e-3]);
    %caxis([min(min(pcsave / 1.0e6)), 4 - min(min(pcsave / 1.0e6))]);
    ylabel(c,'Slip [m]','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
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
    
    
    %--------------------------------------------------------------------------
    % Figure 4, shear and total normal stress
    %--------------------------------------------------------------------------
    fig4=figure(4);
    ha = tight_subplot(2,1,0.040,[.15 .075],[.2 .16]);
    % set(fig4,'Position',[100 700 450 500]);      
    set(fig4, 'Units', 'inches', 'Position', [0   10    6.2500    6.9444]);
    %% shear stress plot
    axes(ha(1));
    % Compute shear stress
    % friction_coeff = tauS ./ sisave; 
    [X, Y] = meshgrid(time,x);
    h=pcolor(X,Y,tauS(:, ind) ./ 1e6);
    shading interp;
    hold on; 
    
    caxis([tau0 / 1.0e6 - 1, tau0 / 1.0e6 + 1]);
    c = colorbar;
    % Map = bluewhitered(400);
    colormap(ha(1),Map);
    crange = [tau0 / 1.0e6 - 1, tau0 / 1.0e6 + 1];
    caxis(crange);
    
    set(h, 'EdgeColor', 'none')
    ylabel('x [m]','FontSize',fontsize, 'interpreter', 'latex');
    % colormap(ha(1),brewermap([],'YlOrRd'));

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
    ylabel(c,'Shear Stress [MPa]','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
    set(gca, 'TickLength', [.01 .01],...
    'TickDir','in',...
    'XMinorTick', 'on','YMinorTick', 'on','FontName',...
    'Avenir','FontSize',fontsize) 
    box on; 
    set(gca,'layer','top')
    set(gca,'LineWidth',1);
    set(gca,'xtick',xticks);
    set(gca,'ytick',yticks);
    set(ha(1),'XTickLabel','');
    set(gca,'Color','none');
    box on; 


    %% normal stress plot
    axes(ha(2));
    %[X, Y] = meshgrid(time,x);
    h=pcolor(XX,YY, sisave ./ 1.0e6);
    shading interp;
    hold on; 
    set(h, 'EdgeColor', 'none')
    ylabel('x [m]','FontSize',fontsize, 'interpreter', 'latex');
    % colormap(ha(2),brewermap([],'YlOrRd'));
    caxis([si0 / 1e6 - 4, si0 / 1e6 + 4]);
    crange = [si0 / 1e6 - 4, si0 / 1e6];
    % Map = bluewhitered(400);
    colormap(ha(2),Map(1:200,:));
    caxis(crange);
    xlim(Trange);
    ylim(Xrange);
    c=colorbar;
    % c.Ruler.TickLabelFormat = '%.3f';
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

    set(c,'LineWidth',1);
    % caxis([0, 5e-3]);
    %caxis([min(min(pcsave / 1.0e6)), 4 - min(min(pcsave / 1.0e6))]);
    ylabel(c,{'Normal Stress',  '[MPa]'}, 'FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
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
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_PcVcontour_long.png');
                disp(savename);
                % saveas(figure(1),savename);
                print(figure(1) ,savename, '-dpng', '-r350');
                
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_mfPhicontour_long.png');
                disp(savename);
                % saveas(figure(2),savename);
                print(figure(2) ,savename, '-dpng', '-r350');
                
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_VSlipcontour_long.png');
                disp(savename);
                % saveas(figure(2),savename);
                print(figure(3) ,savename, '-dpng', '-r350');
                
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_TauNcontour_long.png');
                disp(savename);
                % saveas(figure(2),savename);
                print(figure(3) ,savename, '-dpng', '-r350');
            else
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_PmVcontour_long.png');
                disp(savename);
                % saveas(figure(1),savename);
                print(figure(1) ,savename, '-dpng', '-r500');
                
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_mfPhicontour_long.png');
                disp(savename);
                % saveas(figure(2),savename);
                print(figure(2) ,savename, '-dpng', '-r500');
                
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_VSlipcontour_long.png');
                disp(savename);
                % saveas(figure(2),savename);
                print(figure(3) ,savename, '-dpng', '-r500');
                
                savename = strcat(pwd, '/../dsvg_plots2/', prename, '_TauNcontour_long.png');
                disp(savename);
                % saveas(figure(2),savename);
                print(figure(4) ,savename, '-dpng', '-r500');
            end
        end
    end
end