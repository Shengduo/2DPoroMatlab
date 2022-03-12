function stacy_contours_undim(prename, saveflag, pcflag, subtraction_flag)
    % pwdd = '/Users/shengduoliu/Documents/Elias_Matlab/cand_simus/stacy_replicate/Poroelastic_properties_comparison';
    %saveflag = 1;
    %pcflag = 0;

    % Filename to load
    %prename = 'Original_gamma_0_pflag_3_c_4e-08';
    %gamma = 1.7e-4;
    %Vr = 1.0e-6;
    filename = strcat('../outputMats/', prename, '.mat');
    load(filename, 'G', 'si0', 'a', 'b', 'pcsave', 'psave', 'thetasave', 'x', 'tsaveplot', ...
        'Vsave', 'dphi0', 'gamma', 'Vr', 'L', 'dsave');
    fontsize = 24;
    % Xrange
    Xrange = [-6, 6];
    xticks = 0:2:4;
    Trange = [0, 5];
    yticks = [-4, 0, 4];
    crange = [-13, 1];
    if subtraction_flag == 0
        pcsave = pcsave + 1.912e5;
        psave = psave + 1.912e5;
    end
    if pcflag == 1
        psave = pcsave;
    end
    timeend = find(tsaveplot > 1.05*2000);
    if size(timeend, 2) == 0
        timeend = size(pcsave, 2) - 1;
    else
        timeend = timeend(1);
    end
    pcsave = pcsave(:, 1:timeend);
    psave = psave(:, 1:timeend);
    tsaveplot = tsaveplot(:, 1:timeend);
    thetasave = thetasave(:, 1:timeend);
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

    % Non-dimensionalize fault length
    L_nu = G * L / (b(1) - a(1)) / si0;
    % L_nu = 1;


    % Non-dimensionalize time, time to diffuse by 1 nucleation length
    t_ = L_nu * L_nu / 0.2;
    %--------------------------------------------------------------------------
    % % First the color map
    [XX,YY] = meshgrid(tsaveplot, x);


    ha = tight_subplot(2,1,0.040,[.15 .075],[.2 .16]);
    fig1=figure(1);
    set(fig1,'Position',[700 700 450 500]);      

    %% P plot
    axes(ha(1));
    %[X, Y] = meshgrid(time,x);
    h=pcolor(XX ./ t_, YY ./ L_nu, psave/si0);
    shading interp;
    hold on; 
    set(h, 'EdgeColor', 'none')
    ylabel('$X / L_{nu}$','FontSize',fontsize, 'interpreter', 'latex')
    %colormap(ha(1),brewermap([],'REDBLUE'));
    % imagesc(peaks(250));
    caxis([-1, 1]);
    Map = bluewhitered(400);
    colormap(ha(1),Map([1:2:199 200:end],:));
    % caxis(crange);
    xlim(Trange);
    ylim(Xrange);
    c=colorbar;

    %--------------------------------------------------------------------------
    P = polyfit(mask_pc(1,:), tsaveplot, 2);
    fittime = polyval(P, mask_pc(1,:));
    %plot(fittime, mask_pc(1,:), '--k', 'linewidth', 1.5);
    plot(tsaveplot ./ t_, mask_pc(1,:) ./ L_nu, '--k', 'linewidth', 1.5);
    hold on; grid on;
    P = polyfit(mask_pc(2,:), tsaveplot, 2);
    fittime = polyval(P, mask_pc(2,:));
    % plot(fittime, mask_pc(2,:), '--k', 'linewidth', 1.5);
    plot(tsaveplot ./ t_, mask_pc(2,:) ./ L_nu, '--k', 'linewidth', 1.5);
    %--------------------------------------------------------------------------

    %plot(time(1:10:end),r_press_0_5(1:10:end),'--k','Linewidth',1.5);
    %plot(time(1:10:end),-r_press_0_5(1:10:end),'--k','Linewidth',1.5);
    set(c,'LineWidth',1);
    caxis([-0.5, 1]);
    ylabel(c,'$ \delta p_m / (\sigma_0 - p_0)$','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
    set(gca, 'TickLength', [.01 .01],...
    'TickDir','in',...
    'XMinorTick', 'on','YMinorTick', 'on','FontName',...
    'Avenir','FontSize',fontsize) 
    box on; 
    set(c, 'ylim', [-0.5, 1]);
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
    h=pcolor(X ./ t_, Y ./ L_nu, V_log);
    shading interp;
    hold on; 
    set(h, 'EdgeColor', 'none')
    ylabel('$X / L_{nu}$','FontSize',fontsize, 'interpreter', 'latex');
    colormap(ha(2),brewermap([],'YlOrRd'));

    % Mask for 0.5 MPa
    %--------------------------------------------------------------------------
    P = polyfit(mask_pc(1,:), tsaveplot, 2);
    fittime = polyval(P, mask_pc(1,:));
    %plot(fittime, mask_pc(1,:), '--k', 'linewidth', 1.5);
    plot(tsaveplot ./ t_, mask_pc(1,:) ./ L_nu, '--k', 'linewidth', 1.5);
    hold on; grid on;
    P = polyfit(mask_pc(2,:), tsaveplot, 2);
    fittime = polyval(P, mask_pc(2,:));
    % plot(fittime, mask_pc(2,:), '--k', 'linewidth', 1.5);
    plot(tsaveplot ./ t_, mask_pc(2,:) ./ L_nu, '--k', 'linewidth', 1.5);
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
    xlabel('$t / t_{nu}$','FontSize',fontsize, 'interpreter', 'latex');
    set(gca,'LineWidth',1);
    box on; 


%     %--------------------------------------------------------------------------
%     % Figure 2, dphi and slip rate
%     %--------------------------------------------------------------------------
%     fig2=figure(2);
%     ha = tight_subplot(2,1,0.040,[.15 .075],[.2 .16]);
% 
% 
% 
%     set(fig2,'Position',[100 700 450 500]);      
% 
%     %% log V plot
%     ind = [];
% 
%     % Remove negative slip rates for log scale plot:
% 
%     for i = 1:length(tsaveplot)
%          V_neg = [];
%          V_vector = Vsave(:,i);
%          V_neg = V_vector(V_vector<0);
%          if isempty(V_neg)
%              ind(end+1) = i;
%          end
%     end
% 
%     V = Vsave(:,ind);
%     time = tsaveplot(ind);
%     % r_press_0_5 = r_press_0_5(ind);
% 
%     axes(ha(1));
%     V_log=log10(V);
%     [X, Y] = meshgrid(time,x);
%     h=pcolor(X,Y,V_log);
%     shading interp;
%     hold on; 
%     set(h, 'EdgeColor', 'none')
%     ylabel('x [m]','FontSize',fontsize);
%     colormap(ha(1),brewermap([],'YlOrRd'));
% 
%     % Mask for 0.5 MPa
%     %--------------------------------------------------------------------------
%     P = polyfit(mask_pc(1,:), tsaveplot, 2);
%     fittime = polyval(P, mask_pc(1,:));
%     %plot(fittime, mask_pc(1,:), '--k', 'linewidth', 1.5);
%     plot(tsaveplot, mask_pc(1,:), '--k', 'linewidth', 1.5);
%     hold on; grid on;
%     P = polyfit(mask_pc(2,:), tsaveplot, 2);
%     fittime = polyval(P, mask_pc(2,:));
%     % plot(fittime, mask_pc(2,:), '--k', 'linewidth', 1.5);
%     plot(tsaveplot, mask_pc(2,:), '--k', 'linewidth', 1.5);
%     %--------------------------------------------------------------------------
% 
%     ylim(Xrange);
%     xlim(Trange);
%     c=colorbar;
%     set(c,'LineWidth',1);
%     caxis([-13 1]);
%     ylabel(c,'Log Slip Rate [m/s]','FontName','Avenir','FontSize',fontsize);
%     set(gca, 'TickLength', [.01 .01],...
%     'TickDir','in',...
%     'XMinorTick', 'on','YMinorTick', 'on','FontName',...
%     'Avenir','FontSize',fontsize) 
%     box on; 
%     set(gca,'layer','top')
%     set(gca,'LineWidth',1);
%     set(gca,'xtick',xticks);
%     set(gca,'ytick',yticks);
%     set(ha(1),'XTickLabel','');
%     set(gca,'Color','none');
%     box on; 
% 
% 
%     %% dphi plot
%     axes(ha(2));
%     %[X, Y] = meshgrid(time,x);
%     h=pcolor(XX,YY,dphi0 - gamma * log(Vr * thetasave ./ L));
%     shading interp;
%     hold on; 
%     set(h, 'EdgeColor', 'none')
%     ylabel('x [m]','FontSize',fontsize)
%     colormap(ha(2),brewermap([],'YlOrRd'));
%     xlim(Trange);
%     ylim(Xrange);
%     c=colorbar;
% 
%     %--------------------------------------------------------------------------
%     P = polyfit(mask_pc(1,:), tsaveplot, 2);
%     fittime = polyval(P, mask_pc(1,:));
%     %plot(fittime, mask_pc(1,:), '--k', 'linewidth', 1.5);
%     plot(tsaveplot, mask_pc(1,:), '--k', 'linewidth', 1.5);
%     hold on; grid on;
%     P = polyfit(mask_pc(2,:), tsaveplot, 2);
%     fittime = polyval(P, mask_pc(2,:));
%     % plot(fittime, mask_pc(2,:), '--k', 'linewidth', 1.5);
%     plot(tsaveplot, mask_pc(2,:), '--k', 'linewidth', 1.5);
%     %--------------------------------------------------------------------------
% 
%     %plot(time(1:10:end),r_press_0_5(1:10:end),'--k','Linewidth',1.5);
%     %plot(time(1:10:end),-r_press_0_5(1:10:end),'--k','Linewidth',1.5);
%     set(c,'LineWidth',1);
%     caxis([0, 20]);
%     %caxis([min(min(pcsave / 1.0e6)), 4 - min(min(pcsave / 1.0e6))]);
%     ylabel(c,'Dilatancy \phi^{pl}','FontName','Avenir','FontSize',fontsize);
%     set(gca, 'TickLength', [.01 .01],...
%     'TickDir','in',...
%     'XMinorTick', 'on','YMinorTick', 'on','FontName',...
%     'Avenir','FontSize',fontsize) 
%     box on; 
%     set(gca,'layer','top')
%     set(gca,'xtick',xticks);
%     set(gca,'ytick',yticks);
%     xlabel('Time [s]','FontSize',fontsize);
%     set(gca,'LineWidth',1);
%     box on; 
%     
%     %--------------------------------------------------------------------------
%     % Figure 3, slip and slip rate
%     %--------------------------------------------------------------------------
%     fig3=figure(3);
%     ha = tight_subplot(2,1,0.040,[.15 .075],[.2 .16]);
% 
%     set(fig3,'Position',[100 700 450 500]);      
% 
%     %% log V plot
%     ind = [];
% 
%     % Remove negative slip rates for log scale plot:
% 
%     for i = 1:length(tsaveplot)
%          V_neg = [];
%          V_vector = Vsave(:,i);
%          V_neg = V_vector(V_vector<0);
%          if isempty(V_neg)
%              ind(end+1) = i;
%          end
%     end
% 
%     V = Vsave(:,ind);
%     time = tsaveplot(ind);
%     % r_press_0_5 = r_press_0_5(ind);
% 
%     axes(ha(1));
%     V_log=log10(V);
%     [X, Y] = meshgrid(time,x);
%     h=pcolor(X,Y,V_log);
%     shading interp;
%     hold on; 
%     set(h, 'EdgeColor', 'none')
%     ylabel('x [m]','FontSize',fontsize);
%     colormap(ha(1),brewermap([],'YlOrRd'));
% 
%     % Mask for 0.5 MPa
%     %--------------------------------------------------------------------------
%     P = polyfit(mask_pc(1,:), tsaveplot, 2);
%     fittime = polyval(P, mask_pc(1,:));
%     %plot(fittime, mask_pc(1,:), '--k', 'linewidth', 1.5);
%     plot(tsaveplot, mask_pc(1,:), '--k', 'linewidth', 1.5);
%     hold on; grid on;
%     P = polyfit(mask_pc(2,:), tsaveplot, 2);
%     fittime = polyval(P, mask_pc(2,:));
%     % plot(fittime, mask_pc(2,:), '--k', 'linewidth', 1.5);
%     plot(tsaveplot, mask_pc(2,:), '--k', 'linewidth', 1.5);
%     %--------------------------------------------------------------------------
% 
%     ylim(Xrange);
%     xlim(Trange);
%     c=colorbar;
%     set(c,'LineWidth',1);
%     caxis([-13 1]);
%     ylabel(c,'Log Slip Rate [m/s]','FontName','Avenir','FontSize',fontsize);
%     set(gca, 'TickLength', [.01 .01],...
%     'TickDir','in',...
%     'XMinorTick', 'on','YMinorTick', 'on','FontName',...
%     'Avenir','FontSize',fontsize) 
%     box on; 
%     set(gca,'layer','top')
%     set(gca,'LineWidth',1);
%     set(gca,'xtick',xticks);
%     set(gca,'ytick',yticks);
%     set(ha(1),'XTickLabel','');
%     set(gca,'Color','none');
%     box on; 
% 
% 
%     %% Slip plot
%     axes(ha(2));
%     %[X, Y] = meshgrid(time,x);
%     h=pcolor(XX,YY,dsave(:, 1:size(thetasave, 2)));
%     shading interp;
%     hold on; 
%     set(h, 'EdgeColor', 'none')
%     ylabel('x [m]','FontSize',fontsize)
%     colormap(ha(2),brewermap([],'YlOrRd'));
%     xlim(Trange);
%     ylim([-250, 250]);
%     c=colorbar;
% 
%     %--------------------------------------------------------------------------
%     P = polyfit(mask_pc(1,:), tsaveplot, 2);
%     fittime = polyval(P, mask_pc(1,:));
%     %plot(fittime, mask_pc(1,:), '--k', 'linewidth', 1.5);
%     plot(tsaveplot, mask_pc(1,:), '--k', 'linewidth', 1.5);
%     hold on; grid on;
%     P = polyfit(mask_pc(2,:), tsaveplot, 2);
%     fittime = polyval(P, mask_pc(2,:));
%     % plot(fittime, mask_pc(2,:), '--k', 'linewidth', 1.5);
%     plot(tsaveplot, mask_pc(2,:), '--k', 'linewidth', 1.5);
%     %--------------------------------------------------------------------------
% 
%     %plot(time(1:10:end),r_press_0_5(1:10:end),'--k','Linewidth',1.5);
%     %plot(time(1:10:end),-r_press_0_5(1:10:end),'--k','Linewidth',1.5);
%     set(c,'LineWidth',1);
%     caxis([0, 20e-3]);
%     %caxis([min(min(pcsave / 1.0e6)), 4 - min(min(pcsave / 1.0e6))]);
%     ylabel(c,'Slip [m]','FontName','Avenir','FontSize',fontsize);
%     set(gca, 'TickLength', [.01 .01],...
%     'TickDir','in',...
%     'XMinorTick', 'on','YMinorTick', 'on','FontName',...
%     'Avenir','FontSize',fontsize) 
%     box on; 
%     set(gca,'layer','top')
%     set(gca,'xtick',xticks);
%     set(gca,'ytick',[-150, 0, 150]);
%     xlabel('Time [s]','FontSize',fontsize);
%     set(gca,'LineWidth',1);
%     box on; 
    if subtraction_flag == 0
        if saveflag == 1
            if pcflag == 1
                savename = strcat(pwd, '/plots/', prename, '_PcVcontour_long.png');
                disp(savename);
                % saveas(figure(1),savename);
                print('-painters', '-dsvg' ,savename);
                
                savename = strcat(pwd, '/plots/', prename, '_mVPhicontour_long.png');
                disp(savename);
                % saveas(figure(2),savename);
                print('-painters', '-dsvg' ,savename);
                
                savename = strcat(pwd, '/plots/', prename, '_VSlipcontour_long.png');
                disp(savename);
                % saveas(figure(2),savename);
                print('-painters', '-dsvg' ,savename);
            else
                savename = strcat(pwd, '/plots/', prename, '_PmVcontour_long.png');
                disp(savename);
                % saveas(figure(1),savename);
                print('-painters', '-dsvg' ,savename);
                
                savename = strcat(pwd, '/plots/', prename, '_mVPhicontour_long.png');
                disp(savename);
                print('-painters', '-dsvg' ,savename);
                % saveas(figure(2),savename);
            end
        end
    else
        if saveflag == 1
            if pcflag == 1
                savename = strcat(pwd, '../dsvg_plots_nodim1/', prename, '_PcVcontour_long.png');
                disp(savename);
                % saveas(figure(1),savename);
                print(figure(1) ,savename, '-dpng', '-r350');
                
%                 savename = strcat(pwd, '/dsvg_plots1/', prename, '_mVPhicontour_long.png');
%                 disp(savename);
%                 % saveas(figure(2),savename);
%                 print(figure(2) ,savename, '-dpng', '-r350');
%                 
%                 savename = strcat(pwd, '/dsvg_plots_nodim1/', prename, '_VSlipcontour_long.png');
%                 disp(savename);
%                 % saveas(figure(2),savename);
%                 print(figure(3) ,savename, '-dpng', '-r350');
            else
                savename = strcat(pwd, '/../dsvg_plots_nodim1/', prename, '_PmVcontour_long.png');
                disp(savename);
                % saveas(figure(1),savename);
                print(figure(1) ,savename, '-dpng', '-r500');
                
%                 savename = strcat(pwd, '/dsvg_plots_nodim1/', prename, '_mVPhicontour_long.png');
%                 disp(savename);
%                 % saveas(figure(2),savename);
%                 print(figure(2) ,savename, '-dpng', '-r500');
%                 
%                 savename = strcat(pwd, '/dsvg_plots_nodim1/', prename, '_VSlipcontour_long.png');
%                 disp(savename);
%                 % saveas(figure(2),savename);
%                 print(figure(3) ,savename, '-dpng', '-r500');
            end
        end
    end
end