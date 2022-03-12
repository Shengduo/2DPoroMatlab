%clc,clear;
close all;

% Array of all injection-rate cases
In_hours_ = [1, 24, 24*7, 24*30];
ylims = [1e-24, 1e2];
%for iii = 1:1:size(In_hours_, 2)
for iii = 2:1:2
    % load .mat file
    filename = 'Original_gamma_0_pflag_3_c_4e-08.mat';
    % filename = strcat('In_hours_',num2str(In_hours_(iii)),'_nogamma_NS1.mat');
    %load(filename);
    slices = 1*(1:1:size(Vsave, 2));
    jetcustom = cool(size(slices, 2)); 
    Vdyn = 0.1;
    % Non-dimensionalize fault length
    L_nu = G * L / (b(1) - a(1)) / si0;
    jjj = 1;
    while max(Vsave(:,jjj+20)) < 1e-6
        semilogy(x./L_nu, Vsave(:, jjj)./Vdyn, 'linewidth', 1.5, 'color', jetcustom(jjj, :));
        ylim(ylims);
        xlim([-260/L_nu, 260/L_nu]);
        xlabel('X / L_{nu}');
        ylabel ('V / V_{dyn}');
        hold on; grid on;
        title(strcat('Simulated Time t =  ', num2str(tsaveplot(jjj), '%.2f'), 's'));
        set(gca, 'FontSize', 15);
        pause(0.4);
        jjj = jjj + 300;
    end
    
    while jjj < size(Vsave, 2) && max(Vsave(:,jjj)) < 0.1
        semilogy(x./L_nu, Vsave(:, jjj)./Vdyn, 'linewidth', 1.5, 'color', jetcustom(jjj, :));
        ylim(ylims);
        xlim([-260/L_nu, 260/L_nu]);
        xlabel('X / L_{nu}');
        ylabel ('V / V_{dyn}');
        hold on; grid on;
        title(strcat('Simulated Time t =  ', num2str(tsaveplot(jjj), '%.2f'), 's'));
        set(gca, 'FontSize', 15);
        pause(0.4);
        jjj = jjj + 3000;
    end
    
    while jjj < size(Vsave, 2)
        semilogy(x./L_nu, Vsave(:, jjj)./Vdyn, 'linewidth', 1.5, 'color', jetcustom(jjj, :));
        ylim(ylims);
        xlim([-260/L_nu, 260/L_nu]);
        xlabel('X / L_{nu}');
        ylabel ('V / V_{dyn}');
        hold on; grid on;
        title(strcat('Simulated Time t =  ', num2str(tsaveplot(jjj), '%.2f'), 's'));
        set(gca, 'FontSize', 15);
        pause(0.4);
        jjj = jjj + 1000;
    end
    
    colormap(jet(size(slices, 2)));
    cb = colorbar;
    caxis([slices(1), slices(end)]);
    ylabel(cb,'Time Step Number');
    set(gca, 'FontSize', 15);
end

%hp4 = get(subplot(2,2,4),'Position');
%c = colorbar('Position' , [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)*2.1]);

%title('log_{10}(Slip Rate)');