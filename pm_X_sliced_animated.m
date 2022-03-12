clc,clear;
close all;

% Array of all injection-rate cases
In_hours_ = [1, 24, 24*7, 24*30];

%for iii = 1:1:size(In_hours_, 2)
for iii = 2:1:2
    
    % load .mat file
    filename = 'P_hold_Original_gamma_0_pflag_3_c_4e-08.mat';
    load(filename);
    ind = find(tsaveplot > 1450);
    psave = psave(:, 1:ind(1));
    pcsave = pcsave(:, 1:ind(1));
    sigrsave = sigrsave(:, 1:ind(1));
    sigrnsave = 4 * (psave - 1/4 * sigrsave - 1/2 * pcsave);
    
    % Non-dimensionalize time
    t_ = L/Vr;
    
    % Non-dimensionalize fault length
    L_nu = G * L / (b(1) - a(1)) / si0;
    % 1.5* Every time
    jjj = 1;
    while jjj < size(pcsave, 2)
    %for jjj = round(0.2*size(sigrsave, 2)):1:size(sigrsave, 2)
    %for jjj = round(0.975*size(sigrsave, 2)):1:size(sigrsave, 2)
        plot(x ./ L_nu, sigrsave(:, jjj)/si0, 'linewidth', 2.5);
        
        hold on; grid on;
        plot(x ./ L_nu, pcsave(:, jjj)/si0, '--', 'linewidth', 2.0);
        plot(x ./ L_nu, sigrnsave(:, jjj)/si0, 'linewidth', 1.5);
        xlim([-10, 10]);
        ylim([-0.5, 1]);
        %ylim([-0.0002, 0.0002])
        xlabel('X / L_{nu}');
        ylabel ('\delta p / (\sigma_0 - p_0)');
        hold on; grid on;
        title(strcat('Simulated Time t =  ', num2str(tsaveplot(jjj), '%.1f'), ' s'));
        legend('\delta p^+','\delta p_c', '\delta p^-');
        set(gca, 'FontSize', 20);
        hold off;
        pause(0.1);
        if tsaveplot(jjj+1) - tsaveplot(jjj) < 1
            jjj = jjj + 20;
        elseif tsaveplot(jjj+1) - tsaveplot(jjj) > 5
            jjj = jjj + 10;
        elseif tsaveplot(jjj) < 800
            jjj = jjj + 10;
        else 
            jjj = jjj + 40;
        end
    end
    
end

%hp4 = get(subplot(2,2,4),'Position');
%c = colorbar('Position' , [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)*2.1]);

%title('log_{10}(Slip Rate)');