clc,clear;
close all;

% load .mat file
filename = 'Reduced_gamma_0_pflag_3_c_4e-08.mat';
load(filename);

ind = find((tsaveplot > 1000 & tsaveplot < 1400));
figure(1);

% Save as a file
frame = 1;
for i = 1:300:size(ind, 2)
    pp = pcsave(:,ind(i));
    tt = tsaveplot(ind(i));
    theta = thetasave(:,ind(i));
    dphi_x = dphi0 - gamma * log(Vr * theta ./ L);
    L_nu = G * L / (b(1) - a(1)) / si0;
    
    subplot(2,1,1);
    plot(x, pp/si0, 'linewidth', 2.0); grid on; hold on;
    scatter(x, pp/si0, 'filled');
    xlabel('x [m]'); ylabel('\delta p_c / (\sigma_0 - p_0)');
    xlim([-2*L_nu,2*L_nu]); ylim([-0.2,1.0]);
    title(strcat('\delta p_c at t = ', num2str(tt, '%.1f'), ' s'));
    set(gca, 'fontsize', 20); hold off;

    subplot(2,1,2);
    plot(x, dphi_x, 'linewidth', 2.0); grid on; hold on;
    scatter(x, dphi_x, 'filled');
    xlabel('x [m]'); ylabel('\phi^{pl}');
    xlim([-2*L_nu,2*L_nu]);  ylim([0,5.0e-3]);
    title(strcat('\phi^{pl} at t = ', num2str(tt, '%.1f'), ' s'));
    set(gca, 'fontsize', 20); hold off;
    pause(0.001);
    
    % Save to gif
    F = getframe(gcf);
    im = frame2im(F);
    [I, map] = rgb2ind(im, 256);
    if frame == 1
        imwrite(I, map, 'p_phi.gif', 'Loopcount', inf, 'DelayTime', 0.05);
        frame = frame + 1;
    else
        imwrite(I, map, 'p_phi.gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
    
end
