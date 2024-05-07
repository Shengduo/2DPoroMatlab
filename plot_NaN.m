% clc,clear; 
close all;
% load("../outputMats/shit5.mat");
fsize = 20;


for i = 1:1:5
    fig = figure(); 
    set(fig, 'Units', 'inches', 'Position', [0    10    7.7778* 3    5.8333]);
    subplot(1, 4, 1);
    semilogy(x, Vsave(:, runnerplot - i), 'linewidth', 2.0);
    xlabel("$x$ [m]", 'fontsize', fsize, 'Interpreter', 'latex');
    ylabel("$\log(V)$", 'fontsize', fsize, 'Interpreter', 'latex');
    set(gca, "fontsize", fsize); 
    
    subplot(1, 4, 2)
    plot(x, dsave(:, runnerplot - i), 'linewidth', 2.0);
    xlabel("$x$ [m]", 'fontsize', fsize, 'Interpreter', 'latex');
    ylabel("slip [m]", 'fontsize', fsize, 'Interpreter', 'latex');
    set(gca, "fontsize", fsize); 

    subplot(1, 4, 3)
    plot(x, sisave(:, runnerplot - i) / 1.e6, 'linewidth', 2.0);
    xlabel("$x$ [m]", 'fontsize', fsize, 'Interpreter', 'latex');
    ylabel("$p_m$ [MPa]", 'fontsize', fsize, 'Interpreter', 'latex');
    set(gca, "fontsize", fsize); 

    subplot(1, 4, 4)
    plot(x, tauS(:, runnerplot - i) / 1.e6, 'linewidth', 2.0);
    xlabel("$x$ [m]", 'fontsize', fsize, 'Interpreter', 'latex');
    ylabel("$\tau$ [MPa]", 'fontsize', fsize, 'Interpreter', 'latex');

    set(gcf, "color", "white");
    set(gca, "fontsize", fsize); 
    
    sgtitle("t =\ " + num2str(tsaveplot(runnerplot - i)) + " [s]", ...
        'fontsize', fsize, 'Interpreter', 'latex');
    print(fig , "../plots/shit5_runner_" + num2str(i) + ".png", '-dpng', '-r350');
end