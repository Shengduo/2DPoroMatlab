clc,clear;
close all;
matfiles = ["Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
            "FluxTime_7.5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
            "FluxTime_5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"];

plotTs = [925.2, 1586.62, 3552.53];

for i = 1:1:length(plotTs)
    SliprateFrictionAtTimes(matfiles(i), 1, 0, 1, plotTs(i));
    close all;
end
