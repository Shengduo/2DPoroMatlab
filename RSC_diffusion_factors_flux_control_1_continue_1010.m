clc,clear;
close all;

% File prefixes
prefixs = ["FluxTime_0.0001_1010_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"]; 
% New Terminating times
term_times = [8080]; 


%for i = 1:1:1
for iii = 1:1:size(prefixs, 2)
    RC_diffusion_factors_flux_control_elastic_1_1010_continue( ...
        prefixs(iii), ...
        term_times(iii));
end