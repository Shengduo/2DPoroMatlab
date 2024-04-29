clc,clear;
close all;

% File prefixes
prefixs = ["Elastic_Flag2_FluxTime_0.0001_NewFH_0_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
           "Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"]; 
% New Terminating times
term_times = [4040, 4040]; 


%for i = 1:1:1
for iii = 2:1:size(prefixs, 2)
    RC_diffusion_factors_flux_control_elastic_1_continue( ...
        prefixs(iii), ...
        term_times(iii));
end