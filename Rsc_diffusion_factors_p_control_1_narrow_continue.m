clc,clear;
close all;

% File prefixes
prefixs = ["temp_cbulk_e-6"]; 
% New Terminating times
term_times = [2020]; 


%for i = 1:1:1
for iii = 1:1:size(prefixs, 2)
    RC_diffusion_factors_p_control_elastic_1_narrow_continue( ...
        prefixs(iii), ...
        term_times(iii));
end