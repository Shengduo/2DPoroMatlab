clc,clear;
close all;
% Injection over one [hour]
%delete(gcp('nocreate'));
%parpool(1);

% Dilatancy coefficient gamma varies as the array below
gammas = [0];

namePrefix = "WesterlyGranite_Reverted";

% kpa_cxs = 8.7584e-11 corresponds to cx = 2e-1 m^2/s, we want 2e-5, with bulk
% c_bulk = 2.2e-3 m^2/s
kpa_cxs = [8.7584e-15]; 
% No flash heating
FHFlag = 0;


for jjj = 1:1:size(gammas, 2)
    for iii = 1:1:size(kpa_cxs, 2)
        Regularized_cluster_new_0827_reverted(gammas(jjj), FHFlag, 3, kpa_cxs(iii), namePrefix);
    end
end