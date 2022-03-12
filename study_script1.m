clc,clear;
close all;
% Injection over one [hour]
%delete(gcp('nocreate'));
%parpool(1);
% Final mass injected varies as the below array
in_mass = [0];

% Dilatancy coefficient gamma varies as the array below
gammas = [0, 1.7e-4];

% Bulk diffusivity
ccs = [4.0e-7, 4.0e-8];
%for i = 1:1:1
for iii = 1:1:size(in_mass, 2)
    for jjj = 1:1:size(gammas, 2)
        for kkk = 1:1:size(ccs, 2)
            poro_run_regularized_cluster2(in_mass(iii), gammas(jjj), ccs(kkk));
        end
    end
end