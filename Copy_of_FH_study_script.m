clc,clear;
close all;
% Injection over one [hour]
%delete(gcp('nocreate'));
%parpool(1);
% Final mass injected varies as the below array
nuus = [0.262];

% Dilatancy coefficient gamma varies as the array below
gammas = [0];

% Bulk diffusivity
ccs = [4.0e-7];
%for i = 1:1:1
for iii = 1:1:size(nuus, 2)
    for jjj = 1:1:size(gammas, 2)
        for kkk = 1:1:size(ccs, 2)
            FH_regularized_cluster(nuus(iii), gammas(jjj), ccs(kkk));
        end
    end
end