clc,clear;
close all;
% Find corresponding B and Cs, 
nu0 = 0.24;
nuu0 = 0.35;
B0 = 0.85;
G = 10e9;

% Injection over one [hour]
%delete(gcp('nocreate'));
%parpool(1);
% Final mass injected varies as the below array
nuus = [0.35];
% nuus = [0.262, 0.252, 0.242];

% Dilatancy coefficient gamma varies as the array below
gammas = [0];

% Bulk diffusivity
ccs = [1e-8];

% Factors of layer y-mobility
factors = [1];
% factors = [1];

% No flash heating
FHFlag = 0;

% Compare average and maximum pore pressure influence
poreflags = [6];

%for i = 1:1:1
for iii = 1:1:size(nuus, 2)
    for jjj = 1:1:size(gammas, 2)
        for kkk = 1:1:size(ccs, 2)
            for lll = 1:1:size(factors, 2)
                for mmm = 1:1:size(poreflags, 2)
                    [BC, fval] = findBCKeepKappaCmass(nuus(iii), nuu0, nu0, ccs(kkk), B0, G);
                    Regularized_cluster_cKConst_siReg_Time(nuus(iii), gammas(jjj), ccs(kkk), FHFlag, poreflags(mmm), factors(lll), BC);
                end
            end
        end
    end
end