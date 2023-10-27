clc,clear;
close all;
% Injection over one [hour]
%delete(gcp('nocreate'));
%parpool(1);
% Final mass injected varies as the below array
nuus = [0.35];

% Dilatancy coefficient gamma varies as the array below
gammas = [0];

% Bulk diffusivity
ccs = [1.0e-8];

% No flash heating
FHFlag = 0;

% All diffusion factors
kappacx_factors = [1.];
kappacy_factors = [1.];
bulkc_factors = [1.];

% Injection flux
baseFlux = 1.0e-4;
fluxes = baseFlux .* [0.5];


%for i = 1:1:1
for iii = 1:1:size(nuus, 2)
    for jjj = 1:1:size(gammas, 2)
        for kkk = 1:1:size(ccs, 2)
            for lll = 1:1:size(kappacx_factors, 2)
                for mmm = 1:1:size(kappacy_factors, 2)
                    for nnn = 1:1:size(bulkc_factors, 2)
                        for fff = 1:1:size(fluxes, 2)
                            factors = [kappacx_factors(lll), ...
                                       kappacy_factors(mmm), ...
                                       bulkc_factors(nnn)];
                            Regularized_cluster_diffusion_factors_flux_control(nuus(iii), ...
                                gammas(jjj), ccs(kkk), FHFlag, 3, factors, ...
                                fluxes(fff));
                        end
                    end
                end
            end
        end
    end
end