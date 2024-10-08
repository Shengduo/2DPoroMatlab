clc,clear;
close all;
% Injection over one [hour]
%delete(gcp('nocreate'));
%parpool(1);
% Final mass injected varies as the below array
nus = [0.24];
nuus = [0.35];

% Dilatancy coefficient gamma varies as the array below
gammas = [0];

% Bulk diffusivity
ccs = [1.0e-8];
% ccs = [6.1712e-9]; % cmass

% No flash heating
FHFlag = 0;

% All diffusion factors
kappacx_factors = [1.e-4, 1.e-4];
kappacy_factors = [1.e2, 1.e4];
bulkc_factors = [1.e2, 1.e4];

% Injection flux
baseFlux = 1.0e-4;
fluxes = baseFlux .* [1.0];

% Elastic flag, 1 - elastic, 0 - poroelastic, 2 - half-poroelastic
Elastic_Flag = 0;

% Terminating_time
Terminating_times = [2020];

%for i = 1:1:1
for iii = 1:1:size(nuus, 2)
    for jjj = 1:1:size(gammas, 2)
        for kkk = 1:1:size(ccs, 2)
            for lll = 1:1:size(kappacx_factors, 2)
                % for mmm = 1:1:size(kappacy_factors, 2)
                    % for nnn = 1:1:size(bulkc_factors, 2)
                        for fff = 1:1:size(fluxes, 2)
                            factors = [kappacx_factors(lll), ...
                                       kappacy_factors(lll), ...
                                       bulkc_factors(lll)];
                            Rc_diffusion_factors_p_control_elastic_1_narrow( nuus(iii), nus(iii), ...
                                gammas(jjj), ccs(kkk), FHFlag, 3, factors, ...
                                fluxes(fff), Elastic_Flag); 
                            % RC_diffusion_factors_flux_control_elastic_1_1010( ...
                                %, ...
                                % Terminating_times(iii));
                        end
                    % end
                % end
            end
        end
    end
end