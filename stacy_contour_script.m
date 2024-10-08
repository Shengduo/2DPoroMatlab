clc,clear;
close all;
files = ["Original_gamma_0_pflag_3_c_4e-08", ...
    "Original_gamma_0.00017_pflag_3_c_4e-08", ...
    "Original_gamma_0_pflag_3_c_4e-07", ...
    "Original_gamma_0.00017_pflag_3_c_4e-07", ...
    "Reduced_gamma_0_pflag_3_c_4e-08", ...
    "Reduced_gamma_0.00017_pflag_3_c_4e-08", ...
    "Reduced_gamma_0_pflag_3_c_4e-07", ...
    "Reduced_gamma_0.00017_pflag_3_c_4e-07", ...
    "Small_a_Original_gamma_0_pflag_3_c_4e-08", ...
    "Small_a_Original_gamma_0.00017_pflag_3_c_4e-08", ...
    "Small_a_Original_gamma_0_pflag_3_c_4e-07", ...
    "Small_a_Original_gamma_0.00017_pflag_3_c_4e-07", ...
    "Small_a_Reduced_gamma_0_pflag_3_c_4e-08", ...
    "Small_a_Reduced_gamma_0.00017_pflag_3_c_4e-08", ...
    "Small_a_Reduced_gamma_0_pflag_3_c_4e-07", ...
    "Small_a_Reduced_gamma_0.00017_pflag_3_c_4e-07", ...
    "Smaller_a_Original_gamma_0_pflag_3_c_4e-08", ...             % 17
    "Smaller_a_Original_gamma_0.00017_pflag_3_c_4e-08", ...       % 18
    "Smaller_a_Original_gamma_0_pflag_3_c_4e-07", ...             % 19
    "Smaller_a_Original_gamma_0.00017_pflag_3_c_4e-07", ...       % 20
    "Smaller_a_Reduced_gamma_0_pflag_3_c_4e-08", ...              % 21
    "Smallest_a_Reduced_gamma_0.00017_pflag_3_c_4e-08", ...       % 22
    "Smaller_a_Reduced_gamma_0_pflag_3_c_4e-07", ...              % 23
    "Smaller_a_Reduced_gamma_0.00017_pflag_3_c_4e-07", ...        % 24
    "Original_gamma_0_pflag_3_c_4e-09", ...                       % 25
    "a_0.0033_Reduced_gamma_0.00017_pflag_3_c_4e-08"...           % 26
    "Original_gamma_1.7e-05_pflag_3_c_4e-07", ...                 % 27
    "Original_gamma_1.7e-05_pflag_3_c_4e-08", ...                 % 28
    "Reduced_gamma_1.7e-05_pflag_3_c_4e-07", ...                  % 29
    "Reduced_gamma_1.7e-05_pflag_3_c_4e-08", ...                  % 30
    "1Reduced_gamma_0_pflag_3_c_4e-08", ...                       % 31
    "FH_nuu_0.262_gamma_1.7e-05_pflag_3_c_4e-08", ...             % 32
    "FH_nuu_0.262_gamma_1.7e-05_pflag_3_c_4e-07", ...             % 33
    "FH_nuu_0.35_gamma_1.7e-05_pflag_3_c_4e-07", ...              % 34
    "FH_nuu_0.35_gamma_0_pflag_3_c_4e-08", ...                    % 35
    "FH_nuu_0.262_gamma_0_pflag_3_c_4e-07", ...                   % 36
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_2e-08", ...                  % 37
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_2e-07", ...                  % 38
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_3e-07", ...                  % 39
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_3.5e-07", ...                % 40
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_3.8e-07", ...                % 41
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4.2e-07", ...                % 42
    "FH_0_nuu_0.3_gamma_0_pflag_3_c_4e-07", ...                   % 43
    "FH_0_nuu_0.32_gamma_0_pflag_3_c_4e-07", ...                  % 44
    "FH_0_nuu_0.34_gamma_0_pflag_3_c_4e-07", ...                  % 45
    "FH_0_nuu_0.31_gamma_0_pflag_3_c_4e-07", ...                  % 46
    "FH_0_nuu_0.33_gamma_0_pflag_3_c_4e-07", ...                  % 47
    "FH_0_nuu_0.335_gamma_0_pflag_3_c_4e-07", ...                 % 48
    "FH_0_nuu_0.28_gamma_0_pflag_3_c_4e-07", ...                  % 49
    "FH_0_nuu_0.332_gamma_0_pflag_3_c_4e-07", ...                 % 50
    "FH_0_nuu_0.334_gamma_0_pflag_3_c_4e-07", ...                 % 51
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_3.6e-07", ...                % 52
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_3.7e-07", ...                % 53
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-08_factor_0.01", ...      % 54
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-08_factor_0.0001", ...    % 55
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-08_factor_100", ...       % 56
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-08_factor_1e-06", ...     % 57
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_1e-06", ...     % 58
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_0.0001", ...    % 59
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_0.01", ...      % 60
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_100", ...       % 61
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_10", ...        % 62
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_25", ...        % 63
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_50", ...        % 64
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_2", ...         % 65
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_4", ...         % 66
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_6", ...         % 67
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_1", ...         % 68
    "FH_0_nuu_0.35_gamma_0.00017_pflag_3_c_4e-07_factor_1", ...         % 69
    "FH_0_nuu_0.35_gamma_1.7e-05_pflag_3_c_4e-07_factor_1", ...         % 70
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-08_factor_1", ...         % 71
    "FH_0_nuu_0.35_gamma_0.00017_pflag_3_c_4e-08_factor_1", ...         % 72
    "FH_0_nuu_0.35_gamma_1.7e-05_pflag_3_c_4e-08_factor_1", ...         % 73
    "FH_0_nuu_0.262_gamma_0_pflag_3_c_4e-07_factor_1", ...         % 74
    "FH_0_nuu_0.262_gamma_0.00017_pflag_3_c_4e-07_factor_1", ...         % 75
    "FH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_4e-07_factor_1", ...         % 76
    "FH_0_nuu_0.262_gamma_0_pflag_3_c_4e-08_factor_1", ...         % 77
    "FH_0_nuu_0.262_gamma_0.00017_pflag_3_c_4e-08_factor_1", ...         % 78
    "FH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_4e-08_factor_1", ...         % 79
    "FH_0_nuu_0.35_gamma_0_pflag_3_c_4e-07_factor_8"];            % 80
saveflag = 1;
% files = files(68:79);
% files = files([74, 68]);
files = files([77, 71]);
% files = ["NewFH_0_nuu_0.35_gamma_0_pflag_3_c_4e-08_factor_1", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_0.0004_factor_1", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_0.04_factor_1", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_4e-06_factor_1"];
files = ["NewFH_0_nuu_0.35_gamma_0_pflag_3_c_2.6982e-08_factor_1", ...
         "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_2.6982e-07_factor_1", ...
         "NewFH_0_nuu_0.262_gamma_0_pflag_3_c_3.7707e-08_factor_1", ...
         "NewFH_0_nuu_0.262_gamma_0_pflag_3_c_3.7707e-07_factor_1", ...
         "NewFH_0_nuu_0.35_gamma_1.7e-05_pflag_3_c_2.6982e-08_factor_1", ...  % 5
         "NewFH_0_nuu_0.35_gamma_1.7e-05_pflag_3_c_2.6982e-07_factor_1", ...
         "NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_3.7707e-08_factor_1", ...
         "NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_3.7707e-07_factor_1", ...
         "NewFH_0_nuu_0.35_gamma_0.00017_pflag_3_c_2.6982e-08_factor_1", ...
         "NewFH_0_nuu_0.35_gamma_0.00017_pflag_3_c_2.6982e-07_factor_1", ...  % 10
         "NewFH_0_nuu_0.262_gamma_0.00017_pflag_3_c_3.7707e-08_factor_1", ...
         "NewFH_0_nuu_0.262_gamma_0.00017_pflag_3_c_3.7707e-07_factor_1", ...
         "NewFH_0_nuu_0.35_gamma_1.7e-05_pflag_6_c_2.6982e-08_factor_1", ...
         "NewFH_0_nuu_0.35_gamma_1.7e-05_pflag_6_c_2.6982e-07_factor_1", ...
         "NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_6_c_3.7707e-08_factor_1", ...  % 15
         "NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_6_c_3.7707e-07_factor_1"];

files = ["NewFH_0_nuu_0.262_gamma_0_pflag_3_c_3.7707e-07_factor_1e-11", ...
         "NewFH_0_nuu_0.262_gamma_0_pflag_3_c_3.7707e-07_factor_1", ...
         "NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_3.7707e-07_factor_1e-11"];
     
% Keep c_mass and kappa const
files = ["NewFH_0_nuu_0.262_gamma_0_pflag_3_c_4e-08_factor_1_BC_0.43206_2.8461e-08", ...
         "NewFH_0_nuu_0.262_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.43206_2.8461e-07", ...
         "NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_4e-08_factor_1_BC_0.43206_2.8461e-08", ...
         "NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_4e-07_factor_1_BC_0.43206_2.8461e-07"];

% Second round Mar 9
files = ["NewFH_0_nuu_0.302_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.85_4e-07", ...
         "NewFH_0_nuu_0.282_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.72504_3.7831e-07", ...
         "NewFH_0_nuu_0.282_gamma_0_pflag_3_c_4e-08_factor_1_BC_0.72504_3.7831e-08", ...
         "NewFH_0_nuu_0.282_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.58305_3.0586e-07", ...
         "NewFH_0_nuu_0.302_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.68868_3.2827e-07", ...
         "NewFH_0_nuu_0.322_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.76808_3.5407e-07"];
     
files = ["NewFH_0_nuu_0.242_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.16586_3.297e-07", ...
         "NewFH_0_nuu_0.252_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.40477_3.4599e-07", ...
         "NewFH_0_nuu_0.262_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.54141_3.5747e-07", ...
         "NewFH_0_nuu_0.282_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.72504_3.7831e-07", ...
         "NewFH_0_nuu_0.302_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.85_4e-07"];
     
% files = ["NewFH_0_nuu_0.35_gamma_0_pflag_3_c_2.6982e-07_factor_1", ...
%          "NewFH_0_nuu_0.322_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.76808_3.5407e-07", ...
%          "NewFH_0_nuu_0.302_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.68868_3.2827e-07", ...
%          "NewFH_0_nuu_0.282_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.58305_3.0586e-07", ...
%          "NewFH_0_nuu_0.262_gamma_0_pflag_3_c_4e-07_factor_1_BC_0.43206_2.8461e-07"];

% files = files(4:6);

% files = ["NewFH_0_nuu_0.35_gamma_0_pflag_3_c_2.6982e-08_factor_1", ...
%          "NewFH_0_nuu_0.262_gamma_0_pflag_3_c_3.7707e-08_factor_1", ...
%          "NewFH_0_nuu_0.262_gamma_0_pflag_3_c_4e-08_factor_1_BC_0.43206_2.8461e-08"];

% files = ["NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_3_c_3.7707e-07_factor_1", ...
%          "NewFH_0_nuu_0.262_gamma_1.7e-05_pflag_6_c_3.7707e-07_factor_1"];
% % files = ["NewSiRegFH_0_nuu_0.35_gamma_0_pflag_3_c_4e-08_factor_1_BC_0.85_4e-08"];
% files = ["NewSiRegFH_0_nuu_0.35_gamma_0_pflag_6_c_4e-08_factor_1_BC_0.85_4e-08"];
% files = ["NewTimeSiRegFH_0_nuu_0.35_gamma_0_pflag_6_c_4e-08_factor_1_BC_0.85_4e-08"];
% files = ["NewFH_0_nuu_0.35_gamma_0_pflag_3_c_2.6982e-08_factor_1"];

%% Parametric study on diffusion parameters
% files = ["NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_10_10_10", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_10_1_1", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_10_1_10", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_0.1_1", ...   % 5
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_0.1_1_1", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_10_10", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_10", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_0.1", ...
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_0.1_0.1", ... % 10
%          "NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_10_1"];    % 11
% 
% files = files([2]);

%% MassControl injection
% files = ["Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "FluxTime_2.5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 2
%          "FluxTime_5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 3
%          "FluxTime_7.5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 4 
%          "Elastic_FluxTime_0.0001_NewFH_0_nuu_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 5
%          "shit4_1", ... % 6
%          "Elastic_Flag0_FluxTime_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 7
%          "Elastic1_Flag0_FluxTime_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 8
%          "Elastic_Flag2_FluxTime_0.0001_NewFH_0_nuu_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1_1"];    % 9

% % With higher initial slip rate
% files = ["HighVoFluxTime_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "HighVoFluxTime_2.5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 2
%          "HighVoFluxTime_5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 3
%          "HighVoFluxTime_7.5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"];    % 4

% % With VS rate-and-state friction
% files = ["VsFluxTime_2.5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 1
%          "VsFluxTime_5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 2
%          "VsFluxTime_7.5e-05_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ... % 3
%          "VsFluxTime_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"];    % 4
% 
% % Vs but lower a 2024-1-11
% files = ["VsFluxTime_2.5e-05_a_0.005_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "VsFluxTime_5e-05_a_0.005_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "VsFluxTime_7.5e-05_a_0.005_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "VsFluxTime_0.0001_a_0.005_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"];

% files = ["FluxTime_0.0001_1010_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "Elastic_Flag2_FluxTime_0.0001_1010_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "Elastic_Flag2_FluxTime_0.0001_1010_NewFH_0_nu_nuu_0.24_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1"];
% BigTs = [2000, 4000, 4000];
% BigTs = [8000, 4000, 8000 / 3, 2000, 2000, 2000, 2000, 2000, 2000];
% 
% selected_to_plot = [1];
% % 
% BigTs = BigTs(selected_to_plot);
% files = files(selected_to_plot);

%% William injection
% files = ["William_NewFH_verticalFlag_1_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"]; 
% BigTs = BigTs(selected_to_plot);
% files = files([1, 3]);
% BigTs = [3200];
% wideflag = 1;


%% Injection with elastic permeable solid
% files = ["Elastic_Flag2_FluxTime_2.5e-05_NewFH_0_nuu_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "Elastic_Flag2_FluxTime_5e-05_NewFH_0_nuu_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "Elastic_Flag2_FluxTime_7.5e-05_NewFH_0_nuu_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "Elastic_Flag2_FluxTime_0.0001_NewFH_0_nuu_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "Elastic_Flag2_FluxTime_0.0001_NewFH_0_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          ]; 
% 
% BigTs = [8000, 4000, 8000 / 3, 2000, 2000]; 
% 
% selected_to_plot = [1, 2, 3, 4, 5];

% files = ["Elastic_Flag2_FluxTime_0.0001_NewFH_0_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "Elastic_Flag2_FluxTime_7.5e-05_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "Elastic_Flag2_FluxTime_5e-05_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1", ...
%          "Elastic_Flag2_FluxTime_0.0001_NewFH_0_NewFH_0_nu_nuu_0.24_0.24_gamma_0_pflag_3_c_1e-08_factors_1_1_1_continue_4040", ...
%          "Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_continue_4040", ...
%          "Elastic_Flag2_FluxTime_0.0001_NewFH_0_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_continue_4040", ...
%          "FluxTime_0.0001_1010_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_continue_8080", ...
%          "shit5", ...
%          "Elastic_Flag0_FluxTime_0.0001_1010_NewFH_0_nu_nuu_0.24_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_endTime_8080"];
% BigTs = [2000, 8000 / 3, 4000, 4000, 4000, 4000, 8080, 2020, 8080]; 
% 
% selected_to_plot = [1];
% BigTs = BigTs(selected_to_plot);
% files = files(selected_to_plot);

%% Keeping cmass for elastic, permeable
% files = ["Elastic_Flag2_FluxTime_0.0001_NewFH_0_nu_nuu_0.24_0.24_gamma_0_pflag_3_c_6.1712e-09_factors_1_1_1", ...
%          "Elastic_Flag2_FluxTime_0.0001_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_6.1712e-09_factors_1_1_1"];
% BigTs = [2000, 2000]; 
% 
% selected_to_plot = [1, 2];
% BigTs = BigTs(selected_to_plot);
% files = files(selected_to_plot);

%% Keeping c for poroelastic, elastic, with narrower extent
files = ["Elastic_Flag0_pControl_NewFH_0_nu_nuu_0.24_0.35_gamma_0_pflag_3_c_1e-08_factors_0.0001_1_1", ...
         "Elastic_Flag2_pControl_NewFH_0_nu_nuu_0.24_0.24_gamma_0_pflag_3_c_1e-08_factors_0.0001_1_1", ...
         "Elastic_Flag2_pControl_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_0.0001_1_1", ...
         "Elastic_Flag2_FluxTime_0.0001_NewFH_0_nu_nuu_0.35_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_4080"];
BigTs = [2000, 2000, 2000, 4000]; 
narrow_flags = [1, 1, 1, 0];

selected_to_plot = [4];
BigTs = BigTs(selected_to_plot);
files = files(selected_to_plot);
narrow_flags = narrow_flags(selected_to_plot);

%% Fault healing parametric study
files = ["Elastic_Flag0_FluxTime_0.0001_NewFH_0_nu_nuu_0.24_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_2020_1000", ...
         "Elastic_Flag0_FluxTime_0.0001_NewFH_0_nu_nuu_0.24_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_2020_1000000", ...
         "Elastic_Flag0_FluxTime_0.0001_NewFH_0_nu_nuu_0.24_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_2020_1000000000"]; 
BigTs = [2000, 2000, 2000]; 
narrow_flags = [0, 0, 0];

selected_to_plot = [1, 2, 3];
BigTs = BigTs(selected_to_plot);
files = files(selected_to_plot);
narrow_flags = narrow_flags(selected_to_plot);

%% On fault healing
files = ["Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1"];
BigTs = [2000];
narrow_flags = [0];

selected_to_plot = [1];
BigTs = BigTs(selected_to_plot);
files = files(selected_to_plot);
narrow_flags = narrow_flags(selected_to_plot);

%% Westerly Granite
% files = ["WesterlyGranite_gamma_0_pflag_3_kappacx_8.7584e-11"];
% files = ["WesterlyGranite_Reverted_gamma_0_pflag_3_kappacx_8.7584e-15"];
for iiii = 1:1:files.length()
    % files(iiii) = strcat('../outputMats/', files(iiii));
    % stacy_contours(files(iiii), saveflag, 0, 1, BigTs(iiii), 0, narrow_flags(iiii));
    % close all;
    % % 
    % stacy_contours(files(iiii), saveflag, 0, 1, BigTs(iiii), 1, narrow_flags(iiii));
    % close all;

    stacy_LnuVsPRegion(files(iiii), saveflag, 0, 1, BigTs(iiii));
    close all;
    % % stacy_contours_undim(files(iiii), saveflag, 0, 1);
    % % PlotLeaking(files(iiii));
    % % close all;
    % center_VPres(files(iiii), saveflag, BigTs(iiii));
    % close all;
    % stacy_contours(files(iiii), saveflag, 1, 1);
    % close all;
end
