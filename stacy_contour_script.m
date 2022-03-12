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
files = files([3, 7, 43]);
for iiii = 1:1:files.length()
    % stacy_contours(files(iiii), saveflag, 0, 1);
    stacy_contours_undim(files(iiii), saveflag, 0, 1);
    close all;
    center_VPres_nodim(files(iiii), saveflag);
    close all;
    % stacy_contours(files(iiii), saveflag, 1, 1);
    % close all;
end
