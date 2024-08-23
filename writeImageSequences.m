clc,clear;
close all; 
% import the video file 
obj = VideoReader('../mp4files/Flux_0.0001_NewFH_0_nuu_0.35_gamma_0_pflag_3_c_1e-08_factors_1_1_1_logVP_withDim_wide_closeTime_0.mp4'); 
vid = read(obj); 
  
 % read the total number of frames 
frames = obj.NumFrames; 
PATH = "/Users/shengduoliu/Documents/Elias_Matlab/cand_simus/stacy_replicate/Poroelastic_properties_comparison/mp4files/image_seqs_closeTime/";

% file format of the frames to be saved in 
ST ='.jpg'; 
  
% reading and writing the frames  
for x = 1 : 1 : frames 
  
    % converting integer to string 
    Sx = num2str(x); 
  
    % concatenating 2 strings 
    Strc = strcat(PATH, Sx, ST); 
    Vid = vid(:, :, :, x); 
  
    % exporting the frames 
    imwrite(Vid, Strc); 
end