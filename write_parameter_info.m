clc,clear;
close all;
% Specify the filename
gamma = 0;
poreflag = 3;

% Filename to load
prename = 'Reduced_gamma_0.00017_pflag_3_c_4e-08';
filename = strcat(prename, '.mat');
load(filename);

% Write changable parameters into a '.txt' file
txtname = [prename, '_params.txt'];

fileID = fopen(txtname, 'w');

fprintf(fileID, 'Filename: '); fprintf(fileID, filename);
% fprintf(fileID, '\n%40s', 'Injection Mass (Kg): '); fprintf(fileID, num2str(in_mass));
% fprintf(fileID, '\n%40s', 'NT: '); fprintf(fileID, num2str(NT));
fprintf(fileID, '\n---------------------------------------------------------------------------'); 
fprintf(fileID, '\nInitial Condition: '); 
fprintf(fileID, '\n%40s', 'Initial Pre-stress [MPa]: '); fprintf(fileID, num2str(si0 / 1.0e6));
fprintf(fileID, '\n%40s', 'Initial Friction Coefficient: '); fprintf(fileID, num2str(f0));
fprintf(fileID, '\n%40s', 'Initial Friction [MPa]: '); fprintf(fileID, num2str(tau0 / 1.0e6));
fprintf(fileID, '\n%40s', 'Initial Slip Rate [m/s]: '); fprintf(fileID, num2str(Vo));
fprintf(fileID, '\n%40s', 'Initial State Variable: '); fprintf(fileID, num2str(theta_0));

fprintf(fileID, '\n---------------------------------------------------------------------------'); 
fprintf(fileID, '\nReference Friction Parameters: '); 
fprintf(fileID, '\n%40s', 'Reference Friction Coefficient: '); fprintf(fileID, num2str(fr));
fprintf(fileID, '\n%40s', 'Reference Slip Rate [m/s]: '); fprintf(fileID, num2str(Vr));

fprintf(fileID, '\n---------------------------------------------------------------------------'); 
fprintf(fileID, '\nShear Layer Parameters: '); 
fprintf(fileID, '\n%40s', 'Injection Interval [m]: '); 
fprintf(fileID, strcat(num2str(lhs), ', ', num2str(rhs)));
fprintf(fileID, '\n%40s', 'Fault Width 2h [mm]: '); 
fprintf(fileID, num2str(2 * epsi * 1000)); 
fprintf(fileID, '\n%40s', 'Dilatancy Coefficient Gamma: '); fprintf(fileID, num2str(gamma));
fprintf(fileID, '\n%40s', 'Rate and State a, b: '); 
fprintf(fileID, strcat(num2str(a(1)), ', ', num2str(b(1)))); 
fprintf(fileID, '\n%40s', 'Rate and State D_RS [m]: '); 
fprintf(fileID, num2str(L)); 
fprintf(fileID, '\n%40s', 'Fluid Mobility-y kappac [m^3 s/Kg]: '); fprintf(fileID, num2str(kappac));
fprintf(fileID, '\n%40s', 'Fluid Mobility-x kappacx [m^3 s/Kg]: '); fprintf(fileID, num2str(kappacx));

fprintf(fileID, '\n---------------------------------------------------------------------------'); 
fprintf(fileID, '\nBulk Parameters: '); 
fprintf(fileID, '\n%40s', 'Bulk Diffusivity c [m^2/s]: '); fprintf(fileID, num2str(c));
fprintf(fileID, '\n%40s', 'Drained Shear Modulus G [GPa]: '); fprintf(fileID, num2str(G/1.0e9));
K = 2 * G * (1 + nu) / (3 * (1 - 2 * nu));
fprintf(fileID, '\n%40s', "Drained Bulk Modulus K [GPa]: "); fprintf(fileID, num2str(K/1.0e9));
fprintf(fileID, '\n%40s', "Drained Poisson's Ratio nu: "); fprintf(fileID, num2str(nu));
K_u = 2 * G * (1 + nuu) / (3 * (1 - 2 * nuu));
fprintf(fileID, '\n%40s', "Undrained Bulk Modulus K_u [GPa]: "); fprintf(fileID, num2str(K_u/1.0e9));
fprintf(fileID, '\n%40s', "Undrained Poisson's Ratio nu_u: "); fprintf(fileID, num2str(nuu));
fprintf(fileID, '\n%40s', "Biot Coefficient alpha: "); fprintf(fileID, num2str(alpB));
fprintf(fileID, '\n%40s', "Biot Modulus M [GPa]: "); fprintf(fileID, num2str((K_u-K)/alpB^2/1.0e9));
fprintf(fileID, '\n%40s', "Skempton Coefficient B: "); fprintf(fileID, num2str(B));
fprintf(fileID, '\n%40s', "Bulk Fluid Mobility kappa: "); fprintf(fileID, num2str(kappa));
fclose(fileID);