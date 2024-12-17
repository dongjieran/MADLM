% _______________________________________________________________________
% main_MADLM.m
% version 1 (July 1, 2024)
% Required subroutines: MADLM.m, DLMdiff.m,
%                       interp2qq.m, dataspecDLM.m, tav_list.txt,
%                       MADLM_forward.m
% _______________________________________________________________________

% This script simulates the multi-angular bidirectional reflectance factor (BRF) of the dorsiventral leaf
% and retrieves both leaf biochemical and surface structural parameters from the BRF of adaxial and abaxial sides.
%
% The diffuse component (400-2500 nm) with a 1 nm step is calculated using the DLM model (Stuckens et al., 2009).
% Reference:
% Stuckens J, Verstraeten W W, Delalieux S, et al. A dorsiventral leaf radiative transfer model: Development, 
% validation, and improved model inversion techniques. Remote Sensing of Environment, 2009, 113(12): 2560-2573.
% The MATLAB code for the DLM model is available at:
% https://artmotoolbox.com/radiative-transfer-models/84-rtm-leaf/13-dlm-model.html
%
% The specular component is calculated based on the specular reflection function (Bousquet et al., 2005).
% Reference:
% Bousquet, L., Lach¨¦rade, S., Jacquemoud, S., & Moya, I., 2005. Leaf BRDF measurements and model
% for specular and diffuse components differentiation. Remote Sensing of Environment, 98, 201-211.
% _______________________________________________________________________
% MADLM
% Dongjie Ran, Zhongqiu Sun, Shan Lu, Kenji Omasa, 2024 
% An advanced dorsiventral leaf radiative transfer model for simulating multi-angular
% and spectral reflection: considering asymmetry of leaf internal and surface structure.Remote Sensing of Environment,2025, 318: 114531.
% _______________________________________________________________________

clc;clear;
%% Spectral simulation

% Load structure and biochemistry parameters
load('leaf_parameter.txt'); % Leaf parameters

% Parameters description:
% Cab        = leaf_parameter(1);  % Chlorophyll a+b content in ug/cm2
% Cxc        = leaf_parameter(2);  % Carotenoids content in ug/cm2
% Cw         = leaf_parameter(3);  % Equivalent water thickness in g/cm2 or cm
% Cm         = leaf_parameter(4);  % Dry matter content in g/cm2
% Cbrown     = leaf_parameter(5);  % Brown pigments content in arbitrary units
% fair       = leaf_parameter(6);  % Structure parameter containing the proportion of air spaces between two layers.
%                                    The higher this proportion is, the more internal scattering.
% beta_pigm  = leaf_parameter(7);  % Asymmetry parameter for chlorophyll/carotenoids: the proportion of the total leaf
%                                    pigments in the palisade layer.
% beta_wdm   = leaf_parameter(8);  % Fraction of water and dry matter in the palisade layer.
% beta_abaxep= leaf_parameter(9);  % Fraction of pigments in the abaxial epidermis.
% abdiff     = leaf_parameter(10); % 1/90 of diffusion angle of light for abaxial reflectance. 1 is full diffusion;
%                                    0 is collimated light (no diffusion). This controls the differential scattering of adaxially
%                                    and abaxially incident light.
% mu         = leaf_parameter(11); % Roughness term for DHR simulation in the DLM
% uad        = leaf_parameter(12); % Roughness of adaxial leaf surface.
% uab        = leaf_parameter(13); % Roughness of abaxial leaf surface.
% fad        = leaf_parameter(14); % Fraction of refractive index for the adaxial leaf surface.
% fab        = leaf_parameter(15); % Fraction of refractive index for the abaxial leaf surface.

% Load illumination-viewing geometry
geo = load('geometry.txt'); % Load leaf parameter geometry
% SZA = geo(1,:); % Source zenith angle (SZA) in degrees 
% VZA = geo(2,:); % Viewing zenith angle (VZA) in degrees
% VAA = geo(3,:); % Viewing azimuth angle (VAA) in degrees

% Convert degrees to radians for geometric parameters
geomparams = deg2rad(geo);
geomparams = geomparams';
num_angles = size(geo, 2);

% Wavelength range
wave = (400:1:2500)'; % Wavelength range from 400 nm to 2500 nm with 1 nm step
wl = wave - 399; % Adjust to start from 1 for indexing purposes

% Simulate multi-angular BRFs for adaxial and abaxial sides
BRF = MADLM(leaf_parameter, wl, geomparams); % Call the MADLM function to simulate BRF

% Separate BRF data for adaxial and abaxial sides
BRF_adaxial = BRF(:, 1:num_angles); 
BRF_abaxial = BRF(:, num_angles+1:end); 

% Save the BRF data to an Excel file
xlswrite('leaf_spectrum.xlsx', [geo; BRF_adaxial], 'Adaxial_BRF');
xlswrite('leaf_spectrum.xlsx', [geo; BRF_abaxial], 'Abaxial_BRF');

% Plot the results
figure;
plot(wave, BRF_adaxial(:,1), 'b', wave, BRF_abaxial(:,1), 'r', wave, BRF_adaxial(:,11), 'g', wave, BRF_abaxial(:,11), 'm')
title('Multi-angular BRFs of the dorsiventral leaf')
xlabel('Wavelength (nm)')
ylabel('BRF')
axis([400 2500 0 1.2]) % Set the axis limits
legend({'Adaxial BRF (VZA=0¡ã)', 'Abaxial BRF (VZA=0¡ã)', 'Adaxial BRF (VZA=-30¡ã)', 'Abaxial BRF (VZA=-30¡ã)'})

%% Model inversion

% Identify parameter boundary and initial value
%     Cab  Cxc   Cw     Cm     Cbrown fair beta_pigm beta_wdm beta_abaxep abdiff mu    uad  uab  fad  fab
P0 = [40   10    0.01   0.005  0      0.34  0.52      0.44     0.06       0.67   1.19  0.4  0.3  1.47 1.47]; % initial value
LB = [0    0     0      0      0      0     0         0        0          0      0     0.01 0.01 0.6  0.6]; % Lower bounds
UB = [120  30    0.06   0.05   1      1     1         1        0.25       1      2     1    1    3.5  3.5]; % Upper bounds

options = optimset('Algorithm','trust-region-reflective','TolFun',1e-4); 
sol = lsqcurvefit(@MADLM, P0, wl, BRF, LB, UB, options, geomparams);

% Reconstruct BRF spectra using the retrieved parameters
[BRFsim, specular, diff] = MADLM_forward(sol, wl, geomparams);
% BRFsim: simulated total BRFs for the adaxial and abaxial sides
% specular: simulated specular components for the adaxial and abaxial sides
% diff: diffuse component for the adaxial and abaxial sides

% Output adaxial results
xlswrite('leaf_spectrum_reconstruction_adaxial.xlsx', BRFsim(:, 1:num_angles), 'Simulated_BRF');
xlswrite('leaf_spectrum_reconstruction_adaxial.xlsx', specular(:, 1:num_angles), 'specular');
xlswrite('leaf_spectrum_reconstruction_adaxial.xlsx', diff(:, 1), 'diff');

% Output abaxial results
xlswrite('leaf_spectrum_reconstruction_abaxial.xlsx', BRFsim(:, num_angles+1:end), 'Simulated_BRF');
xlswrite('leaf_spectrum_reconstruction_abaxial.xlsx', specular(:, num_angles+1:end), 'specular');
xlswrite('leaf_spectrum_reconstruction_abaxial.xlsx', diff(:, 2), 'diff');

