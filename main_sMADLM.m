% _______________________________________________________________________
% main_sMADLM.m
% version 1 (July 1, 2024)
% Required subroutines: sMADLM.m, DLMdiff.m,
%                       interp2qq.m, dataspecDLM.m, tav_list.txt
% _______________________________________________________________________

% This model allows to retrieve biochemical parameters for one dorsiventral leaf at individual vewing angle 
% or for a single pixel based on the DLM model (Stuckens et al., 2009) and a simplified specular reflection term Rspec.
%
% The diffuse component (400-2500 nm) with a 1 nm step is calculated using the DLM model (Stuckens et al., 2009).
% Reference:
% Stuckens J, Verstraeten W W, Delalieux S, et al. A dorsiventral leaf radiative transfer model: Development, 
% validation, and improved model inversion techniques. Remote Sensing of Environment, 2009, 113(12): 2560-2573.
% The MATLAB code for the DLM model is available at:
% https://artmotoolbox.com/radiative-transfer-models/84-rtm-leaf/13-dlm-model.html
% Inverting the sMADLM model allows the estimation of the
% following 13 parameters:
%
%    - Cab: Chlorophyll a+b content in ¦Ìg/cm2
%    - Cxc: Carotenoids content in ¦Ìg/cm2
%    - Cw: Equivalent water thickness in g/cm2 or cm
%    - Cm: Dry matter content in g/cm2
%    - Cbrown: Brown pigments content in arbitrary units
%    - fair: Airspaces between layers
%    - beta_pigm: Fraction of total pigments in the palisade layer
%    - beta_wdm: Fraction of total water and dry matter in the palisade layer
%    - beta_abaxep: Fraction of total pigments in the abaxial epidermis
%    - abdiff: Abaxial scattering coefficient
%    - mu: Roughness term for DHR simulation in the DLM
%    - vad: The wavelength-independent term that depends on viewing geometry and roughness for BRF simulation 
%           in the sMADLM for the adaxial leaf surfaces
%    - vab: The wavelength-independent term that depends on viewing geometry and roughness for BRF simulation 
%           in the sMADLM for the abaxial leaf surfaces
% _______________________________________________________________________
% sMADLM: a simplified version of MADLM
% Dongjie Ran, Zhongqiu Sun, Shan Lu, Kenji Omasa, 2024 
% An advanced dorsiventral leaf radiative transfer model for simulating multi-angular
% and spectral reflection: considering asymmetry of leaf internal and surface structure.Remote Sensing of Environment.
% _______________________________________________________________________

clc; clear;
%% Model simulation

% Load leaf parameters
leaf_parameters = load('leaf_parameter_sMADLM.txt'); % Leaf parameters

% Wavelength range
wave = (400:1:2500)'; % Wavelength range from 400 nm to 2500 nm with 1 nm step
wl = wave - 399; % Adjust to start from 1 for indexing purposes

% Simulate BRF using sMADLM model
data=dataspecDLM;
BRF = sMADLM(leaf_parameters, wl,data);

% Plot the results
figure;
plot(wave, BRF(:,1), 'r', 'LineWidth', 1.5);
hold on;
plot(wave, BRF(:,2), 'b', 'LineWidth', 1.5);
xlabel('Wavelength (nm)', 'FontSize', 11, 'FontWeight', 'Bold');
ylabel('Bidirectional Reflectance Factor', 'FontSize', 12, 'FontWeight', 'Bold');
title(['C_{ab} = ', num2str(leaf_parameters(1), '%4.2f'), ...
       '  C_{xc} = ', num2str(leaf_parameters(2), '%4.1f'), ...
       '  C_{w} = ', num2str(leaf_parameters(3), '%4.2f'), ...
       '  C_{m} = ', num2str(leaf_parameters(4), '%4.1f'), ...
       '  C_{br} = ', num2str(leaf_parameters(5), '%7.4f'), ...
       '  f_{air} = ', num2str(leaf_parameters(6), '%7.4f'), ...
       '  \beta_{pigm} = ', num2str(leaf_parameters(7), '%4.1f'), ...
       '  \beta_{wdm} = ', num2str(leaf_parameters(8), '%4.2f'), ...
       '  \beta_{abaxep} = ', num2str(leaf_parameters(9), '%4.1f'), ...
       '  \delta = ', num2str(leaf_parameters(10), '%4.2f'), ...
       '  \mu = ', num2str(leaf_parameters(11), '%4.2f'), ...
       '  v_{ad} = ', num2str(leaf_parameters(12), '%4.1f'), ...
       '  v_{ab} = ', num2str(leaf_parameters(13), '%4.2f')], ...
       'FontWeight', 'Bold');
axis([400 2500 0 1]);
legend('Adaxial BRF', 'Abaxial BRF'); 
%%
clc;clear;
% Load measured close-range hyperspectral image data
load('Populus_adaxial.mat'); % Load adaxial leaf data
load('Populus_abaxial.mat'); % Load abaxial leaf data

% Define the wavelength range for close-range hyperspectral images (450-800 nm)
range=(39:205);
wl=wave(range,:);

data=dataspecDLM;
data = spline(data(:,1),data',wl)'; % Resampling to the actual spectral sampling interval of the measurement

% Display the RGB images of the adaxial and abaxial sides
rgb = [149 87 38];
figure,imshow(image_ad(:,:,rgb)*3);
figure,imshow(image_ab(:,:,rgb)*3);

% Extract and smooth spectra at a specific pixel (110, 120)
spectrum_ad= squeeze(image_ad(230, 160, :)); 
spectrum_ab= squeeze(image_ab(230, 160, :));
% Smooth spectra using Savitzky-Golay filter
spectrum_ad_SG=sgolayfilt(spectrum_ad, 2, 19);
spectrum_ab_SG=sgolayfilt(spectrum_ab, 2, 19);

BRFmea=[spectrum_ad_SG(range,:) spectrum_ab_SG(range,:)];

%   Cab Cxc   Cw    Cm     Cbrown   fair  beta_pigm beta_wdm beta_abaxep abdiff mu      vad     vab
P0=[40    10    0.01  0.005  0      0.34  0.52      0.44     0.06         0.67  1.19    1.19    1.19]; % initial value
LB=[0     0     0.01  0.005  0      0     0         0        0            0     0       0       0];    % Lower bounds
UB=[120   30    0.01  0.005  1      1     1         1        0.25         1     2       30      30];   % Upper bounds

options = optimset('Algorithm','trust-region-reflective','TolFun',1e-4);
sol= lsqcurvefit(@sMADLM,P0,wl,BRFmea,LB,UB,options,data);

BRFsim=sMADLM(sol,wl,data);

figure
p=plot(wl, BRFmea(:,1), 'r', wl, BRFsim(:,1), 'r--', wl, BRFmea(:,2), 'b', wl, BRFsim(:,2), 'b--');
set(p(1),'LineWidth',2.5,'Color',[0.75,0.75,0.99])
set(p(2),'LineWidth',1.5,'Color',[0,0,0.7],'LineStyle',':')
set(p(3),'LineWidth',2.5,'Color',[0.99, 0.5, 0.5])
set(p(4),'LineWidth',1.5,'Color',[0.7, 0, 0],'LineStyle',':')
xlabel('Wavelength (nm)','FontSize',11,'FontWeight','Bold')
ylabel('Bidirectional Reflectance Factor','FontSize',12,'FontWeight','Bold')
title(['C_{ab} = ', num2str(sol(1), '%4.2f'), ...
       '  C_{xc} = ', num2str(sol(2), '%4.1f'), ...
       '  C_{w} = ', num2str(sol(3), '%4.2f'), ...
       '  C_{m} = ', num2str(sol(4), '%4.1f'), ...
       '  C_{br} = ', num2str(sol(5), '%7.4f'), ...
       '  f_{air} = ', num2str(sol(6), '%7.4f'), ...
       '  \beta_{pigm} = ', num2str(sol(7), '%4.1f'), ...
       '  \beta_{wdm} = ', num2str(sol(8), '%4.2f'), ...
       '  \beta_{abaxep} = ', num2str(sol(9), '%4.1f'), ...
       '  \delta = ', num2str(sol(10), '%4.2f'), ...
       '  \mu = ', num2str(sol(11), '%4.2f'), ...
       '  v_{ad} = ', num2str(sol(12), '%4.1f'), ...
       '  v_{ab} = ', num2str(sol(13), '%4.2f')], ...
       'FontWeight', 'Bold');
axis([min(wl) max(wl) 0 1])
legend('Measured Adaxial BRF','Inverted Adaxial BRF ','Measured Abaxial BRF','Inverted Abaxial BRF')

