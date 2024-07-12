% _______________________________________________________________________
% sMADLM.m (July 1, 2024)
% a simplified version of MADLM
% This script allows to retrieve biochemical parameters for one dorsiventral leaf at individual vewing angle 
% or for a single pixel based on the DLM model (Stuckens et al., 2009) and a simplified specular reflection term Rspec.
% _______________________________________________________________________
% ***********************************************************************
% Dongjie Ran, Zhongqiu Sun, Shan Lu, Kenji Omasa, 2024 
% An advanced dorsiventral leaf radiative transfer model for simulating multi-angular
% and spectral reflection: considering asymmetry of leaf internal and surface structure. Remote Sensing of Environment.
% ***********************************************************************


function BRF_output=sMADLM(a,wl,data)
% a: Biochemical and structural parameters of MADLM model:
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
%
% wl: Wavelength (nm)
%% Diffuse component calculation
Cab=a(1);
Cxc=a(2);
Cw=a(3);
Cm=a(4);
Cbrown=a(5);
fair=a(6);
beta_pigm=a(7);
beta_wdm=a(8);
beta_abaxep=a(9);
abdiff=a(10);
mu=a(11);
vad=a(12);
vab=a(13);

absparams=[Cab,Cxc,Cw,Cm,Cbrown];
strparams=[fair,beta_pigm,beta_wdm,abdiff,beta_abaxep,mu];

RT=DLMdiff(absparams,strparams,data);% [wavelength Rdiff_adaxial transmittance_adaxial Rdiff_abaxial transmittance_abaxial]
Rdiff_ad=RT(:,2);
Rdiff_ab=RT(:,4);
%% Specular component calculation 
n_ep=data(:,3); %epidermal refractive index
F=((n_ep-1)./(n_ep+1)).^2; % Fresnel factor 
Rspec_ad=vad.*F;
Rspec_ab=vab.*F;
%% Total BRF calculation
BRFad=Rdiff_ad+Rspec_ad;
BRFab=Rdiff_ab+Rspec_ab;
BRF_output=[BRFad BRFab];
end