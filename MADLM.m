% _______________________________________________________________________
% MADLM.m (July 1, 2024)
%
% This script used to forward simulate BRF of the dorsiventral leaf.
% _______________________________________________________________________
% ***********************************************************************
% Dongjie Ran, Zhongqiu Sun, Shan Lu, Kenji Omasa, 2024 
% An advanced dorsiventral leaf radiative transfer model for simulating multi-angular
% and spectral reflection: considering asymmetry of leaf internal and surface structure. Remote Sensing of Environment.
% ***********************************************************************
% This script allows to simulate multi-angular BRFs of the adaxial and abaxial sides based on the DLM model (Stuckens et al., 2009)
% and a specular reflection function (Bousquet et al. 2005)
% Bousquet, L., Lach¨¦rade, S., Jacquemoud, S., & Moya, I., 2005. Leaf BRDF measurements and model
% for specular and diffuse components differentiation. Remote Sensing of Environment, 98, 201-211.
% ***********************************************************************
% Leaf BRF are calculated from 400 nm to 2500 nm (1 nm step)
% _______________________________________________________________________

function BRF_output=MADLM(a,wl,geomparams)
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
%    - uad: Roughness of adaxial leaf surface
%    - uab: Roughness of abaxial leaf surface
%    - fad: Fraction of refractive index for the adaxial leaf surface
%    - fab: Fraction of refractive index for the abaxial leaf surface

% geomparams: Source zenith angle, viewing zenith angle, and viewing azimuth angle in radians
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
uad=a(12);
uab=a(13);
fad=a(14);
fab=a(15);

absparams=[Cab,Cxc,Cw,Cm,Cbrown];
strparams=[fair,beta_pigm,beta_wdm,abdiff,beta_abaxep,mu];
datas=dataspecDLM;
RT=DLMdiff(absparams,strparams,datas);% [wavelength Rdiff_adaxial transmittance_adaxial Rdiff_abaxial transmittance_abaxial]
Rdiff_ad=RT(:,2);
Rdiff_ab=RT(:,4);
%% Specular component calculation
 
BRFspecad=Cal_spec(geomparams,fad,uad);
BRFspecab=Cal_spec(geomparams,fab,uab);

function rspec=Cal_spec(geomparams,f,sigma)
data=dataspecDLM;
n_ep=data(:,3);
fn=n_ep'.*f; %refractive index

% illumination-viewing geometry
thetas=geomparams(:,1); % source zenith angle, radian
thetav=geomparams(:,2); % viewing zenith angle, radian
fav=geomparams(:,3); % viewing azimuth angle, radian
thetaa=0.5.*acos(cos(thetas).*cos(thetav)+sin(thetas).*sin(thetav).*cos(fav));
cosalfa=(cos(thetas)+cos(thetav))./(2.*cos(thetaa));

E1=2.*cosalfa.*cos(thetav)./cos(thetaa);
E2=2.*cosalfa.*cos(thetas)./cos(thetaa);
G=(min(1,min(E1,E2)));

D=exp(-(tan(acos(cosalfa))).^2./sigma.^2)./(sigma.^2.*cosalfa.^4);

g=sqrt(abs(fn.^2+(cos(thetaa)).^2-1));
F=0.5.*((g-cos(thetaa))./(g+cos(thetaa))).^2.*(1+((cos(thetaa).*(g+cos(thetaa))-1)./(cos(thetaa).*(g-cos(thetaa))+1)).^2);

BRDFspec=F.*D.*G./(4.*pi.*cos(thetas).*cos(thetav));
BRFs=pi.*BRDFspec';
rspec=BRFs;
end
%% Total BRF calculation
BRFad=Rdiff_ad+BRFspecad;
BRFab=Rdiff_ab+BRFspecab;
BRF_output=[BRFad(wl,:) BRFab(wl,:)];
end