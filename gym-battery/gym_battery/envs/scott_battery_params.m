% current sign is positive for discharge/ negetive for charge


%% Scott's parameters
epsilon_sn=0.6;   % average
epsilon_sp=0.50;   % average positive active volume fraction
epsilon_e_n=0.3;  % Liquid volume fraction
epsilon_e_p=0.3;
F=96485.3329;   % Faraday constant
Rn=10e-6 ;   % Negative active particle radius
Rp=10e-6;
R=8.314 ;  % Universial gas constant
T=298.15;
Ar_n=1;      % Current collector area
Ar_p=1;
Ln=100e-6  ;   % Negative electrode thickness
Lp=100e-6;
Lsep=25e-6;
Lc=Ln+Lp+Lsep;
Ds_n=3.9e-14 ;    % Negative solid phase diffusion coefficient
Ds_p=(1e-13)*0.5;  %Positive solid phase diffusion coefficient
kn=1e-5/F;  %Rate constant of exchange current density
kp=3e-7/F;
De=2.7877e-10;

 stoi_x0=0.0069;
 stoi_x100=0.6760;

 stoi_y0= 0.8228;
 stoi_y100=0.442;

SOC=0; % SOC can change from 0 to 1

stoi_x=(stoi_x100-stoi_x0)*SOC+stoi_x0;
stoi_y=stoi_y0-(stoi_y0-stoi_y100)*SOC;

cs_max_n=3.6e3 * 372 * 1800 /F; % 0.035~0.870=1.0690e+03~ 2.6572e+04
cs_max_p=3.6e3 * 274 * 5010 / F; % Positive electrode  maximum solid-phase concntration 0.872~0.278=  4.3182e+04~1.3767e+04

Rf=1e-3;
as_n=3*epsilon_sn/Rn ;
as_p=3*epsilon_sp/Rp;   %  Active surface area per electrode unit volume
Vn=Ar_n*Ln  ;      % electrode volume
Vp=Ar_p*Lp  ;
t_plus=0.4;
cep=1000;
cen=1000;

% % load pulse profile
Kup=1300; %1C 3601; 2C 1712; 3C 1083;
I(1)=0;
 for k=2:Kup

    %  I(k)=1;
   I(k)=-25.67*3;
%     % I(k)=25.5;
 end

Kup=length(I);

