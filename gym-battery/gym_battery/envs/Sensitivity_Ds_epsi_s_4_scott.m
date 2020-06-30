clear all
close all
clc

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

%% Negtive electrode three-state state space model for the particle

rfa_n=1/(F*as_n*Rn^5);
rfa_p=1/(F*as_p*Rp^5);
% matrix
An=[0 1 0; 0 0 1;0 -(3465*(Ds_n^2)/Rn^4) -(189*Ds_n/Rn^2)];
Bn=[0;0;-1];
Cn=rfa_n*[10395*Ds_n^2  1260*Ds_n*Rn^2   21*Rn^4];
Dn=[0];

% Approximate Discretization
Ts=1;
[n,m]=size(An);
A_dn=eye(n)+An*Ts;
B_dn=Bn*Ts;

a_n=A_dn;
b_n=B_dn;
c_n=Cn;
d_n=Dn;

%% discharge
xn(:,1)=[stoi_x *cs_max_n/(rfa_n*10395*Ds_n^2);0;0];  % stoi_x100 should be changed if the initial soc is not equal to 50%
xp(:,1)=[stoi_y *cs_max_p/(rfa_p*10395*Ds_p^2);0;0];  % initial positive electrode ion concentration

%% Positive electrode three-state state space model for the particle

% matrix
Ap=1*[0 1 0; 0 0 1;0 -(3465*(Ds_p^2)/Rp^4) -(189*Ds_p/Rp^2)];
Bp=[0;0;-1];
Cp=rfa_p*[10395*Ds_p^2  1260*Ds_p*Rp^2   21*Rp^4];
Dp=[0];

% Approximate Discretization

[n,m]=size(Ap);
A_dp=eye(n)+Ap*Ts;
B_dp=Bp*Ts;
a_p=A_dp;
b_p=B_dp;
c_p=Cp;
d_p=Dp;

%% sensitivity realization in time domain for epsilon_sp from third order pade(you can refer to my slides)
coef= (3)/(F*Rp^6 *as_p^2 *Ar_p*Lp);

Sepsi_A=[0 1 0 ;0 0 1; 0  -(3465*Ds_p^2)/Rp^4  -(189*Ds_p)/Rp^2 ];
Sepsi_B=[0;0;1];
Sepsi_C=coef*[10395*Ds_p^2   1260*Ds_p*Rp^2  21*Rp^4];
Sepsi_D=[0];

[n,m]=size(Sepsi_A);
Sepsi_A_dp=eye(n)+Sepsi_A*Ts;
Sepsi_B_dp=Sepsi_B*Ts;
Sepsi_a_p=Sepsi_A_dp;
Sepsi_b_p=Sepsi_B_dp;
Sepsi_c_p=Sepsi_C;
Sepsi_d_p=Sepsi_D;

Sepsi_p(:,1)=[0;0;0];
%% sensitivity realization in time domain for epsilon_sn from third order pade(you can refer to my slides)
coefn= (3)/(F*Rn^6 *as_n^2 *Ar_n*Ln);

Sepsi_A_n=[0 1 0 ;0 0 1; 0  -(3465*Ds_n^2)/Rn^4  -(189*Ds_n)/Rn^2 ];
Sepsi_B_n=[0;0;1];
Sepsi_C_n=coefn*[10395*Ds_n^2   1260*Ds_n*Rn^2  21*Rn^4];
Sepsi_D_n=[0];

[n,m]=size(Sepsi_A_n);
Sepsi_A_dn=eye(n)+Sepsi_A_n*Ts;
Sepsi_B_dn=Sepsi_B_n*Ts;
Sepsi_a_n=Sepsi_A_dn;
Sepsi_b_n=Sepsi_B_dn;
Sepsi_c_n=Sepsi_C_n;
Sepsi_d_n=Sepsi_D_n;

Sepsi_n(:,1)=[0;0;0];

%% sensitivity realization in time domain for D_sp from third order pade

coefDsp= (63*Rp)/(F*as_p*Ar_p*Lp*Rp^8);

Sdsp_A=[0 1 0 0;0 0 1 0; 0  0 0 1; -(12006225*Ds_p^4)/Rp^8  -1309770*Ds_p^3/Rp^6   -42651*Ds_p^2/Rp^4   -378*Ds_p/Rp^2];
Sdsp_B=[0;0;0;1];
Sdsp_C=coefDsp*[38115*Ds_p^2   1980*Ds_p*Rp^2   43*Rp^4  0];
Sdsp_D=[0];

[n,m]=size(Sdsp_A);
Sdsp_A_dp=eye(n)+Sdsp_A*Ts;
Sdsp_B_dp=Sdsp_B*Ts;
Sdsp_a_p=Sdsp_A_dp;
Sdsp_b_p=Sdsp_B_dp;
Sdsp_c_p=Sdsp_C;
Sdsp_d_p=Sdsp_D;

Sdsp_p(:,1)=[0;0;0; 0];
%% sensitivity realization in time domain for D_sn from third order pade

coefDsn= (63*Rn)/(F*as_n*Ar_n*Ln*Rn^8);

Sdsn_A=[0 1 0 0;0 0 1 0; 0  0 0 1; -(12006225*Ds_n^4)/Rn^8  -1309770*Ds_n^3/Rn^6   -42651*Ds_n^2/Rn^4   -378*Ds_n/Rn^2];
Sdsn_B=[0;0;0;1];
Sdsn_C=coefDsn*[38115*Ds_n^2   1980*Ds_n*Rn^2   43*Rn^4  0];
Sdsn_D=[0];

[n,m]=size(Sdsn_A);
Sdsn_A_dn=eye(n)+Sdsn_A*Ts;
Sdsn_B_dn=Sdsn_B*Ts;
Sdsn_a_n=Sdsn_A_dn;
Sdsn_b_n=Sdsn_B_dn;
Sdsn_c_n=Sdsn_C;
Sdsn_d_n=Sdsn_D;

Sdsn_n(:,1)=[0;0;0; 0];

%%

for k=1:Kup
    %% electrolyte

  %% Negetive electrode three-state state space model for the particle
  
    I_Qdyn1(k)=I(k);
 
    Jn=I(k)/Vn; 
    ut_n=Jn;
    i(k)=k*Ts;
    yn(k)=c_n*xn(:,k);
    xn(:,k+1)=a_n*xn(:,k)+b_n*ut_n;

   j0_n(k)=kn*F*((cen)*(cs_max_n-yn(k))*yn(k))^(0.5) ;
    k_n(k)=Jn/(2*as_n*j0_n(k));
    eta_n(k)=R*T*log(k_n(k)+(k_n(k)^2+1)^0.5)/(F*0.5);
   
    x(k)=yn(k)/(cs_max_n); 

    %% Un(k)= interp1(x_un,y_un,x(k),'spline'); 
    
    %% Positive electrode three-state state space model for the particle
    Jp=-I(k)/Vp; % current densitivity input for cathode is negative
    ut_p=Jp; 
    
    j1(k)=k*Ts;
 
    yp(k)=c_p*xp(:,k);
    xp(:,k+1)=a_p*xp(:,k)+b_p*ut_p;

    j0_p(k)=F*kp*((cep)*(cs_max_p-yp(k))*yp(k))^(0.5); 
    k_p(k)=Jp/(2*as_p*j0_p(k));
    eta_p(k)=R*T*log(k_p(k)+(k_p(k)^2+1)^0.5)/(F*0.5);
   
    y(k)=yp(k)/(cs_max_p);%yp is surface concentration
    
  %% state space realization for epsilon_sp
 out_Sepsi_p(k)=Sepsi_c_p*Sepsi_p(:,k);
Sepsi_p(:,k+1)=Sepsi_a_p*Sepsi_p(:,k)+Sepsi_b_p*I(k);  % current input for positive electrode is negative, ...
                                                       % therefore the sensitivity output should be multiplied by -1; 
  %% state space realization for epsilon_sn
 out_Sepsi_n(k)=Sepsi_c_n*Sepsi_n(:,k);
Sepsi_n(:,k+1)=Sepsi_a_n*Sepsi_n(:,k)+Sepsi_b_n*I(k);

%% state space realization for D_sp
 out_Sdsp_p(k)=Sdsp_c_p*Sdsp_p(:,k);
Sdsp_p(:,k+1)=Sdsp_a_p*Sdsp_p(:,k)+Sdsp_b_p*I(k);

%% state space realization for D_sn
 out_Sdsn_n(k)=Sdsn_c_n*Sdsn_n(:,k);
Sdsn_n(:,k+1)=Sdsn_a_n*Sdsn_n(:,k)+Sdsn_b_n*I(k);
%% OcP slope
docvp_dCsep(k) = 0.07645*(-54.4806/cs_max_p)* ...
((1.0./cosh(30.834-54.4806*y(k))).^2) ...
+2.1581*(-50.294/cs_max_p)*((cosh(52.294-50.294*y(k))).^(-2)) ...
+0.14169*(19.854/cs_max_p)*((cosh(11.0923-19.8543*y(k))).^(-2)) ...
-0.2051*(5.4888/cs_max_p)*((cosh(1.4684-5.4888*y(k))).^(-2)) ...
-0.2531/0.1316/cs_max_p*((cosh((-y(k)+0.56478)/0.1316)).^(-2)) ...
 - 0.02167/0.006/cs_max_p*((cosh((y(k)-0.525)/0.006)).^(-2));

docvn_dCsen(k)= -1.5*(120.0/cs_max_n)*exp(-120.0*x(k))  ...
 +(0.0351/(0.083*cs_max_n))*((cosh((x(k)-0.286)/0.083)).^(-2)) ...
 -(0.0045/(cs_max_n*0.119))*((cosh((x(k)-0.849)/0.119)).^(-2)) ...
 -(0.035/(cs_max_n*0.05))*((cosh((x(k)-0.9233)/0.05)).^(-2)) ...
 -(0.0147/(cs_max_n*0.034))*((cosh((x(k)-0.5)/0.034)).^(-2)) ...
 -(0.102/(cs_max_n*0.142))*((cosh((x(k)-0.194)/0.142)).^(-2)) ...
 -(0.022/(cs_max_n*0.0164))*((cosh((x(k)-0.9)/0.0164)).^(-2)) ...
 -(0.011/(cs_max_n*0.0226))*((cosh((x(k)-0.124)/0.0226)).^(-2)) ...
 +(0.0155/(cs_max_n*0.029))*((cosh((x(k)-0.105)/0.029)).^(-2)); 

%%
rho1p_1(k)=-sign(I(k))*(-3*R*T)/(0.5*F*Rp*as_p) *((1+1/(k_p(k))^2)^(-0.5));
rho1p(k)=R*T/(0.5*F) * (1/(k_p(k)+(k_p(k)^2+1)^0.5)) * (1+k_p(k)/((k_p(k)^2+1)^0.5)) *(-3*Jp/(2*as_p^2*j0_p(k)*Rp));
rho2p(k)=(R*T)/(2*0.5*F)  * (cep*cs_max_p-2*cep*yp(k))/(cep*yp(k)*(cs_max_p-yp(k))) * (1+1/(k_p(k))^2)^(-0.5);


rho1n_1(k)=sign(I(k))*(-3*R*T)/(0.5*F*Rn*as_n) *((1+1/(k_n(k))^2)^(-0.5));
rho1n(k)=R*T/(0.5*F) * (1/(k_n(k)+(k_n(k)^2+1)^0.5)) * (1+k_n(k)/((k_n(k)^2+1)^0.5)) *(-3*Jn/(2*as_n^2*j0_n(k)*Rn));

rho2n(k)=(-R*T)/(2*0.5*F) *(cen*cs_max_n-2*cen*yn(k))/(cen*yn(k)*(cs_max_n-yn(k))) *(1+1/(k_n(k))^2)^(-0.5);

%% sensitivity of epsilon_sp    epsilon_sn
 sen_out_spsi_p(k)=((rho1p(k)+(rho2p(k)+docvp_dCsep(k))*-out_Sepsi_p(k)));
 sen_out_spsi_n(k)=((rho1n(k)+(rho2n(k)+docvn_dCsen(k))*out_Sepsi_n(k)));

out_deta_p_desp(k)=rho1p(k)+(rho2p(k))*-out_Sepsi_p(k);
out_deta_n_desn(k)=rho1n(k)+(rho2n(k))*out_Sepsi_n(k);

out_semi_linear_p(k)=(docvp_dCsep(k))*out_Sepsi_p(k);
out_semi_linear_n(k)=(docvn_dCsen(k))*out_Sepsi_n(k);

%% sensiviity of Dsp Dsn
sen_out_ds_p(k)=(((rho2p(k)+docvp_dCsep(k))*-out_Sdsp_p(k)))*Ds_p;
sen_out_ds_n(k)=(((rho2n(k)+docvn_dCsen(k))*out_Sdsn_n(k)))*Ds_n;

 
 dV_dDsp(k)=sen_out_ds_p(k);
 dV_dDsn(k)=sen_out_ds_n(k);
 
 dCse_dDsp(k)=-out_Sdsp_p(k)*Ds_p;
 dCse_dDsn(k)=out_Sdsn_n(k)*Ds_n;
 
k=k+1;
end



figure(1)
plot(j1,sen_out_spsi_p*epsilon_sp,'--','LineWidth',1.4)
hold on
