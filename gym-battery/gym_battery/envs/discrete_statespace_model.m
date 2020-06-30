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


% plot(j1,sen_out_spsi_p*epsilon_sp,'--','LineWidth',1.4)
% hold on
