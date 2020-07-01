

function [var1, var2] = sensitivity_anlysis()


% import Battery Parameters
scott_battery_params

% Import Discrete State Space model
discrete_statespace_model


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


var1 = j1;
var2 = sen_out_spsi_p*epsilon_sp;




end 


% for k=1:Kup
%     %% electrolyte
% 
%   %% Negetive electrode three-state state space model for the particle
% 
%     I_Qdyn1(k)=I(k);
% 
%     Jn=I(k)/Vn;
%     ut_n=Jn;
%     i(k)=k*Ts;
%     yn(k)=c_n*xn(:,k);
%     xn(:,k+1)=a_n*xn(:,k)+b_n*ut_n;
% 
%    j0_n(k)=kn*F*((cen)*(cs_max_n-yn(k))*yn(k))^(0.5) ;
%     k_n(k)=Jn/(2*as_n*j0_n(k));
%     eta_n(k)=R*T*log(k_n(k)+(k_n(k)^2+1)^0.5)/(F*0.5);
% 
%     x(k)=yn(k)/(cs_max_n);
% 
%     %% Un(k)= interp1(x_un,y_un,x(k),'spline');
% 
%     %% Positive electrode three-state state space model for the particle
%     Jp=-I(k)/Vp; % current densitivity input for cathode is negative
%     ut_p=Jp;
% 
%     j1(k)=k*Ts;
% 
%     yp(k)=c_p*xp(:,k);
%     xp(:,k+1)=a_p*xp(:,k)+b_p*ut_p;
% 
%     j0_p(k)=F*kp*((cep)*(cs_max_p-yp(k))*yp(k))^(0.5);
%     k_p(k)=Jp/(2*as_p*j0_p(k));
%     eta_p(k)=R*T*log(k_p(k)+(k_p(k)^2+1)^0.5)/(F*0.5);
% 
%     y(k)=yp(k)/(cs_max_p);%yp is surface concentration
% 
%   %% state space realization for epsilon_sp
%  out_Sepsi_p(k)=Sepsi_c_p*Sepsi_p(:,k);
% Sepsi_p(:,k+1)=Sepsi_a_p*Sepsi_p(:,k)+Sepsi_b_p*I(k);  % current input for positive electrode is negative, ...
%                                                        % therefore the sensitivity output should be multiplied by -1;
%   %% state space realization for epsilon_sn
%  out_Sepsi_n(k)=Sepsi_c_n*Sepsi_n(:,k);
% Sepsi_n(:,k+1)=Sepsi_a_n*Sepsi_n(:,k)+Sepsi_b_n*I(k);
% 
% %% state space realization for D_sp
%  out_Sdsp_p(k)=Sdsp_c_p*Sdsp_p(:,k);
% Sdsp_p(:,k+1)=Sdsp_a_p*Sdsp_p(:,k)+Sdsp_b_p*I(k);
% 
% %% state space realization for D_sn
%  out_Sdsn_n(k)=Sdsn_c_n*Sdsn_n(:,k);
% Sdsn_n(:,k+1)=Sdsn_a_n*Sdsn_n(:,k)+Sdsn_b_n*I(k);
% %% OcP slope
% docvp_dCsep(k) = 0.07645*(-54.4806/cs_max_p)* ...
% ((1.0./cosh(30.834-54.4806*y(k))).^2) ...
% +2.1581*(-50.294/cs_max_p)*((cosh(52.294-50.294*y(k))).^(-2)) ...
% +0.14169*(19.854/cs_max_p)*((cosh(11.0923-19.8543*y(k))).^(-2)) ...
% -0.2051*(5.4888/cs_max_p)*((cosh(1.4684-5.4888*y(k))).^(-2)) ...
% -0.2531/0.1316/cs_max_p*((cosh((-y(k)+0.56478)/0.1316)).^(-2)) ...
%  - 0.02167/0.006/cs_max_p*((cosh((y(k)-0.525)/0.006)).^(-2));
% 
% docvn_dCsen(k)= -1.5*(120.0/cs_max_n)*exp(-120.0*x(k))  ...
%  +(0.0351/(0.083*cs_max_n))*((cosh((x(k)-0.286)/0.083)).^(-2)) ...
%  -(0.0045/(cs_max_n*0.119))*((cosh((x(k)-0.849)/0.119)).^(-2)) ...
%  -(0.035/(cs_max_n*0.05))*((cosh((x(k)-0.9233)/0.05)).^(-2)) ...
%  -(0.0147/(cs_max_n*0.034))*((cosh((x(k)-0.5)/0.034)).^(-2)) ...
%  -(0.102/(cs_max_n*0.142))*((cosh((x(k)-0.194)/0.142)).^(-2)) ...
%  -(0.022/(cs_max_n*0.0164))*((cosh((x(k)-0.9)/0.0164)).^(-2)) ...
%  -(0.011/(cs_max_n*0.0226))*((cosh((x(k)-0.124)/0.0226)).^(-2)) ...
%  +(0.0155/(cs_max_n*0.029))*((cosh((x(k)-0.105)/0.029)).^(-2));
% 
% %%
% rho1p_1(k)=-sign(I(k))*(-3*R*T)/(0.5*F*Rp*as_p) *((1+1/(k_p(k))^2)^(-0.5));
% rho1p(k)=R*T/(0.5*F) * (1/(k_p(k)+(k_p(k)^2+1)^0.5)) * (1+k_p(k)/((k_p(k)^2+1)^0.5)) *(-3*Jp/(2*as_p^2*j0_p(k)*Rp));
% rho2p(k)=(R*T)/(2*0.5*F)  * (cep*cs_max_p-2*cep*yp(k))/(cep*yp(k)*(cs_max_p-yp(k))) * (1+1/(k_p(k))^2)^(-0.5);
% 
% 
% rho1n_1(k)=sign(I(k))*(-3*R*T)/(0.5*F*Rn*as_n) *((1+1/(k_n(k))^2)^(-0.5));
% rho1n(k)=R*T/(0.5*F) * (1/(k_n(k)+(k_n(k)^2+1)^0.5)) * (1+k_n(k)/((k_n(k)^2+1)^0.5)) *(-3*Jn/(2*as_n^2*j0_n(k)*Rn));
% 
% rho2n(k)=(-R*T)/(2*0.5*F) *(cen*cs_max_n-2*cen*yn(k))/(cen*yn(k)*(cs_max_n-yn(k))) *(1+1/(k_n(k))^2)^(-0.5);
% 
% %% sensitivity of epsilon_sp    epsilon_sn
%  sen_out_spsi_p(k)=((rho1p(k)+(rho2p(k)+docvp_dCsep(k))*-out_Sepsi_p(k)));
%  sen_out_spsi_n(k)=((rho1n(k)+(rho2n(k)+docvn_dCsen(k))*out_Sepsi_n(k)));
% 
% out_deta_p_desp(k)=rho1p(k)+(rho2p(k))*-out_Sepsi_p(k);
% out_deta_n_desn(k)=rho1n(k)+(rho2n(k))*out_Sepsi_n(k);
% 
% out_semi_linear_p(k)=(docvp_dCsep(k))*out_Sepsi_p(k);
% out_semi_linear_n(k)=(docvn_dCsen(k))*out_Sepsi_n(k);
% 
% %% sensiviity of Dsp Dsn
% sen_out_ds_p(k)=(((rho2p(k)+docvp_dCsep(k))*-out_Sdsp_p(k)))*Ds_p;
% sen_out_ds_n(k)=(((rho2n(k)+docvn_dCsen(k))*out_Sdsn_n(k)))*Ds_n;
% 
% 
%  dV_dDsp(k)=sen_out_ds_p(k);
%  dV_dDsn(k)=sen_out_ds_n(k);
% 
%  dCse_dDsp(k)=-out_Sdsp_p(k)*Ds_p;
%  dCse_dDsn(k)=out_Sdsn_n(k)*Ds_n;
% 
% k=k+1;
% end


% figure(1)
% plot(j1,sen_out_spsi_p*epsilon_sp,'--','LineWidth',1.4)
% hold on