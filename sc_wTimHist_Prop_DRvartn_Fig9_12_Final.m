% This script file calculates the deflection of rail beam specifically to
% check the effect of change of zeta on the deflection values both for the
% Proposed and Winkler model(Commented, not needed).
% Both FIG9:(w_N vs x_N at two dr values) and 
% FIG12:(DAF vs zeta(dr) at different v_R) are 
% calculated through this code.
% Also saves the values at the end.


clear;clc;
% Input Parameters
sc_IP_Final1
% Axle Load
P_W=1.25e5;                        % Amplitude of load on rail beam(N)

% K_static Value
ha_k_st=@get_K_eqt_0_0n;
omega_0=0;
K_eqt_0=real(ha_k_st(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho_S,c,K1_b,K_b,K_P,c_P));

% Parameters defined for normalisation
v_cr_con=(4*E_R*I_R*K_eqt_0/(rho_R^2))^0.25;
c_cr_con=2*sqrt(K_eqt_0*rho_R);
lambda_1=(K_eqt_0/4/E_R/I_R)^0.25;

% Input velocity/damping parameters
% v_R_vect=linspace(0,2,1e3);
v_R_vect=[0.25 0.5 1 2]; %ALPHA
vel_vect=v_R_vect.*v_cr_con;

% Damping Ratio
% dr_vect=[0.05 0.3];
dr_vect=linspace(0,2,1e3);

% Input damping parameters
c_R_vect=dr_vect*c_cr_con;

% Input normalised distance/time
nd_vect=-50:0.1:50;                   % Normalised distance in Mallik's graphs
x_vect=nd_vect/lambda_1;
% t_vect=x_vect./vel_vect;

% Maximum static deflection (for normalisation)
ha_st=@get_static_max;
[w_st_max,w_vect_st]=ha_st(x_vect,K_eqt_0,P_W,E_R,I_R);

% Limits of integration
lim_1=-1e4;
lim_2=1e4;

% K(omega) for proposed model
load('K_eqt_Final1');
load('omega_vect_input');
%%%%%%%%%% NOTE: OMEGA_VECT_INPUT IS IN RADIANS/SECOND %%%%%%%%%%%%

% Exponential Window Method
o_I=0.0001;

% Function for calculating Deflection value for Proposed model
ha_In_t=@get_w_rail_omega;
% Function for calculating Deflection value for Winkler model
ha_In_t_1=@get_w_rail_omega_Kcnst;

% Deflection time histories and Deflection amplification Factor of Proposed and Winkler Model
In_t_vect_1=zeros(length(c_R_vect),length(x_vect));
In_t_vect_2=zeros(length(c_R_vect),length(x_vect));
w_tdr_vect_1=zeros(length(c_R_vect),length(x_vect));
W_DR_MAT=zeros(length(v_R_vect),length(c_R_vect));

% In_t_vect_Wink1=zeros(length(c_R_vect),length(x_vect));
% In_t_vect_Wink2=zeros(length(c_R_vect),length(x_vect));
% w_tdr_vect_Wink1=zeros(length(c_R_vect),length(x_vect));
% W_DR_MAT_Wink=zeros(length(v_R_vect),length(c_R_vect));

tic;
for k=1:length(v_R_vect)
    v_R=v_R_vect(1,k);
    vel=vel_vect(1,k);
    t_vect=x_vect./vel;
    if v_R>=0 && v_R<=0.000800160032006401
        w_tdr_vect_1=w_vect_st;
        w_tdr_vect_Wink1=w_vect_st;
    else
        for i=1:length(c_R_vect)
            c_R=c_R_vect(1,i);
            d_R=dr_vect(1,i);
            display(d_R)
            
            In_t_vect_1(i,:)=real(integral(@(omega_vect_1)ha_In_t(omega_vect_1,vel,x_R,t_vect,E_R,I_R,rho_R,c_R,o_I,omega_vect,K_eqt_vect),lim_1,lim_2,'ArrayValued',true));
            In_t_vect_2(i,:)=In_t_vect_1(i,:).*exp(o_I*t_vect);
            
%             In_t_vect_Wink1(i,:)=real(integral(@(omega_vect_1)ha_In_t_1(omega_vect_1,vel,x_R,t_vect,E_R,I_R,rho_R,c_R,o_I,K_eqt_0),lim_1,lim_2,'ArrayValued',true));
%             In_t_vect_Wink2(i,:)=In_t_vect_Wink1(i,:).*exp(o_I*t_vect);
            
%             w_tdr_vect_Wink1(i,:)=fliplr(-(P_W/pi/2).*In_t_vect_Wink2(i,:));
            w_tdr_vect_1(i,:)=fliplr(-(P_W/pi/2).*In_t_vect_2(i,:));    % Deflection when K is function of omega + No Damping
            
            if v_R>=0 && v_R<=0.000800160032006401
                W_DR_MAT(k,i)=w_st_max;
%                 W_DR_MAT_Wink(k,i)=w_st_max;
                
            else
                W_DR_MAT(k,i)=max(abs(w_tdr_vect_1(i,:)));
%                 W_DR_MAT_Wink(k,i)=max(abs(w_tdr_vect_Wink1(i,:)));
                
            end
        end
    end
    
end

t_exec_varun=toc;
% Normalised Values
w_dr_ratio=w_tdr_vect_1./w_st_max;
% w_dr_ratio_Wink=w_tdr_vect_Wink1./w_st_max;

W_ABS_DR_RMAT=W_DR_MAT./w_st_max;
% W_ABS_DR_RMAT_1=W_DR_MAT_Wink./w_st_max;

% Plotting Winkler
% figure
% plot(nd_vect,w_dr_ratio_Wink(1,:),'k')
% hold on
% plot(nd_vect,w_dr_ratio_Wink(2,:),'b')


% Plotting/Saving the values.
%% FIG 9
% figure
% plot(nd_vect,w_dr_ratio(1,:),'k')
% hold on
% plot(nd_vect,w_dr_ratio(2,:),'b')
% fileID=fopen('w_N_vs_x_N_Proposed_zeta=0.05_&_zeta=0.3_at_alpha=2.00.txt','w');
% file_1=[nd_vect;w_dr_ratio(1,:);w_dr_ratio(2,:)];
% fprintf(fileID,'%5s %17s %16s\n','x_N','w_N:zeta_0.05','w_N:zeta_0.3');
% fprintf(fileID,'%6.2f  %15.10f  %15.10f\n', file_1);

%% FIG 12
figure
plot(dr_vect,W_ABS_DR_RMAT(1,:),'k')
hold on
plot(dr_vect,W_ABS_DR_RMAT(2,:),'b')
hold on
plot(dr_vect,W_ABS_DR_RMAT(3,:),'r')
hold on
plot(dr_vect,W_ABS_DR_RMAT(4,:),'g')

fileID=fopen('DAF(w)_vs_zeta_@alpha=0.25_0.50_1.00_2.00.txt','w');
file_2=[dr_vect;W_ABS_DR_RMAT(1,:);W_ABS_DR_RMAT(2,:);W_ABS_DR_RMAT(3,:);W_ABS_DR_RMAT(4,:)];
fprintf(fileID,'%12s %15s %15s %15s %15s\n','dr(zeta)','alpha=0.25','alpha=0.50','alpha=1.0','alpha=2.0');

fprintf(fileID,'%12.8f %15.10f %15.10f %15.10f %15.10f\n', file_2);
