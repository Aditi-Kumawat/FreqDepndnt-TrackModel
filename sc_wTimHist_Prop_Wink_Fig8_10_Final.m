% This script file compares the deflection of rail beam for
% Proposed and Winkler models.
% FIG8:(w_N vs x_N at different v_R(alpha=0.25,0.5,1,2) values and zeta=0.05)  
% FIG10:(DAF vs alpha(v_R) at zeta=0.05) are 
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
w0_IF_rail=sqrt(K_eqt_0/rho_R);
f0_IF_rail=w0_IF_rail/2/pi;
dif=((w0_IF_rail-w0_IF_sleeper)*100)/w0_IF_sleeper;
% Input velocity/damping parameters
v_R_vect=linspace(0,2,5e2);
% v_R_vect=0.5; %ALPHA
vel_vect=v_R_vect.*v_cr_con;
% dr_vect=[0.05 0.3 1.1];

% Damping Ratio
dr_vect=0.05;

% Input damping parameters
c_R_vect=dr_vect*c_cr_con;

% Input normalised distance/time
nd_vect=-6:0.1:6;                   % Normalised distance in Mallik's graphs
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
In_t_vect_1=zeros(length(v_R_vect),length(x_vect));
In_t_vect_2=zeros(length(v_R_vect),length(x_vect));
w_t_vect=zeros(length(v_R_vect),length(x_vect));
W_ABS_MAT=zeros(length(c_R_vect),length(v_R_vect));

In_t_vect_Wink1=zeros(length(v_R_vect),length(x_vect));
In_t_vect_Wink2=zeros(length(v_R_vect),length(x_vect));
w_t_vect_Wink=zeros(length(v_R_vect),length(x_vect));
W_ABS_MAT_Wink=zeros(length(c_R_vect),length(v_R_vect));

tic;
for k=1:length(c_R_vect)
    c_R=c_R_vect(1,k);
    
    for i=1:length(v_R_vect)
        v_R=v_R_vect(1,i);
        display(v_R)
        
        if v_R>=0 && v_R<=0.000800160032006401
            w_t_vect=w_vect_st;
            w_t_vect_Wink=w_vect_st;
        else
            vel=vel_vect(1,i);
            t_vect=x_vect./vel;
            
            In_t_vect_1(i,:)=real(integral(@(omega_vect_1)ha_In_t(omega_vect_1,vel,x_R,t_vect,E_R,I_R,rho_R,c_R,o_I,omega_vect,K_eqt_vect),lim_1,lim_2,'ArrayValued',true));
            In_t_vect_2(i,:)=In_t_vect_1(i,:).*exp(o_I*t_vect);
            
            In_t_vect_Wink1(i,:)=real(integral(@(omega_vect_1)ha_In_t_1(omega_vect_1,vel,x_R,t_vect,E_R,I_R,rho_R,c_R,o_I,K_eqt_0),lim_1,lim_2,'ArrayValued',true));
            In_t_vect_Wink2(i,:)=In_t_vect_Wink1(i,:).*exp(o_I*t_vect);
            
            w_t_vect(i,:)=fliplr(-(P_W/pi/2).*In_t_vect_2(i,:));    % Deflection when K is function of omega + No Damping
            w_t_vect_Wink(i,:)=fliplr(-(P_W/pi/2).*In_t_vect_Wink2(i,:));    % Deflection when K is function of omega + No Damping
        end
        
        if v_R>=0 && v_R<=0.000800160032006401
            W_ABS_MAT(k,i)=w_st_max;
            W_ABS_MAT_Wink(k,i)=w_st_max;
        else
            W_ABS_MAT(k,i)=max(abs(w_t_vect(i,:)));
            W_ABS_MAT_Wink(k,i)=max(abs(w_t_vect_Wink(i,:)));
        end
    end
end

t_exec=toc;

% Normalised values
w_ratio=w_t_vect./w_st_max;
W_ABS_RMAT=W_ABS_MAT./w_st_max;

w_ratio_Wink=w_t_vect_Wink./w_st_max;
W_ABS_RMAT_Wink=W_ABS_MAT_Wink./w_st_max;

idx_1=(W_ABS_RMAT==max(W_ABS_RMAT));
idx_2=(W_ABS_RMAT_Wink==max(W_ABS_RMAT_Wink));
vr1=v_R_vect(idx_1);
vr2=v_R_vect(idx_2);
dif2=(vr2-vr1)/vr1;
% Plotting and Saving
%% FIG 8
% figure
% plot(nd_vect,w_ratio,'k')
% hold on
% plot(nd_vect,w_ratio_Wink,'b')

% fileID=fopen('w_N_vs_x_N_Proposed_Winkler_zeta=0.05_alpha=2.00.txt','w');
% % file_1=[nd_vect;w_ratio;w_ratio_Wink];
% fprintf(fileID,'%5s %17s %16s\n','x_N','w_N:Proposed','w_N:Winkler');
% fprintf(fileID,'%6.2f  %15.10f  %15.10f\n', file_1);

% %% FIG 10
figure
plot(v_R_vect,W_ABS_RMAT)
hold on
plot(v_R_vect,W_ABS_RMAT_Wink)
% fileID=fopen('DAF(w)_vs_alpha(v_R)_Proposed_Winkler_@zeta=0.05.txt','w');
% 
% file_2=[v_R_vect;W_ABS_RMAT;W_ABS_RMAT_Wink];
% fprintf(fileID,'%12s %19s %18s\n','v_R(alpha)','Proposed Model','Winkler Model');
% 
% fprintf(fileID,'%12.8f %18.10f %18.10f\n', file_2);
