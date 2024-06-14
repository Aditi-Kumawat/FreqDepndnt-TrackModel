% PLOTS AND SAVES THE Fig 7 DATA

clear;clc;
% Isolation Frequency
f0=29.456652092120930;
omega_0_1=2*pi*f0;

% Input frequency
omega_vect_req=linspace(0,4000,1e4);
% omega_vect_req=linspace(0,400,1e3);

% Saved K(Omega) and corresponding omega_Vect
load('K_eqt_Final1');
load('omega_vect_input');

% K(Omega) at the required omega
ha_k=@get_Komega_vect_intrpol;
K_eqt_omega_vect=ha_k(omega_vect,K_eqt_vect,omega_vect_req);

% K_static
sc_IP_Final1
ha_k_st=@get_K_eqt_0_0n;
omega_0=0;
K_eqt_0=real(ha_k_st(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho_S,c,K1_b,K_b,K_P,c_P));

% Normalized values
K_eqt_0_vect=K_eqt_0*ones(1,length(omega_vect_req));
K_ratio_vect=real(K_eqt_omega_vect./K_eqt_0_vect);
zero_vect=zeros(1,length(omega_vect_req));
k_static_vect=abs(K_eqt_omega_vect./K_eqt_omega_vect);

% Plotting
figure;
ha_K_R=tight_subplot(1,1,[.05 .13],[0.11 0.1],[0.12 0.08]);
hold on;
axes(ha_K_R(1));
plot(omega_vect_req,k_static_vect,':k','LineWidth',2);
hold on
plot(omega_vect_req,real(K_ratio_vect),'-k','LineWidth',1.5);
hold on
plot(omega_vect_req,zero_vect,'--k','LineWidth',1);
grid on;
% set(gcf, 'Units','inches','Position', [4 4 3.5 2.5]);
set(ha_K_R,'XTickLabelMode','auto');
set(ha_K_R,'FontSize',12);
set(ha_K_R,'YTickLabelMode','auto');
set(ha_K_R,'FontSize',12);

% Saving as text file
% fileID=fopen('k_eqt(omega)_by_k_static_vs_omega_range_0_4000(radpersec).txt','w');
% file_2=[omega_vect_req;real(K_ratio_vect);k_static_vect;zero_vect];
% fprintf(fileID,'%15s %15s %15s %15s\n','omega(rad/s)','K_eqt/Kst','K_static','zero_vect');

% fprintf(fileID,'%15.10f %15.10f %15.10f %15.10f\n', file_2);

