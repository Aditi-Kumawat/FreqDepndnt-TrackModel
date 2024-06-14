% This script file calculates the load on sleeper at a given location on
% the rail beam for the different values of parameters assigned in
% sc_P_t_vect (x_R,vel)
% Using the interpolation function get_P_t, this code calculates the load
% on sleeper for the same parameters (as defined in sc_P_t_vect) but for the
% multiple wheel loads(N_w).
% Uses the results generated from sc_P_t_vect.m i.e."P_t_vect" and
% "t_vect".

clear;clc;
P_1_1=200.4e3;
P_1_2=117.7e3;



P_W_vect_2=[P_1_1 P_1_2];                        % Amplitude of load on rail beam(N)
P_W_vect_1=(P_W_vect_2);
% P_W_vect=[1.81e5 1.8e5 1.22e5 1.22e5 1.22e5 1.22e5 1.22e5 1.22e5 1.17e5 1.6e5];                        % Amplitude of load on rail beam(N)
N_W=24;
v=201; %velocity in kmph
vel=v*1000/60/60;

delta_vect_1=[0 3 7.5 3 6.35 2.56 15.44 2.56 5.8 2.56 15.44 2.56 5.8 2.56 15.44 2.56 5.8 2.56 15.44 2.56 5.8 2.56 15.44 2.56];  %distance b/W wheels

delta_vect=zeros(1,length(delta_vect_1));

for d=1:length(delta_vect)
    for e=1:d
        delta_vect(1,d)=delta_vect(1,d)+delta_vect_1(1,(d+1-e));
    end
end
sc_IP_multi
ha_k_st=@get_K_eqt_0_val2;
omega_0=0;
K_eqt_0=real(ha_k_st(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho_S,c,K1_b,K_b,K_P,c_P));

ha_P=@get_P_t;
t_vect=-5:0.01:5;
ha_w=@get_wrailKcnst_forP1;
w_mat_1=zeros(length(P_W_vect_1),length(t_vect));


for k=1:length(P_W_vect_1)
    P_W=P_W_vect_1(k);
    display(P_W)
    w_mat_1(k,:)=ha_w(P_W,vel,t_vect,rho_R,x_R,E_R,I_R,dr,K_eqt_0);
end
w_mat=[w_mat_1(1,:);w_mat_1(1,:);w_mat_1(1,:);w_mat_1(1,:);...
    w_mat_1(2,:);w_mat_1(2,:);w_mat_1(2,:);w_mat_1(2,:);...
    w_mat_1(2,:);w_mat_1(2,:);w_mat_1(2,:);w_mat_1(2,:);...
    w_mat_1(2,:);w_mat_1(2,:);w_mat_1(2,:);w_mat_1(2,:);...
    w_mat_1(2,:);w_mat_1(2,:);w_mat_1(2,:);w_mat_1(2,:);...
    w_mat_1(2,:);w_mat_1(2,:);w_mat_1(2,:);w_mat_1(2,:)];

% save('w_mat_VAL_Kvary_forP.mat','w_mat')
% load('w_mat_VAL_Kvary_forP')
%%

% t_net_vect=linspace()
% t_net_vect=-0.08:0.01:2.6;
% x_net_vect=t_net_vect*vel;
% load('t_Santos_Filt')
load('t_INCT')

% x_net_vect_1=x_final_185.';

% x_net_vect=x_net_vect_1+3;
t_net_vect=t_INCT-3.727;

w_net_vect=zeros(1,length(t_net_vect));

tic;
for i=1:length(t_net_vect)
    t=t_net_vect(i);
    for j=1:N_W
        w_t_vect=w_mat(j,:);
        delta=delta_vect(j);
        w_net_vect(1,i)=w_net_vect(1,i)+ha_P(t_vect,w_t_vect,t-delta/vel);
    end
end
t_exec=toc;
w_net_mm=w_net_vect*1e3;
load('w_INCT');
%% FIG
figure;
ha_fig=tight_subplot(1,1,[.05 .13],[0.12 0.1],[0.11 0.07]);
hold on;
axes(ha_fig(1));
% x_o=x_net_vect-3;
% plot((t_Santos),fliplr(w_net_vect),'-k','LineWidth',1.5);

plot(t_INCT,w_net_mm,'-k','LineWidth',2);
% figure
hold on
% figure
plot(t_INCT,w_INCT,'-.b','LineWidth',2);

% grid on;
ylabel('$\bf w~(mm)$','Interpreter','latex','FontSize',17);
xlabel('$\bf x~(m)$','Interpreter','latex','FontSize',17);

set(ha_fig,'YTickLabelMode','auto');
set(ha_fig,'XTickLabelMode','auto');

set(ha_fig,'FontSize',15);
xlim([3.5 6.5])
ylim([-1.5 0.2])
%%
ha_col=@colors;
figure
ha_plot = tight_subplot(1, 1,[.12 .05],[0.2 0.07],[0.11 0.02]);
hold on;
set(ha_plot,'FontSize',10, 'Box', 'on','LineWidth',1,...
    'TickLabelInterpreter', 'latex','TickLength',[0.01, 0.01]);
set(ha_plot, 'Color', 'w')
axes(ha_plot(1));
hold on
plot(t_INCT,w_INCT,'color',ha_col('light blue'),'LineStyle','-','LineWidth',2);
hold on
plot(t_INCT,w_net_mm*1.0,':k','LineWidth',1.2);
ylabel('Deflection~(mm)','Interpreter','latex','FontSize',12);
xlabel('Distance~(m)','Interpreter','latex','FontSize',12);

set(ha_plot,'YTickLabelMode','auto');
set(ha_plot,'XTickLabelMode','auto');
legend({strcat('Santos et al (2016)'),strcat('Proposed Model')},...
    'Location','southeast',...
    'FontSize',10,'box','off','Interpreter','latex')
x0=3;
y0=1;
width=5;
height=2.5;
set(gcf,'Units','inches','position',[x0,y0,width,height])
xlim([3 7])
% ylim([-10 2])


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
cd 'save_fig'
print -dpdf -painters validation_SANTOS_201kmph
cd ..
% fileID=fopen('w_rail_multi_SANTOS_INCT_201kmph.txt','w');
% 
% % %
% file_1=[t_INCT;w_INCT;w_net_mm];
% fprintf(fileID,'%6s %16s %16s\n','t(s)','w_INCT','w_net_mm');
% fprintf(fileID,'%6.2f %15.10f %15.10f\n', file_1);


box on