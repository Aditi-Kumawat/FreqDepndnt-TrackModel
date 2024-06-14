%%% generating and saving contour plots %%%
tic
x_vect=omega_vect_1;
y_vect=v_R_vect;
Z=disp_3D(:,:,1);
% Z=abs(P_dyn_vect);
% Z=f_mat;
[X,Y]=meshgrid(x_vect,y_vect);
% level_vect=linspace(0,1,31);
level_vect=0:(1/31):1;

C=cell(1,1);h=cell(1,1);
figure;
ha_plot = fns_contour_plotting.tight_subplot(1,1,[.05 .03],[0.165 0.07],[0.13 0.07]);
hold on;
set(ha_plot,'FontSize',10, 'Box', 'on','LineWidth',1,...
    'TickLabelInterpreter', 'latex','TickLength',[0.01, 0.01]);
set(ha_plot, 'Color', 'w')
[C{1,1},h{1,1}]=contour(X,Y,Z);
ylim([y_vect(1) y_vect(end)]);
% set(ha(1),'YTick',(-1:0.25:1));
ylabel({'Velocity~Ratio,~$\alpha$'},'FontSize',11,'Interpreter','latex');
xlabel({'Frequency,~$\omega$~(rad/s)'},'FontSize',11,'Interpreter','latex');

set(h{1},'LineWidth',2);
set(h{1},'LevelList',level_vect);
set(h{1},'Fill','on');
% colormap default
colormap(hot)

hc=colorbar;
% col_pos=get(hc,'Position');
% set(hc,'Position',[col_pos(1)+0.22 col_pos(2) 0.7*col_pos(3) col_pos(4)]);
ylabel(hc,'$\widehat{w}_{N}(\omega,\alpha)$','FontSize',11,'Interpreter','latex')
set(hc,'FontSize',10,...
    'TickLabelInterpreter', 'latex');
hold on;
set(ha_plot,'XTickLabelMode','auto');
set(ha_plot,'YTickLabelMode','auto');
set(gcf,'Units','inches', 'Position', [3 3 4.5 3]);
grid off
set(gcf,'renderer','Painters')
saveas(gcf,'contour_Proposed_model.pdf')
% %%%%%%%%%%%%%%% Inserting right-hand-axis %%%%%%%%%%%%%%%%%%%%%
% ax1_pos = get(ha(1),'Position');
% ax2 = axes('Position',ax1_pos,'YAxisLocation','right','Color','none');
% ylim([-pi pi]);
% set(ax2,'YTick',[-3 0 3]);
% set(ax2,'YTickLabel',[-3 0 3]);
% set(ax2,'fontsize',18);
% ylab2=ylabel('another y-axis','FontSize',32);
% ylab2_pos=get(ylab2,'Position');
% set(ylab2, 'Units', 'Normalized', 'Position', [1.08, 0.5, 0]);
% set(ax2,'XTickLabel','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% title(strcat('Contour plot of FT of deflection'),'FontSize',24);
% set(gcf, 'Units','inches','Position', [4 4 4.5 3]);
% texec=toc;
% plot_filename_str=strcat('plot');
% saveas(gca,strcat(plot_filename_str,'.fig'));
% saveas(gca,strcat(plot_filename_str,'.bmp'));