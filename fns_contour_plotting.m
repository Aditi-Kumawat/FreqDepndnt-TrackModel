classdef fns_contour_plotting
    methods (Static)
        % This function file calculates the value of K_eqt_0 for the given input
        % parameters

        function [K_eqt_0,fn_phi]=get_K_eqt_0_0n(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho,c,K1_b,K_b,K_P,c_P)
            % omega_vect is in rad/s

            K_RP_vect=K_P+1i*c_P.*omega_0;

            fn_phi=zeros(1,length(omega_0));
            for k=1:length(omega_0)

                omega=omega_0(1,k);
                %     Disp=['omega = ',num2str(omega)];
                %     disp(Disp)
                rho_b=(rho*omega^2*L_S^4)/(E_S*I_S);   %non-dimensionalised inertial term
                c_b=(c*omega*L_S^4)/(E_S*I_S);         %non-dimensionalised damping term

                a1_b=K1_b;
                a2_b=(K_b-rho_b+1i*c_b);

                m1=sqrt((a1_b+sqrt(a1_b^2-4*a2_b))*0.5);
                m2=sqrt((a1_b-sqrt(a1_b^2-4*a2_b))*0.5);
                m3=-m1;
                m4=-m2;

                M=zeros(8,8);
                M(1,:)=[m1^2 m2^2 m3^2 m4^2 0 0 0 0];
                M(2,:)=[m1^3 m2^3 m3^3 m4^3 0 0 0 0];
                M(3,:)=[exp(a_b*m1) exp(a_b*m2) exp(a_b*m3) exp(a_b*m4) -exp(a_b*m1) -exp(a_b*m2) -exp(a_b*m3) -exp(a_b*m4)];
                M(4,:)=[m1*exp(a_b*m1) m2*exp(a_b*m2) m3*exp(a_b*m3) m4*exp(a_b*m4) -m1*exp(a_b*m1) -m2*exp(a_b*m2) -m3*exp(a_b*m3) -m4*exp(a_b*m4)];
                M(5,:)=[m1^3*exp(a_b*m1) m2^3*exp(a_b*m2) m3^3*exp(a_b*m3) m4^3*exp(a_b*m4) -m1^3*exp(a_b*m1) -m2^3*exp(a_b*m2) -m3^3*exp(a_b*m3) -m4^3*exp(a_b*m4)];
                M(6,:)=[m1^2*exp(a_b*m1) m2^2*exp(a_b*m2) m3^2*exp(a_b*m3) m4^2*exp(a_b*m4) -m1^2*exp(a_b*m1) -m2^2*exp(a_b*m2) -m3^2*exp(a_b*m3) -m4^2*exp(a_b*m4)];
                M(7,:)=[0 0 0 0 m1^3*exp(m1/2) m2^3*exp(m2/2) m3^3*exp(m3/2) m4^3*exp(m4/2)];
                M(8,:)=[0 0 0 0 m1*exp(m1/2) m2*exp(m2/2) m3*exp(m3/2) m4*exp(m4/2)];

                %Right-hand column matrix
                % % % % % % %     R=[0;0;0;0;(P*L^3/E/I);0;0;0]; %%%%%%%%%%
                phi_c_vect=[0;0;0;0;(L_S^3/E_S/I_S);0;0;0];  %phi_c_vect :column vector

                % A=[A1;A2;A3;A4;B1;B2;B3;B4];
                A=M\phi_c_vect;
                A1=A(1,1);A2=A(2,1);A3=A(3,1);A4=A(4,1);
                % B1=A(5,1);B2=A(6,1);B3=A(7,1);B4=A(8,1);
                fn_phi(1,k)=A1*exp(m1*x_b)+A2*exp(m2*x_b)+A3*exp(m3*x_b)+A4*exp(m4*x_b);   %x_b=a_b, I don't remember why I wrote it as x_b
            end

            K_sub_vect=-(1./fn_phi)./S;
            K_eqt_0=(K_RP_vect.*K_sub_vect)./(K_RP_vect+K_sub_vect);
            % K_eqt_vect=K_sub_vect;

        end

        %%
        % Just w(omega) for Proposed Model

        function Int=get_w_omega_contour(omega_vect_1,vel,x_R,E_R,I_R,rho_R,c_R,omega_vect,K_eqt_vect)
            ha_k=@fns_contour_plotting.get_Komega_vect_intrpol;
            K_omega_vect=ha_k(omega_vect,K_eqt_vect,omega_vect_1);
            Int=((vel^3).*exp(-1i*(omega_vect_1).*(x_R/vel)))./(E_R*I_R*(omega_vect_1).^4+vel^4.*K_omega_vect-rho_R.*(omega_vect_1).^2.*vel^4+1i*(omega_vect_1)*c_R*vel^4);
        end
        %%
        % This function file calculates the value of K(omega) at any given value of
        % omega through linear interpolation.
        function K_omega_vect=get_Komega_vect_intrpol(omega_vect,K_eqt_vect,omega_vect_1)

            o_lim_1=omega_vect(1,1);
            o_lim_2=omega_vect(1,end);
            o_0=omega_vect(1,2)-omega_vect(1,1);
            K_omega_vect=zeros(1,length(omega_vect_1));

            for j=1:length(omega_vect_1)
                omega=omega_vect_1(1,j);
                if omega>=o_lim_1 && omega<o_lim_2
                    i_num=((omega-o_lim_1)/o_0)+1;
                    i=floor(i_num);
                    K_omega_vect(1,j)=K_eqt_vect(i)+((K_eqt_vect(i+1)-K_eqt_vect(i))/(omega_vect(i+1)-omega_vect(i)))*(omega-omega_vect(i));
                else
                    K_omega_vect(1,j)=0;
                end
            end
        end
        %%
        function plot_contour(omega_vect_1,v_R_vect,disp_3D)
            %%% generating and saving contour plots %%%
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
        end
        %%
        function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

            % tight_subplot creates "subplot" axes with adjustable gaps and margins
            %
            % ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
            %
            %   in:  Nh      number of axes in hight (vertical direction)
            %        Nw      number of axes in width (horizontaldirection)
            %        gap     gaps between the axes in normalized units (0...1)
            %                   or [gap_h gap_w] for different gaps in height and width
            %        marg_h  margins in height in normalized units (0...1)
            %                   or [lower upper] for different lower and upper margins
            %        marg_w  margins in width in normalized units (0...1)
            %                   or [left right] for different left and right margins
            %
            %  out:  ha     array of handles of the axes objects
            %                   starting from upper left corner, going row-wise as in
            %                   going row-wise as in
            %
            %  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
            %           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
            %           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

            % Pekka Kumpulainen 20.6.2010   @tut.fi
            % Tampere University of Technology / Automation Science and Engineering


            if nargin<3; gap = .02; end
            if nargin<4 || isempty(marg_h); marg_h = .05; end
            if nargin<5; marg_w = .05; end

            if numel(gap)==1;
                gap = [gap gap];
            end
            if numel(marg_w)==1;
                marg_w = [marg_w marg_w];
            end
            if numel(marg_h)==1;
                marg_h = [marg_h marg_h];
            end

            axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh;
            axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

            py = 1-marg_h(2)-axh;

            ha = zeros(Nh*Nw,1);
            ii = 0;
            for ih = 1:Nh
                px = marg_w(1);

                for ix = 1:Nw
                    ii = ii+1;
                    ha(ii) = axes('Units','normalized', ...
                        'Position',[px py axw axh], ...
                        'XTickLabel','', ...
                        'YTickLabel','');
                    px = px+axw+gap(2);
                end
                py = py-axh-gap(1);
            end
        end
    end
end