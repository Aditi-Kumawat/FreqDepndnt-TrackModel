% This script file calculates the deflection of rail beam in time domain (w_t_vect_4),
% at a given values of x_R at a given velocity(vel)/velocity ratio, as the wheel moves over the rail beam.
% and the max value of deflection (W_ABS_MAT_3n_Kcnst) at all the input
% values of v_R_vect..

% Input parameters: sc_IPn_3 --> K_S=1e7 N/m^2
% Rail beam UIC60
% K_eqt is constant= K_eqt_0 --> Winkler Model
function w_t_vect=get_wrail_forP1(P_W,vel,t_vect,K_eqt_0,rho_R,x_R,E_R,I_R,dr,omega_vect,K_eqt_vect)
v_cr_con_3n=(4*E_R*I_R*K_eqt_0/(rho_R^2))^0.25;
c_cr_con=2*sqrt(K_eqt_0*rho_R);
c_R=dr*c_cr_con;
v_R=vel/v_cr_con_3n;
x_vect=t_vect*vel;
ha_st=@get_static_max;
[w_st_max_val1,w_vect_st]=ha_st(x_vect,K_eqt_0,P_W,E_R,I_R);

lim_1=-1e4;
lim_2=1e4;

o_I=0.0001;
ha_1=@get_w_rail_omega;

if v_R>=0 && v_R<=0.000800160032006401
    w_t_vect=w_vect_st;
else    
    In_t_vect=real(integral(@(omega_vect_1)ha_1(omega_vect_1,vel,x_R,t_vect,E_R,I_R,rho_R,c_R,o_I,omega_vect,K_eqt_vect),lim_1,lim_2,'ArrayValued',true));
    
    In_t_vect_1=In_t_vect.*exp(o_I*t_vect);
    
    w_t_vect=fliplr(-(P_W/pi/2).*In_t_vect_1);    % Deflection when K is function of omega + No Damping
end




