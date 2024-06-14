% This function file calculates the integrand (int_4) for the inverse fourier
% transform of the deflection function i.e rail deflection as a function of omega
% (W(omega)), when K is function of
% omega and K(omega) has been calculated using linear interpolation..
function int_4=get_w_rail_omega(omega_vect_1,vel,x_R,t_vect,E_R,I_R,rho_R,c_R,o_I,omega_vect,K_eqt_vect)
ha_k=@get_Komega_vect_intrpol;
K_omega_vect=ha_k(omega_vect,K_eqt_vect,omega_vect_1);
% I_1=((vel^3).*exp(-1i*(omega_vect-1i*o_I).*(x_R/vel)))./(E_R*I_R*(omega_vect-1i*o_I).^4+vel^4.*K_omega_vect-rho_R.*(omega_vect-1i*o_I).^2.*vel^4);
I_1=((vel^3).*exp(-1i*(omega_vect_1-1i*o_I).*(x_R/vel)))./(E_R*I_R*(omega_vect_1-1i*o_I).^4+vel^4.*K_omega_vect-rho_R.*(omega_vect_1-1i*o_I).^2.*vel^4+1i*(omega_vect_1-1i*o_I)*c_R*vel^4);
int_4=I_1.*exp(1i*omega_vect_1*t_vect(1:end));
end