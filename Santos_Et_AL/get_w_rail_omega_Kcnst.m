% This function file calculates the integrand for the inverse fourier
% transform of the deflection function (W(omega)), when K is constant equal
% to K_eqt_0..

% sc_rail_3n_Kcnst.m
function int_4=get_w_rail_omega_Kcnst(omega_vect_1,vel,x_R,t_vect,E_R,I_R,rho_R,c_R,o_I,K_eqt_0)
I_1=((vel^3).*exp(-1i*(omega_vect_1-1i*o_I).*(x_R/vel)))./...
    (E_R*I_R*(omega_vect_1-1i*o_I).^4+vel^4.*K_eqt_0-rho_R.*...
    (omega_vect_1-1i*o_I).^2.*vel^4+1i*(omega_vect_1-1i*o_I)*c_R*vel^4);
int_4=I_1.*exp(1i*omega_vect_1*t_vect(1:end));
end