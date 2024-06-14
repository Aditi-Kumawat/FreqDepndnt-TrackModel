% Just w(omega) for unit loading for Winkler Model

% sc_W_Kconst.m

function w_omega=get_w_omega_contour_Wink(omega_vect_1,vel,x_R,E_R,I_R,rho_R,c_R,K_eqt_0)
w_omega=((vel^3).*exp(-1i*(omega_vect_1).*(x_R/vel)))./(E_R*I_R*(omega_vect_1).^4+vel^4.*K_eqt_0-rho_R.*(omega_vect_1).^2.*vel^4+1i*(omega_vect_1)*c_R*vel^4);

end