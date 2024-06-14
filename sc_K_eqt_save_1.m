close all;clear;clc;
sc_IP_test2
omega_0=-1e3:0.02:1e3;
K_eqt_1=real(get_K_eqt_0_val1(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho_S,c,K1_b,K_b));
