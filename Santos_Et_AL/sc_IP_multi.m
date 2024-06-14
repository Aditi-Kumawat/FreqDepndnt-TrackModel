% This script file provides the input parameters to other codes..
%%%%%%%%%%%%%%% RAILPROPERTIES %%%%%%%%%%%%%%%%%
E_R=200e9;                          % Young's modulus of beam material(N/m^2)
I_R=3055e-8;
rho_R=60;                           % Density of rail beam in kg/m
x_R=0;
dr=0.3;
BS_R=E_R*I_R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_S= 2.52;                           %sleeper length(m)

B_S=0.3;                           %sleeper width(m)
H_S=0.22;                           %sleeper thickness(m)
% Used for calculating the multi wheel load profile at k_s=5.6MPa
S=0.6;                              %distance between two sleeper beams(m) (1998_Auersch L_Vehicle track interaction and Soil dynamics)
G_L=1.5;                            %gauge length(m)
% rho_S=300/(L_S);
rho_S1=1894; %in kg/m^3
rho_S=rho_S1*B_S*H_S; %in kg/m

% rho=(2400*B_S*H_S); 
a=(L_S-G_L)/2;                      %distance of load from end(m)
a_b=a/L_S;                          %normalized rail beam distance from sleeper end
x_b=a_b;

E_S=30e9;                           %Young's modulus of beam material(N/m^2)
I_S=B_S*(H_S^3)/12;                     %MOI of beam C/S about neutral axis(m^4)
Bending_Stiffness_S=E_S*I_S;
% G_B=62.5e6;                         %shear modulus of ballast(N/m^2)
% G_SB=104.17e6;                    %shear modulus of sub-ballast(N/m^2)
G_B=4.33e7;                         %shear modulus of ballast(N/m^2)
G_SB=8.15e7;                    %shear modulus of sub-ballast(N/m^2)
H_B=0.35;                          %height of ballast(m)
% H_B=0;                          %height of ballast(m)
H_SB=0.55;                        %height of sub-ballast(m)

K=96e6; %<THE VALUE>                          %subgrade modulus(N/m^2)/subgrade modulus(first parameter) of soil per unit beam length
% K=2.7e6;
K1=((G_B*H_B)+(G_SB*H_SB));   %shear parameter of soil,(N)

% rho=(2500*B_S*H_S);                 %density in kg/m or(Ns^2/m^2)

w0_1=2*pi*168.067;                         %natural frequency of sleeper, rad/s(from natural_frequency.m)
w0_IF_sleeper=sqrt(K/rho_S);

f0=(1/2/pi)*w0_IF_sleeper;                         %natural frequency of sleeper, hz(from natural_frequency.m)
zeta=0.0004;                        %damping ratio
% zeta=0; 
c=2*rho_S*w0_1*zeta;                    %coefficient of viscous damping per unit length

K1_b=(K1*L_S^2)/(E_S*I_S);          %non-dimensionalised shear parameter
K_b=(K*L_S^4)/(E_S*I_S);            %non-dimensionalised soil reaction parameter

K_P_star=6.2e8;                     %Rail-pad stiffness per unit length (N/m) (1984_Grassie and Cox_The Dynamic Response of Railway Track With Flexible Sleepers to High Frequency Vertical Excitation)
c_P_star=2.25e5;                      %Rail-pad viscous damping constant per unit length(Ns/m) (1984_Grassie and Cox_The Dynamic Response of Railway Track With Flexible Sleepers to High Frequency Vertical Excitation)

K_P=K_P_star/S;                     %Rail-pad stiffness per unit length (N/m^2)
c_P=c_P_star/S;                     %Rail-pad viscous damping constant per unit length(Ns/m^2)