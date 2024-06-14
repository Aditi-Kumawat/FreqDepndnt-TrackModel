% This script file provides the input parameters to other codes.
%%%%%%%%%%%%%%% RAILPROPERTIES %%%%%%%%%%%%%%%%%
E_R=210e9;                          % Young's modulus of beam material(N/m^2)
I_R=3055e-8;                     % MOI of beam C/S about neutral axis(m^4)
rho_R=60;                           % Density of rail beam in kg/m
x_R=0;
dr=0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_S= 2.52;                           %sleeper length(m)
B_S=0.3;                           %sleeper width(m)
H_S=0.233;                           %sleeper thickness(m)
S=0.67;                              %distance between two sleeper beams(m) 
G_L=1.5;                            %gauge length(m)
rho_S=300/(L_S);
% rho=(2400*B_S*H_S);
a=(L_S-G_L)/2;                      %distance of load from end(m)
a_b=a/L_S;                          %normalized rail beam distance from sleeper end
x_b=a_b;

E_S=38.45e9;                           %Young's modulus of beam material(N/m^2)
I_S=B_S*(H_S^3)/12;                     %MOI of beam C/S about neutral axis(m^4)
Bending_Stiffness_S=E_S*I_S;
G_B=62.5e6;                         %shear modulus of ballast(N/m^2)
G_SB=104.17e6;                    %shear modulus of sub-ballast(N/m^2)
% H_B=1.4;                          %height of ballast(m)
% H_SB=0;                          %height of ballast(m)
                        %height of sub-ballast(m)
H_B=0.3;
H_SB=0.25;
K=4.078e6;                          %subgrade modulus(N/m^2)/subgrade modulus(first parameter) of soil per unit beam length
% K1=666875;
K1=((G_B*H_B)+(G_SB*H_SB))*B_S;   %shear parameter of soil,(N)

% rho=(2500*B_S*H_S);                 %density in kg/m or(Ns^2/m^2)

w0_1=2*pi*196.202;                         %natural frequency of sleeper, rad/s(from natural_frequency.m)
w0_IF_sleeper=sqrt(K/rho_S);

f0=(1/2/pi)*w0_IF_sleeper;                         %natural frequency of sleeper, hz(from natural_frequency.m)
zeta=0.0004;                        %damping ratio
% zeta=0;
c=2*rho_S*w0_1*zeta;                    %coefficient of viscous damping per unit length

K1_b=(K1*L_S^2)/(E_S*I_S);          %non-dimensionalised shear parameter
K_b=(K*L_S^4)/(E_S*I_S);            %non-dimensionalised soil reaction parameter

K_P_star=5e8;                     %Rail-pad stiffness per unit length (N/m) (1984_Grassie and Cox_The Dynamic Response of Railway Track With Flexible Sleepers to High Frequency Vertical Excitation)
c_P_star=2.5e5;                      %Rail-pad viscous damping constant per unit length(Ns/m) (1984_Grassie and Cox_The Dynamic Response of Railway Track With Flexible Sleepers to High Frequency Vertical Excitation)
%
K_P=K_P_star/S;                     %Rail-pad stiffness per unit length (N/m^2)
c_P=c_P_star/S;                     %Rail-pad viscous damping constant per unit length(Ns/m^2)