% This function file calculates the value of K_eqt_0 for the given input
% parameters

function [K_eqt_0,fn_phi]=get_K_eqt_0_val1(omega_0,L_S,S,a_b,x_b,E_S,I_S,rho,c,K1_b,K_b)
% omega_vect is in rad/s

% K_RP_vect=K_P+1i*c_P.*omega_0;

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

K_eqt_0=-(1./fn_phi)./S;
% K_eqt_0=(K_RP_vect.*K_sub_vect)./(K_RP_vect+K_sub_vect);
% K_eqt_vect=K_sub_vect;

end
