% This function file calculates the magnitude of the static deflection (w_vect_st) and
% maximum static deflection (w_st_max) for the rail beam on winkler
% foundation(K_eqt_0=K(omega=0)).
function [w_st_max,w_vect_st]=get_static_max(x_vect,K_eqt_0,P_W,E_R,I_R)
w_vect_st=zeros(1,length(x_vect));
% BM_vect_st=zeros(1,length(x_vect));
K_b=K_eqt_0/E_R/I_R;

m_1=exp(pi*1i/4)*K_b^0.25;
m_2=exp(-pi*1i/4)*K_b^0.25;

M=zeros(2,2);
M(1,:)=[m_1 m_2];
M(2,:)=[m_1^3 m_2^3];

%Right-hand column matrix
R=[0;P_W/(2*E_R*I_R)];

% A=[A2;A4];
A=M\R;
A1=A(1,1);A2=A(2,1);

for k=1:length(x_vect)
    x_b=x_vect(1,k);
    if x_b<0
        w_vect_st(1,k)=A1*exp(m_1*x_b)+A2*exp(m_2*x_b);
%         BM_vect_st(1,k)=E_R*I_R*(m_1^2*A1*exp(m_1*x_b)+m_2^2*A2*exp(m_2*x_b));
    else
        w_vect_st(1,k)=A1*exp(-m_1*x_b)+A2*exp(-m_2*x_b);
%         BM_vect_st(1,k)=E_R*I_R*(m_1^2*A1*exp(-m_1*x_b)+m_2^2*A2*exp(-m_2*x_b));  
    end
end
w_st_max=max(abs(w_vect_st));
% BM_st_max=max(abs(BM_vect_st));
end