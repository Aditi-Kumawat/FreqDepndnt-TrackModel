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