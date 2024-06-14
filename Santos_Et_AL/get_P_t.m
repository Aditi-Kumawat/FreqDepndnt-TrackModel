% This function file finds the value of sleeper load (P_t) at any given value
% time (t) using linear interpolation.
% Uses the "P_t_vect" and "t_vect" generated for any given velocity from
% the code sc_P_t_vect.m.
function w_t=get_P_t(t_vect,w_t_vect,t_ny)

    t_lim_1=t_vect(1,1);
    t_lim_2=t_vect(1,end);
    d_t0=t_vect(1,2)-t_vect(1,1);
    
    if t_ny>=t_lim_1 && t_ny<t_lim_2
        i_num=((t_ny-t_lim_1)/d_t0)+1;
        i=floor(i_num);
        w_t= w_t_vect(i)+ ((w_t_vect(i+1)-w_t_vect(i))/(t_vect(i+1)-t_vect(i)))*(t_ny-t_vect(i));
    else
        w_t=0;
    end
    
end