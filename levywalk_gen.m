function [out1, out2, out3] = levywalk_gen(mu,kpareto,thpareto,sigpareto,nu,Nt,Np,final_t,max_step_size,gen_vel)

% generate a Levy walk with stepsize probability r^-mu and t = r^-nu
% inner loop on each particle 1) step through a time grid with step size
% delta_t, 2) find the current time cur_t=i*delta_t, 3) while the Levy step
% time lev_t (initialized to zero) is less than cur_t 4) generate step_r and
% step_t, 5) set lev_t = lev_t + step_t, 6) init_r = final_r, 7) final_r =
% init_r + step_r, 8) end while, 9) linear interpolate between init_r and
% final_r to find step size during each delta_t

pos_x = zeros(Np,Nt);
pos_y = zeros(Np,Nt);
time_r = zeros(Np,Nt);

% time loop from 0 to final_t

for j_p = 1:Np
    if mod(j_p,Np/100)==0
        j_p
    end
    
    step_t = 0;
    cur_t = 0;
    i_t = 0;
        
    while cur_t <= final_t
        
        step_x = (random('gp',kpareto,sigpareto,thpareto,1,1)).*sign(random('uniform',-1,1,1));
        if abs(step_x) > max_step_size
            step_x = max_step_size;
        end
        
        step_y = (random('gp',kpareto,sigpareto,thpareto,1,1));
        if abs(step_y) > max_step_size
            step_y = max_step_size;
        end
        
        step_t = gen_vel.*step_y^(1/nu);
        
        cur_t = cur_t + step_t;
        i_t = i_t + 1;
        
        pos_x(j_p,i_t) = step_x;
        pos_y(j_p,i_t) = step_y;
        time_r(j_p,i_t) = step_t;
        
    end
end

out1 = pos_x;
out2 = pos_y;
out3 = time_r;

end