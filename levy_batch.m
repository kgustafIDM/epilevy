
% number of timesteps to initialize the matrix of steps
Nt = 100;
% number of trajectories, particles
Np = 1000;
% the final time in the units of interest (days for example)
final_t  = 3650;

%generalized velocity, units are just dist/time^nu
gen_vel = 2;

% maximum single step size, in km for example, set by the boundaries of the
% region where the random walkers are allowed
max_step_size = 1000;

mu_space1 = 2:0.1:2.5;
nu_space1 = 1.1:0.1:2;

mu_space2 = 1.5:0.1:3;
nu_space2 = 0.1:0.1:0.4;

gamma = zeros(length(mu_space1),length(nu_space1));

for iii = 1:length(mu_space1)
    for jjj = 1:length(nu_space1)
        mu = mu_space1(iii)
        nu = nu_space1(jjj)
        muprime = mu - 1 + 1/nu;
        muprimenu = muprime*nu;
        smpn = muprimenu - 2*nu;
        if muprimenu >= 2
            if smpn >= 1
                gamma(iii,jjj) = 1;
            elseif smpn < 1 && smpn > 0
                gamma(iii,jjj) = 2 - smpn;
            end
        elseif muprimenu < 2 && muprimenu > 1
            if smpn >= 1
                gamma(iii,jjj) = muprimenu - 1;
            elseif smpn < 1 && smpn > 0
                gamma(iii,jjj) = 2*nu;
            end
        end
        if gamma(iii,jjj) > 0
            gammarun = gamma(iii,jjj);
            kpareto = 1/(mu-1);
            % this "threshold parameter" essentially sets the scale of the random walk
            % through the coupling to sigma that I introduce to eliminate the 1 + term
            % in the parentheses of the Pareto distribution
            thpareto = 1;
            % this coupling removes the 1 + term in the Pareto distribution
            sigpareto = kpareto*thpareto;
            % the time-scale exponent for the Lévy walk coupling
            [pos_x,pos_y,time_r] = levywalk_gen(mu,kpareto,thpareto,sigpareto,nu,Nt,Np,final_t,max_step_size,gen_vel);
            
            % process the sequence of random numbers produced in the Lévy walk fxn
            dispLW_distr
            runname = sprintf('distr_dispLW_%.2f_%.2f',mu,nu);

            save([runname,'.mat']);

        end
    end
end

gamma = zeros(length(mu_space2),length(nu_space2));

for iii = 13:length(mu_space2)
    for jjj = 1:length(nu_space2)
        mu = mu_space2(iii)
        nu = nu_space2(jjj)
        muprime = mu - 1 + 1/nu;
        muprimenu = muprime*nu;
        smpn = muprimenu - 2*nu;
        if muprimenu >= 2
            if smpn >= 1
                gamma(iii,jjj) = 1;
            elseif smpn < 1 && smpn > 0
                gamma(iii,jjj) = 2 - smpn;
            end
        elseif muprimenu < 2 && muprimenu > 1
            if smpn >= 1
                gamma(iii,jjj) = muprimenu - 1;
            elseif smpn < 1 && smpn > 0
                gamma(iii,jjj) = 2*nu;
            end
        end
        if gamma(iii,jjj) > 0
            gammarun = gamma(iii,jjj);
            kpareto = 1/(mu-1);
            % this "threshold parameter" essentially sets the scale of the random walk
            % through the coupling to sigma that I introduce to eliminate the 1 + term
            % in the parentheses of the Pareto distribution
            thpareto = 1;
            % this coupling removes the 1 + term in the Pareto distribution
            sigpareto = kpareto*thpareto;
            % the time-scale exponent for the Lévy walk coupling
            [pos_x,pos_y,time_r] = levywalk_gen(mu,kpareto,thpareto,sigpareto,nu,Nt,Np,final_t,max_step_size,gen_vel);
            
            % process the sequence of random numbers produced in the Lévy walk fxn
            dispLW_distr
            runname = sprintf('distr_dispLW_%.2f_%.2f',mu,nu);

            save([runname,'.mat']);

        end
    end
end