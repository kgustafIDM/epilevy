
tic
% number of trajectories, particles
Np = 1000;

mu_space = 2.1:0.1:2.6;
nu_space = 0.1:0.1:1.1;
detect_r = 1:1:100;

gamma = zeros(length(mu_space),length(nu_space));

% make the data noisier by decreasing sampling of Levy walk ensembles
sampsize = 1000;

tstep = 100; 
tstepmin = 10; 
tstepmaxm = 10000; 

nposbins = 1; 
npdfbins = 100;

pdf_slopes = [];
disp_slopes = [];
fp_binary = [];

for iiii = 1:length(mu_space)
    for jjjj = 1:length(nu_space)
        mu = mu_space(iiii);
        nu = nu_space(jjjj);
        muprime = mu - 1 + 1/nu;
        muprimenu = muprime*nu;
        smpn = muprimenu - 2*nu;
        if muprimenu >= 2
            if smpn >= 1
                gamma(iiii,jjjj) = 1;
            elseif smpn < 1 && smpn > 0
                gamma(iiii,jjjj) = 2 - smpn;
            end
        elseif muprimenu < 2 && muprimenu > 1
            if smpn >= 1
                gamma(iiii,jjjj) = muprimenu - 1;
            elseif smpn < 1 && smpn > 0
                gamma(iiii,jjjj) = 2*nu;
            end
        end
        if gamma(iiii,jjjj) > 0
            
            runname = sprintf('distr_dispLW_%.2f_%.2f',mu,nu);
            runname

            load([runname,'.mat'],'Np', 'disp_r', 'pos_r', 'tgrid', 'tgridsize', 'times1', 'time_r', 'max_step_size')
            tstepmax = min(tstepmaxm,tgridsize);
  
            slopes_t = cuml_pdf(sampsize,npdfbins,Np,pos_r,time_r,times1,nposbins,max_step_size,tgrid,tstepmin,tstep,tstepmax);
            d_slope = disp_slope(sampsize,Np,disp_r,tgrid,tstepmin,tstep,tstepmax);
            firstpass = first_pass(sampsize,Np,detect_r,disp_r,tgrid,tgridsize,tstepmin,tstep,tstepmax);
            firstpass_sum = sum(firstpass,1);
            
            clear Np disp_r pos_r tgrid tgridsize times1 time_r max_step_size;
        end
        
        pdf_slopes = [pdf_slopes; slopes_t];
        disp_slopes = [disp_slopes; d_slope];
        fp_binary = [fp_binary; firstpass_sum];
 
    end
end
toc