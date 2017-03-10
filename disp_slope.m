function output1 = disp_slope(sampsize,Np,disp_r,tgrid,tmin,tstep,tmax)

figure('Visible','off');
slp = zeros(1,tmax/tstep);

sampr = randperm(Np,sampsize);

disp_r_sample = disp_r(sampr,:);
stddpos_r = std(disp_r_sample,0,1);

jj=0;
for ii = tmin:tstep:tmax
    jj = jj+1;
    startid = 2;%max(1,ii-2*floor(tmax/5));
    %[slp(jj), ~] = logfit(tgrid(startid:ii),stddpos_r(startid:ii).^2, 'loglog');
    logt = log10(tgrid(startid:ii));
    logvar = log10(stddpos_r(startid:ii).^2);
    p1 = robustfit(logt,logvar);
    slp(jj) = p1(2);
end

output1 = slp;

end
