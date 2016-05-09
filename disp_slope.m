function output1 = disp_slope(sampsize,disp_r,tgrid,tmin,tstep,tmax)

figure('Visible','off');
slp = zeros(1,tmax/tstep);

sampr = randi(size(disp_r,1),1,sampsize);

disp_r_sample = disp_r(sampr,:);
stddpos_r = std(disp_r_sample,0,1);

jj=0;
for ii = tmin:tstep:tmax
    jj = jj+1;
    [slp(jj) intcp] = logfit(tgrid(1:ii),stddpos_r(1:ii).^2, 'loglog');
end

output1 = slp;

end
