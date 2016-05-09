function [output1] = cuml_pdf(npdfbins,Np,pos_r,time_r,times1,nposbins,maxipos,tgrid,tmin,tstep,tmax)
% short script to find slope of loglog decay of step size PDF

%maxipos = max(max(pos_r)); % in a real situation, this could be the maximum dimension of the region

yabins = maxipos/npdfbins:maxipos/npdfbins:maxipos;

posbins = [];
for ll = 1:nposbins
    posbins(:,ll) = 1+(ll-1)*Np/nposbins:ll*Np/nposbins;
end
slopes = [];
umma = [];

for ii = tmin:tstep:tmax
    checktime1 = tgrid(ii-1);
    checktime2 = tgrid(ii);
    fz_id_prev = 10;
    for kk = 1:Np
        
        for jj = 1:size(times1,2)
            checksum = tgrid(end) > times1(kk,jj:end);
%             scheck = sum(checksum) 
%             checktime2
%             t1 = times1(kk,jj)
%             checktime1
            if checktime2 > times1(kk,jj) && checktime1 < times1(kk,jj)
                umma = [umma pos_r(kk,jj)];
            elseif sum(checksum)==0
                fff = find(~pos_r(kk,:),1,'first');
                % the fraction of the final step before it goes off the
                % grid
                step_fraction = pos_r(kk,fff-1)*(max(tgrid)-checktime2)/time_r(kk,fff-1);
                umma = [umma step_fraction];
                break
            end
        end
    end
    figure('Visible','off');
    % find normalized histogram of step sizes
    nnz(umma);
    if sum(umma>0)
        [yim, yam] = histcounts(umma(:),yabins,'Normalization','pdf');
        % find the first zero in the histogram to avoid fitting the slope
        %to very sparse data
        firstzero_id = find(~yim(2:end),1,'first');
        % since the first zero in the histogram can oscillate, while the
        % PDF tends to elongate, do not exclude previously nonzero bins
        fz_id = max(firstzero_id,fz_id_prev);
        % make a note of the first nonzero bin for this iteration, to use
        % in the next iteration
        fz_id_prev = max(fz_id_prev,firstzero_id);
        % here, 1.5 is an adhoc scaling parameter that has worked in
        % validation, for keeping a significant portion of the tail while
        % excluding most zero elements that ruin the fit
        % Moreover, I have chosen the fourth bin as a rule of thumb for
        % the lower bound of the fitting range
        red_yam = yam(4:min(int8(1.5*fz_id),size(yam,2)));
        red_yim = yim(4:min(int8(1.5*fz_id),size(yam,2)));
        if nnz(red_yim) > 1
            [slp intcp] = logfit(red_yam, red_yim, 'loglog');
        else
            slp = 0;
        end
    else
        slp = 0;
    end
    %title([ii slp size(yam,2) ]); pause(0.2);
    close
    slopes = [slopes -1.*slp];
end

output1 = slopes;

end % of function