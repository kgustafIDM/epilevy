function [output1] = cuml_pdf(sampsize,npdfbins,Np,pos_r,time_r,times1,nposbins,maxipos_give,tgrid,tmin,tstep,tmax)
% short script to find slope of loglog decay of step size PDF

% maxipos_emp = max(max(pos_r)); % in a real situation, this could be the maximum dimension of the region
% maxipos = min(maxipos_emp,maxipos_give);
% yabins = maxipos/npdfbins:maxipos/npdfbins:maxipos;

posbins = [];
for ll = 1:nposbins
    posbins(:,ll) = 1+(ll-1)*Np/nposbins:ll*Np/nposbins;
end
slopes = [];
umma = []; % initialized outside the loop purposefully - it is a cumulative tracker of step sizes

sampr = randperm(Np,sampsize);
fz_id_def = 20;

% initialize the left side of the timeframe
checktime1 = tgrid(1);

for ii = tmin:tstep:tmax
    checktime2 = tgrid(ii);
    for kk = sampr
        for jj = 1:size(times1,2)
            % check whether this time interval stretches beyond the
            % timeframe
            checksum = checktime2 < times1(kk,end);
            if checktime1 <= times1(kk,jj) && checktime2 >= times1(kk,jj)
                umma = [umma pos_r(kk,jj)];
            elseif sum(checksum)==0
                fff = find(~pos_r(kk,:),1,'first');
                % the fraction of the final step before it goes off the
                % grid is computed by linear interpolation
                step_fraction = pos_r(kk,fff-1)*(max(tgrid)-checktime2)/time_r(kk,fff-1);
                umma = [umma step_fraction];
                break
            end
        end
    end
    
    %figure('Visible','off');
    % find normalized histogram of step sizes
    %     nnz(umma);
    %
    %     maxipos_emp = max(max(umma)); % in a real situation, this could be the maximum dimension of the region
    %     maxipos = min(maxipos_emp,maxipos_give);
    %     yabins = min(umma(:)):maxipos/npdfbins:maxipos;
    
    if sum(umma>0)
        % using plfit seems to eliminate the need for histogram smoothing
        % sometimes, but other times it blows up
%         [alpha, xmin, L] = plfit(umma);
%         slp = alpha
                [yim, yam] = histcounts(umma(:),100,'Normalization','pdf');
        % find the first zero in the histogram to avoid fitting the slope
        %to very sparse data
        yimrun = runmean(yim,5);
        %         firstzero_id = find(~yimrun,1,'first');
        %         if isempty(firstzero_id)
        %             firstzero_id = 1;
        %         end
        % since the first zero in the histogram can oscillate, while the
        % PDF tends to elongate, do not exclude previously nonzero bins
        %fz_id = max(firstzero_id,fz_id_def);
        %peak_size = max(yimrun(fz_id:end));
%         max_id = 1;%fz_id + find(yimrun(fz_id:end)==peak_size,1);
%         tail_zero = find(~yim,1);
%         if isempty(tail_zero)
%             tail_zeroid = size(yim,2);
%         else
%             tail_zeroid = find(~yim,1);
%         end
        % make a note of the first nonzero bin for this iteration, to use
        % in the next iteration
        % I have chosen the peak value, indexed by max_id, the
        % start the fitting of the tail.
        red_yam = yam(10:90);
        red_yimrun = yimrun(10:90);
        if nnz(red_yimrun) > 1
            logy1 = log10(red_yam);
            logy2 = log10(red_yimrun);
            p1 = robustfit(logy1,logy2);
            %slp = -1*alpha;
            slp = -1*p1(2);
        else
            slp = 0;
        end
    else
        slp = 0;
    end
    %title([ii slp nnz(umma) ]); pause(0);
    %hold on; loglog(yam(2:end),yimrun); loglog(yam(2:end),yim);
    slopes = [slopes slp];
    checktime1 = checktime2;
end

output1 = slopes;

end % of function