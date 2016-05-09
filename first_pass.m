function output = first_pass(detect_r,Np,disp_r,tgrid,tgridsize,tmin,tstep,tmax)
% short script to find first-passage time

fp_time = zeros(length(detect_r),Np);

fp_matrix = [];

for kk = 1:length(detect_r)
    reduce_disp_r = disp_r(:,tmin:tstep:tmax);
    fp_r = zeros(Np,tgridsize)+1;
    [xa xb] = find(reduce_disp_r >= detect_r(kk));
    if size(xb>0)
        fp_b = [zeros(1,min(xb)-1) ones(1,size(reduce_disp_r,2)-min(xb)+1)];
    else
        fp_b = zeros(1,size(reduce_disp_r,2));
    end
    
    for jj = 1:Np
        for ii = 1:tgridsize
            if disp_r(jj,ii) >= detect_r(kk)
                fp_r(jj,ii) = 0;
            end
        end
    end
    
    firstpass = zeros(Np,1);
    for jj = 1:Np
        fff = find(~fp_r(jj,:),1,'first');
        if isempty(fff)
            firstpass(jj) = tgrid(tgridsize);
        else
            firstpass(jj) = tgrid(fff);
        end
    end
    
    fp_time(kk,:) = firstpass;
    fp_matrix = [fp_matrix; fp_b];
    
end

% 
% mean_fpassage_bins = [];
% for ii = 1:size(walk_bins,1)
%     mean_fpassage_bins(ii) = mean(fpassage(walkbins(ii,:)));
% end

output = fp_matrix;

%%
% pass_r = [];
% 
% for jj = 1:size(dpos_x,2)
%     jj
%     ii
%     ii = 1;
%     passage = 0;
%     while passage == 0
%         ii = ii + 1;
%         if disp_r(ii,jj) >= detect_r
%             pass_r(jj) = time(ii);
%             passage = 1;
%         elseif ii >= length(time)
%             passage = 1;
%         end
%     end
% end

end % of function