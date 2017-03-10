function [moddelg, modelg2, modelpow, aX3cdflessfit, aX3cdfmorefit, likliplcdfless, likliplcdfmore, cpow] = fit_gravlevy(rho,tau1,tau2,alevy,EVDpop,EVDroaddist)

LR_params;

% acquire the nonzero distances and populations between those locations
% y should be count of occurrences of unique distances

% [uEVDrddmean, iEVDrddmean, iuEVDrddmean] = unique(nzEVDrdm);

onlyunique = 0;

% yEVDrdd = [];
% 
% for k = 1:size(uEVDrddmean) % start at k=2 to remove 0, which represents distances below resolution
%     yEVDrdd(k) = sum(uEVDrddmean(k)==nzEVDrdm); 
% end
% probyEVDrdd = yEVDrdd./sum(yEVDrdd);

% for X, you could select only one example of the unique distances
% in each of the vectors, but that is not representative of the data when
% asking for the probability of a case to go somewhere
% XEVDrdd = [EVDpop.mean(nzidEVDrdm(:),1),EVDpop.mean(nzidEVDrdm(:),2),nzEVDrdm(:)];
% here, I think I have to select the unique occurences to put into the
% nonlinear model fitting
if onlyunique==0
    if strcmp(chiefpop,'median')
        nzidEVDrdm = find(EVDroaddist.median>0);
        pop1 = EVDpop.median(nzidEVDrdm(:),1);
        pop2 = EVDpop.median(nzidEVDrdm(:),2);
        dist = EVDroaddist.median(nzidEVDrdm);
    elseif strcmp(chiefpop,'max')
        nzidEVDrdm = find(EVDroaddist.max>0);
        pop1 = EVDpop.max(nzidEVDrdm(:),1);
        pop2 = EVDpop.max(nzidEVDrdm(:),2);
        dist = EVDroaddist.max(nzidEVDrdm);
    elseif strcmp(chiefpop,'mean')
        nzidEVDrdm = find(EVDroaddist.mean>0);
        pop1 = EVDpop.mean(nzidEVDrdm(:),1);
        pop2 = EVDpop.mean(nzidEVDrdm(:),2);
        dist = EVDroaddist.mean(nzidEVDrdm);
    elseif strcmp(chiefpop,'min')
        nzidEVDrdm = find(EVDroaddist.min>0);
        pop1 = EVDpop.min(nzidEVDrdm(:),1);
        pop2 = EVDpop.min(nzidEVDrdm(:),2);
        dist = EVDroaddist.min(nzidEVDrdm);
    end
    XEVDrdd = [pop1,pop2,dist];
% elseif onlyunique==1
    % DONT USE THIS
%     pop1 = EVDpop.mean(nzidEVDrdm(iEVDrddmean),1);
%     pop2 = EVDpop.mean(nzidEVDrdm(iEVDrddmean),2);
%     dist = nzEVDrdm(iEVDrddmean);
%     XEVDrdd_Park = [pop1,pop2,dist];
end

% XEVDrddnorm = [EVDpop.mean(nzidEVDrdm(iEVDrddmean),1)./mean(EVDpop.mean(nzidEVDrdm(iEVDrddmean),1)),...
%                 EVDpop.mean(nzidEVDrdm(iEVDrddmean),2)./mean(EVDpop.mean(nzidEVDrdm(iEVDrddmean),2)),...
%                 nzEVDrdm(iEVDrddmean)./mean(nzEVDrdm(iEVDrddmean))];

% simple linear fit to the histogram of distances
linkdist_km = XEVDrdd(:,3)./1000;
%
nbinsEVD = 20;

[countX3,edgesX3,binX3] = histcounts(linkdist_km,nbinsEVD,'Normalization','count');
[probX3,edgesX3,binX3] = histcounts(linkdist_km,nbinsEVD,'Normalization','probability');

% start the tail at the maximum of the distribution
minbinEVD = find(countX3==max(countX3));
if minbinEVD == size(countX3,2)
    minbinEVD = floor(size(countX3,2)/2);
end

ecentEVD = (edgesX3(2:end)-edgesX3(1:end-1))./2+edgesX3(1:end-1);

l10centers = log10(ecentEVD);
l10counts = log10(countX3+1);

[cpow, Spow] = polyfit(l10centers(minbinEVD:end),l10counts(minbinEVD:end),1);
fittedX = linspace(min(l10centers),max(l10centers),100);
[fittedY,deltaY] = polyval(cpow,fittedX,Spow);

% figure; plot(l10centers,l10counts,'o')
% hold on; plot(fittedX,fittedY,'k-');
% title(['linear polyfit(), ',num2str(nbinsEVD),' bins, slope= ',num2str(coefpow(1))]);
% xlabel('log10(distance(km))');  ylabel('count');

%% binning to find generalized gravity law fit

%   TEST CODE
% rexp = 2.1;
% scale = 5;
% scale2 = 5;
% 
% lock_exponent = 1.8;
% 
% expt_gravity=[];
% expt_power=[];
% for jj = 1:size(probX3,2)
%     gig = find(binX3==jj);
%     expt_gravity(jj) = mean(scale.*XEVDrdd_Park(gig,1).*XEVDrdd_Park(gig,2)./XEVDrdd_Park(gig,3).^2);
%     expt_power(jj) = mean(scale2./XEVDrdd_Park(gig,3).^rexp);
% end

binmean_pop1 = zeros(size(probX3));
binmean_pop2 = zeros(size(probX3));
binmean_dist = zeros(size(probX3));

for jj = 1:size(probX3,2)
    gig = find(binX3==jj);
    if isempty(gig)
        0;
    else
        binmean_pop1(jj) = mean(XEVDrdd(gig,1));
        binmean_pop2(jj) = mean(XEVDrdd(gig,2));
        binmean_dist(jj) = mean(XEVDrdd(gig,3));
    end
end

logcounts = log10(countX3'+1);

opts = statset('Display','iter','TolFun',1e-8,'MaxIter',200);
logfitX = [log10(binmean_pop1+1); log10(binmean_pop2+1); log10(binmean_dist+1)];

modelfun_gravity = @(b,x)b(1) + x(:,1).*b(2) + x(:,2).*b(3) + x(:,3).*b(4);
betagrav = [1 1 1 -2];
moddelg = fitnlm(logfitX(:,minbinEVD:end)',logcounts(minbinEVD:end),modelfun_gravity,betagrav,'Options',opts);
moddelg.Coefficients;

modelfun_gravity2 = @(b,x)b(1) + x(:,1).*b(2) + x(:,2).*b(3) - rho.*x(:,3);
betagrav2 = [1 1 1];
modelg2 = fitnlm(logfitX(:,minbinEVD:end)',logcounts(minbinEVD:end),modelfun_gravity2,betagrav2,'Options',opts);
modelg2.Coefficients;

modelfun_power = @(b,x)b(1) + x(:,3).*b(2);
betapow = [1 -2];
modelpow = fitnlm(logfitX(:,minbinEVD:end)',logcounts(minbinEVD:end),modelfun_power,betapow,'Options',opts);
modelpow.Coefficients;

lesslink = find(linkdist_km<xcutvalue);
morelink = find(linkdist_km>xcutvalue);
midlink = find(linkdist_km>xminvalue & linkdist_km<xmaxvalue);
% 
% linkdist_less = [linkdist_less; linkdist_km(lesslink)];
% linkdist_more = [linkdist_more; linkdist_km(morelink)];
linkdist_less = linkdist_km(lesslink);
linkdist_more = linkdist_km(morelink);

if size(unique(linkdist_less),1)<2
    aX3cdfless = 0;
    likliplcdfless = 0;
else
    [aX3cdfless,bmX3cdfless,likliplcdfless] = plfit(linkdist_less,'xmin',xminvalue);
end
aX3cdflessfit = aX3cdfless;
% aX3cdfless = 0.5;
if size(unique(linkdist_more),1)<2
    aX3cdfmore = 0;
    likliplcdfmore = 0;
else
    [aX3cdfmore,bmX3cdfmore,likliplcdfmore] = plfit(linkdist_more,'xmin',xminvalue);
%     plplot(linkdist_more,xminvalue,aX3cdfmore); drawnow; pause(1);
end
aX3cdfmorefit = aX3cdfmore
% aX3cdfmore = 2.1;
% plplot(linkdist_more,xminvalue,aX3cdfmorefit); drawnow; pause(1); 
% Xgrav = [ones(size(binmean_pop1)); log10(binmean_pop1); log10(binmean_pop2); log10(binmean_dist)];
% logGrav = Xgrav';
% 
% [bgrav, bgravint, bgravr, bgravrin, bgravstats]  = regress(logP,logGrav);    % Removes NaN data
% 
% Xpow = [ones(size(binmean_pop1)); log10(binmean_dist)];
% bpow = regress(log10(ecentEVD'),Xpow');    % Removes NaN data
% 
% expt_gravity_norm = expt_gravity./sum(expt_gravity);
% expt_power_norm = expt_power./sum(expt_power);
% 
% % figure; plot(expt_gravity_norm,NX3,'.');
% % hold on; plot(expt_power_norm,NX3,'o');
% 
% gravitybinned = scale.*expt_gravity_norm;
% powerlawbinned = 1.*scale.*expt_power_norm;
% 
% figure; loglog(ecentEVD,probX3,'o','DisplayName','pdf data');
% hold on; plot(ecentEVD,gravitybinned,'o','DisplayName','gravity');
% hold on; plot(ecentEVD,powerlawbinned,'o','DisplayName','power law');
% xlabel('distance (m)'); ylabel('scaled probability density');
% 
% figure; plot(ecent,expt_gravity_norm./NX3,'o');
% hold on; plot(ecent,5.*expt_power_norm./NX3,'o');
end