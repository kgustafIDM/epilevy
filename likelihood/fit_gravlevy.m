%% acquire the nonzero distances and populations between those locations
% y should be count of occurrences of unique distances

nzidEVDrdm = find(EVDroaddist.mean>0);
nzEVDrdm = EVDroaddist.mean(nzidEVDrdm);
[uEVDrddmean, iEVDrddmean, iuEVDrddmean] = unique(nzEVDrdm);

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
    pop1 = EVDpop.median(nzidEVDrdm(:),1);
    pop2 = EVDpop.median(nzidEVDrdm(:),2);
    dist = nzEVDrdm(:);
    XEVDrdd_Park = [pop1,pop2,dist];
elseif onlyunique==1
    pop1 = EVDpop.mean(nzidEVDrdm(iEVDrddmean),1);
    pop2 = EVDpop.mean(nzidEVDrdm(iEVDrddmean),2);
    dist = nzEVDrdm(iEVDrddmean);
    XEVDrdd_Park = [pop1,pop2,dist];
end

% XEVDrddnorm = [EVDpop.mean(nzidEVDrdm(iEVDrddmean),1)./mean(EVDpop.mean(nzidEVDrdm(iEVDrddmean),1)),...
%                 EVDpop.mean(nzidEVDrdm(iEVDrddmean),2)./mean(EVDpop.mean(nzidEVDrdm(iEVDrddmean),2)),...
%                 nzEVDrdm(iEVDrddmean)./mean(nzEVDrdm(iEVDrddmean))];

% simple linear fit to the histogram of distances
linkdist_km_Park = XEVDrdd_Park(:,3)./1000;
%%
nbinsEVD = 20;

[countX3,edgesX3,binX3] = histcounts(linkdist_km_Park,nbinsEVD,'Normalization','count');
[probX3,edgesX3,binX3] = histcounts(linkdist_km_Park,nbinsEVD,'Normalization','probability');

% start the tail at the maximum of the distribution
minbinEVD = find(countX3==max(countX3));
if minbinEVD == size(countX3,2)
    minbinEVD = floor(size(countX3,2)/2);
end

ecentEVD = (edgesX3(2:end)-edgesX3(1:end-1))./2+edgesX3(1:end-1);

l10centers = log10(ecentEVD);
l10counts = log10(countX3+1);

[coefpow, Spow] = polyfit(l10centers(minbinEVD:end),l10counts(minbinEVD:end),1);
fittedX = linspace(min(l10centers),max(l10centers),100);
[fittedY,deltaY] = polyval(coefpow,fittedX,Spow);

% figure; plot(l10centers,l10counts,'o')
% hold on; plot(fittedX,fittedY,'k-');

%% binning to find generalized gravity law fit

rexp = 2.1;
scale = 5;
scale2 = 5;

lock_exponent = 1.8;

expt_gravity=[];
expt_power=[];
for jj = 1:size(probX3,2)
    gig = find(binX3==jj);
    expt_gravity(jj) = mean(scale.*XEVDrdd_Park(gig,1).*XEVDrdd_Park(gig,2)./XEVDrdd_Park(gig,3).^2);
    expt_power(jj) = mean(scale2./XEVDrdd_Park(gig,3).^rexp);
end

binmean_pop1 = zeros(size(probX3));
binmean_pop2 = zeros(size(probX3));
binmean_dist = zeros(size(probX3));

for jj = 1:size(probX3,2)
    gig = find(binX3==jj);
    if isempty(gig)
        0;
    else
        binmean_pop1(jj) = mean(XEVDrdd_Park(gig,1));
        binmean_pop2(jj) = mean(XEVDrdd_Park(gig,2));
        binmean_dist(jj) = mean(XEVDrdd_Park(gig,3));
    end
end

logcounts = log10(countX3'+1);

opts = statset('Display','iter','TolFun',1e-8,'MaxIter',200);
logfitX = [log10(binmean_pop1+1); log10(binmean_pop2+1); log10(binmean_dist+1)];

modelfun_gravity = @(b,x)b(1) + x(:,1).*b(2) + x(:,2).*b(3) + x(:,3).*b(4);
betagrav = [1 1 1 -2];
mdlg = fitnlm(logfitX(:,minbinEVD:end)',logcounts(minbinEVD:end),modelfun_gravity,betagrav,'Options',opts);
mdlg.Coefficients;

modelfun_gravity2 = @(b,x)b(1) + x(:,1).*b(2) + x(:,2).*b(3) - gravgamma0.*x(:,3);
betagrav2 = [1 1 1];
mdlg2 = fitnlm(logfitX(:,minbinEVD:end)',logcounts(minbinEVD:end),modelfun_gravity2,betagrav2,'Options',opts);
mdlg2.Coefficients;

modelfun_power = @(b,x)b(1) + x(:,3).*b(2);
betapow = [1 -2];
mdlpow = fitnlm(logfitX(:,minbinEVD:end)',logcounts(minbinEVD:end),modelfun_power,betapow,'Options',opts);
mdlpow.Coefficients;

xminvalue = 50;
xmaxvalue = 500;

lesslink_Park = find(linkdist_km_Park<xmaxvalue);

if size(unique(linkdist_km_Park(lesslink_Park)),1)<2
    aX3cdf = 0;
    likliplcdf = 0;
else
    [aX3cdf,bmX3cdf,likliplcdf] = plfit(linkdist_km_Park(lesslink_Park),'xmin',xminvalue);
end

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
